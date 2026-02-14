#!/bin/bash
set -o nounset
set -e # Exit immediately if a command exits with a non-zero status

# =============================================================================
# Configuration
# =============================================================================
# Input files and parameters based on user request:
# RAGA.sh -r Col-CEN_v1.2.fasta -c CRR591673.fastq -t 80
REF="Col-CEN_v1.2.fasta"
CCS="CRR591673.fastq"
THREADS="80"

# Default parameters from RAGA.sh / RAGA-same.sh
NPR="3"         # Number of Polishing Rounds
DFI="99"        # Minimum alignment identity [0, 100]
DFL="20000"     # Minimum alignment length
PER="0.9"       # Target PacBio HiFi read align length >= *% of its own length
PER_ALT="0.5"   # Target longAlt read align length >= *% of its own length

# Derived filenames
REF_BASE=$(basename "$REF")
CCS_BASE=$(basename "$CCS")
REF_NAME="${REF_BASE%.*}"
CCS_NAME="${CCS_BASE%.*}"

echo "Starting RAGA Workflow..."
echo "Reference: $REF"
echo "CCS Reads: $CCS"
echo "Threads: $THREADS"
echo "-----------------------------------------------------------------------------"

# =============================================================================
# Step 1: Initial Assembly
# =============================================================================
# Purpose: Generate an initial de novo assembly of the target genome using HiFi reads.
# Software: hifiasm, awk

echo "[Step 1] Initial Assembly"
mkdir -p Initial_assembly
cd Initial_assembly

# Link CCS reads to current directory
ln -sf "../$CCS" "$CCS_BASE"

# Run hifiasm for de novo assembly (homozygous mode -l0)
# -l0: Disables duplication purging for homozygous genomes
hifiasm -o contigs -l0 "$CCS_BASE" -t "$THREADS"

# Convert GFA output to FASTA
awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
ln -sf contigs.bp.p_ctg.fa initial.fa

# Clean up link
rm "$CCS_BASE"
cd ..

INITIAL_ASM="Initial_assembly/initial.fa"
INITIAL_ASM_BASE=$(basename "$INITIAL_ASM")
INITIAL_ASM_NAME="${INITIAL_ASM_BASE%.*}"

echo "[Step 1] Done. Initial assembly at $INITIAL_ASM"
echo "-----------------------------------------------------------------------------"

# =============================================================================
# Step 2: Generate Alternative Long Reads
# =============================================================================
# Purpose: Use the reference genome to guide the assembly and generate "longAlt" reads
#          for gap filling and improvement.
# Software: minimap2, seqkit, racon, ragtag.py, nucmer, delta-filter, show-coords, hifiasm, samtools

echo "[Step 2] Generate Alternative Long Reads"
OUT_DIR="Alternative_reads"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# Link input files
ln -sf "../$REF" "${REF_NAME}_racon0.fa"
ln -sf "../$INITIAL_ASM" "$INITIAL_ASM_BASE"
ln -sf "../$CCS" "$CCS_BASE"

# 2.1 Filter CCS reads relevant to the target assembly
# Software: minimap2 (alignment), awk (filtering), seqkit (extraction)
echo "  - Filtering CCS reads..."
minimap2 -x map-hifi "../$INITIAL_ASM" "../$CCS" -t "$THREADS" > "${INITIAL_ASM_NAME}_${CCS_NAME}.paf"
# Keep reads with >= 90% alignment coverage
awk '$10/$2>=0.90{print $1}' "${INITIAL_ASM_NAME}_${CCS_NAME}.paf" | sort | uniq > "${CCS_NAME}_f.id"
seqkit grep -f "${CCS_NAME}_f.id" "../$CCS" > "${CCS_NAME}_f.fq"

# 2.2 Polish the Reference Genome using CCS reads
# Software: minimap2 (alignment), racon (polishing)
echo "  - Polishing reference genome ($NPR rounds)..."
for i in $(seq 1 "$NPR"); do
    prev=$((i-1))
    # Align reads to current reference iteration
    minimap2 -ax map-hifi "${REF_NAME}_racon${prev}.fa" "${CCS_NAME}_f.fq" -t "$THREADS" > "${REF_NAME}_racon${prev}.sam"

    # Polish using Racon
    if [ "$prev" == "0" ]; then
         # First round uses the symlinked original reference
         racon "${CCS_NAME}_f.fq" "${REF_NAME}_racon${prev}.sam" "../${REF_NAME}_racon${prev}.fa" -t "$THREADS" > "${REF_NAME}_racon${i}.fa"
    else
         racon "${CCS_NAME}_f.fq" "${REF_NAME}_racon${prev}.sam" "${REF_NAME}_racon${prev}.fa" -t "$THREADS" > "${REF_NAME}_racon${i}.fa"
    fi
    echo "    Round $i done."
done

# Get length of polished reference
seqkit fx2tab -l -n "${REF_NAME}_racon${NPR}.fa" | awk '{print $1"\t0\t"$5}' > "${REF_NAME}_racon${NPR}_len.txt"

# 2.3 Scaffolding & Gap Location
# Software: ragtag.py (homology-based scaffolding)
echo "  - Scaffolding target assembly against polished reference..."
ragtag.py scaffold "${REF_NAME}_racon${NPR}.fa" "../$INITIAL_ASM" -t "$THREADS" -o .
ln -sf ragtag.scaffold.fasta "${INITIAL_ASM_NAME}_ragtag.fa"

# Identify gaps (U = unplaced/gap)
awk '$5=="U"{if($2-1000000>=0){print $1"\t"$2-1000000"\t"$3+1000000} else {print $1"\t0\t"$3+1000000}}' ragtag.scaffold.agp > gap_around.bed

# 2.4 Locate expanded gap area on Reference
# Software: nucmer, delta-filter, show-coords (MUMmer suite)
echo "  - Locating gap areas on reference..."
nucmer -p "${REF_NAME}_racon${NPR}To${INITIAL_ASM_NAME}_ragtag" "${REF_NAME}_racon${NPR}.fa" "${INITIAL_ASM_NAME}_ragtag.fa" -t "$THREADS"
delta-filter -i "$DFI" -l "$DFL" "${REF_NAME}_racon${NPR}To${INITIAL_ASM_NAME}_ragtag.delta" > "${REF_NAME}_racon${NPR}To${INITIAL_ASM_NAME}_ragtag_f.delta"
show-coords -b -B "${REF_NAME}_racon${NPR}To${INITIAL_ASM_NAME}_ragtag_f.delta" | awk '{if($9<$10){print $1"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} else{print $1"\t"$8"\t"$10"\t"$9"\t"$11"\t"$12}}' > "${REF_NAME}_racon${NPR}To${INITIAL_ASM_NAME}_ragtag_f.coords"

# Intersect gaps with alignments to find coordinates on reference
awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($1 in a && $2"_RagTag" in a && $4>=b[$1] && $3<=c[$1]){print $2"\t"$5"\t"$6}}' gap_around.bed "${REF_NAME}_racon${NPR}To${INITIAL_ASM_NAME}_ragtag_f.coords" > "${REF_NAME}_racon${NPR}_gapA.bed"

# Extract "longAlt" reference sequences for these gaps
mkdir -p longAlt_ref
flag=0
while read id start end; do
    flag=$((flag + 1))
    echo -e "$id\t$start\t$end" > "longAlt_ref/gapA_${flag}.bed"
    seqkit subseq --bed "longAlt_ref/gapA_${flag}.bed" "${REF_NAME}_racon${NPR}.fa" -j "$THREADS" | seqkit replace -p ":." -j "$THREADS" > "longAlt_ref/gapA_${flag}.fa"
    rm "longAlt_ref/gapA_${flag}.bed"
done < "${REF_NAME}_racon${NPR}_gapA.bed"

# Visualization (Optional but part of workflow)
pl_locate.pl -f "${REF_NAME}_racon${NPR}_gapA.bed" -l "${REF_NAME}_racon${NPR}_len.txt" -o "${REF_NAME}_racon${NPR}_gapA.svg" || echo "Warning: SVG generation failed"

# 2.5 Extract CCS reads for Gap Regions
# Software: minimap2, awk, seqkit
echo "  - Extracting CCS reads for gap regions..."
minimap2 -x map-hifi "${REF_NAME}_racon${NPR}.fa" "../$CCS" -t "$THREADS" > "${REF_NAME}_racon${NPR}To${CCS_NAME}.paf"

mkdir -p ccsAlt_tgt
flag=0
while read id start end; do
    flag=$((flag + 1))
    echo -e "$id\t$start\t$end" > "ccsAlt_tgt/gapA_${flag}.bed"
    # Filter reads mapping to this region with high coverage/quality
    awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($6 in a && $9>=b[$6] && $8<=c[$6] && $10/$2>="'$PER'" && $13~/tp:A:P/){print $1}}' "ccsAlt_tgt/gapA_${flag}.bed" "${REF_NAME}_racon${NPR}To${CCS_NAME}.paf" | sort | uniq > "ccsAlt_tgt/gapA_${flag}.id"
    seqkit grep -f "ccsAlt_tgt/gapA_${flag}.id" "../$CCS" -j "$THREADS" > "ccsAlt_tgt/gapA_${flag}.fq"
    rm "ccsAlt_tgt/gapA_${flag}.bed" "ccsAlt_tgt/gapA_${flag}.id"
done < "${REF_NAME}_racon${NPR}_gapA.bed"

# 2.6 Partial Assembly for each Gap
# Software: hifiasm
echo "  - Running partial assemblies..."
mkdir -p partial_assembly
NUM_GAPS=$(wc -l < "${REF_NAME}_racon${NPR}_gapA.bed")
for i in $(seq 1 "$NUM_GAPS"); do
    if [[ -s "longAlt_ref/gapA_${i}.fa" ]] && [[ -s "ccsAlt_tgt/gapA_${i}.fq" ]]; then
        hifiasm --primary -o "partial_assembly/gapA_${i}" --ul "longAlt_ref/gapA_${i}.fa" "ccsAlt_tgt/gapA_${i}.fq" -t "$THREADS"
    else
        echo "    Skipping gap $i (empty inputs)"
    fi
done

# 2.7 Generate Final LongAlt Reads
# Software: awk, seqkit, minimap2, samtools
echo "  - Generating final longAlt reads..."
mkdir -p longAlt_tgt

# Convert GFA to FASTA for each partial assembly
for i in $(seq 1 "$NUM_GAPS"); do
    if [[ -s "longAlt_ref/gapA_${i}.fa" ]] && [[ -s "ccsAlt_tgt/gapA_${i}.fq" ]]; then
        awk '/^S/{print ">"$2;print $3}' "partial_assembly/gapA_${i}.p_ctg.gfa" > "longAlt_tgt/gapA_${i}.p_ctg.fa"
        awk '/^S/{print ">"$2;print $3}' "partial_assembly/gapA_${i}.a_ctg.gfa" > "longAlt_tgt/gapA_${i}.a_ctg.fa"
        cat "longAlt_tgt/gapA_${i}.p_ctg.fa" "longAlt_tgt/gapA_${i}.a_ctg.fa" > "longAlt_tgt/gapA_${i}_pa.fa"
        seqkit seq -n "longAlt_tgt/gapA_${i}_pa.fa" -j "$THREADS" > "longAlt_tgt/gapA_${i}_pa.id"
        rm "longAlt_tgt/gapA_${i}.p_ctg.fa" "longAlt_tgt/gapA_${i}.a_ctg.fa"
    fi
done

# Filter by Read Depth
for i in $(seq 1 "$NUM_GAPS"); do
    if [[ -s "longAlt_ref/gapA_${i}.fa" ]] && [[ -s "ccsAlt_tgt/gapA_${i}.fq" ]]; then
        minimap2 -ax map-hifi "longAlt_tgt/gapA_${i}_pa.fa" "ccsAlt_tgt/gapA_${i}.fq" -t "$THREADS" | samtools view -bS - | samtools sort -o "longAlt_tgt/gapA_${i}_s.bam" -
        samtools depth -a "longAlt_tgt/gapA_${i}_s.bam" > "longAlt_tgt/gapA_${i}_dep.txt"
        awk '$3==0' "longAlt_tgt/gapA_${i}_dep.txt" | awk '{print $1}' | sort | uniq > "longAlt_tgt/gapA_${i}_dep0.id"
        awk 'ARGIND==1{a[$1]} ARGIND==2{if(!($1 in a)){print $1}}' "longAlt_tgt/gapA_${i}_dep0.id" "longAlt_tgt/gapA_${i}_pa.id" > "longAlt_tgt/gapA_${i}_pa_f1.id"
        seqkit grep -f "longAlt_tgt/gapA_${i}_pa_f1.id" "longAlt_tgt/gapA_${i}_pa.fa" >> "longAlt_tgt/gapA_pa_f1_t.fa"
    fi
done

if [ -f "longAlt_tgt/gapA_pa_f1_t.fa" ]; then
    seqkit rename "longAlt_tgt/gapA_pa_f1_t.fa" > "longAlt_tgt/gapA_pa_f1.fa"
    rm "longAlt_tgt/gapA_pa_f1_t.fa"
else
    touch "longAlt_tgt/gapA_pa_f1.fa"
fi

# Filter by Coverage (vs Target Assembly)
minimap2 -x asm5 "../$INITIAL_ASM" "longAlt_tgt/gapA_pa_f1.fa" -t "$THREADS" > "longAlt_tgt/${INITIAL_ASM_NAME}TogapA_pa_f1.paf"
awk '$10/$2>="'$PER_ALT'" && $10/$2<1{print $1}' "longAlt_tgt/${INITIAL_ASM_NAME}TogapA_pa_f1.paf" | sort | uniq > "longAlt_tgt/gapA_pa_f.id"
seqkit grep -f "longAlt_tgt/gapA_pa_f.id" "longAlt_tgt/gapA_pa_f1.fa" > "longAlt_tgt/gapA_pa_f2.fa"

# Filter by Length
seqkit seq -m "$DFL" -j "$THREADS" "longAlt_tgt/gapA_pa_f2.fa" | seqkit fx2tab -l -j "$THREADS" | awk '{print ">longAlt_"NR; print $2}' > "longAlt_tgt/gapA_pa_f3.fa"
ln -sf "longAlt_tgt/gapA_pa_f3.fa" longAlt_tgt.fa

# Visualization stats
pl_lenDis.pl -f longAlt_tgt.fa -split-length "$DFL" -o longAlt_tgt_lenDis.svg || echo "Warning: SVG generation failed"

# Clean up symlinks in current dir
rm "../${REF_NAME}_racon0.fa"
rm "../$INITIAL_ASM_BASE"
rm "../$CCS_BASE"
mv "../${INITIAL_ASM_NAME}.fa.fai" . 2>/dev/null || true

cd ..
LONG_ALT="Alternative_reads/longAlt_tgt.fa"
echo "[Step 2] Done. LongAlt reads generated at $LONG_ALT"
echo "-----------------------------------------------------------------------------"

# =============================================================================
# Step 3: RAGA Optimized Assembly
# =============================================================================
# Purpose: Final assembly using original CCS reads + generated LongAlt reads
# Software: hifiasm, awk

echo "[Step 3] Optimized Assembly"
mkdir -p Optimized_assembly
cd Optimized_assembly

# Link CCS reads
ln -sf "../$CCS" "$CCS_BASE"

# Run hifiasm with ultra-long reads (--ul) integration
hifiasm -o contigs --ul "../$LONG_ALT" "$CCS_BASE" -t "$THREADS"

# Convert GFA to FASTA
awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
ln -sf contigs.bp.p_ctg.fa optimized.fa

# 3.3 Chromosome Scaffolding (Added to ensure chromosome-level output)
# Software: ragtag.py
echo "  - Scaffolding optimized assembly against reference to achieve chromosome-level output..."
ragtag.py scaffold "../$REF" optimized.fa -t "$THREADS" -o .
ln -sf ragtag.scaffold.fasta optimized_chromosome_level.fa

# Clean up
rm "$CCS_BASE"
cd ..

OPTIMIZED_ASM="Optimized_assembly/optimized.fa"
FINAL_ASM="Optimized_assembly/optimized_chromosome_level.fa"
echo "[Step 3] Done. Final optimized assembly (contigs) at $OPTIMIZED_ASM"
echo "        Final chromosome-level assembly at $FINAL_ASM"
echo "============================================================================="
echo "RAGA Workflow Completed Successfully."
