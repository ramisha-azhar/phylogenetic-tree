#!/home/ramisha_azhar/anaconda3/bin/python3

import os #joining paths etc 
import glob #finding files with patters ".fastq.gz"
import subprocess #for running external bioinformatic pipelines
import gzip #compressing to .gz
from concurrent.futures import ThreadPoolExecutor #running steps in parallel 
from Bio import SeqIO  #reading fasta


# ============================================
# CONFIGURATION
# ============================================

BASE_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/fasta_files"
READS_DIR = os.path.join(BASE_DIR, "data")
GENOME_DIR = os.path.join(BASE_DIR, "genome")
OUTPUT_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results"
TRIM_DIR = os.path.join(OUTPUT_DIR, "trimmed_reads")
THREADS = 5  # adjust based on CPU cores

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TRIM_DIR, exist_ok=True)

# ============================================
# HELPER FUNCTION
# ============================================

def find_files(directory, patterns):
    """Find all files in 'directory' matching any of the given patterns."""
    files = []
    for pattern in patterns:
        files.extend(glob.glob(os.path.join(directory, pattern)))
    return sorted(files)

def find_reference(genome_dir):
    """Find reference genome file (FASTA/FNA)."""
    for ext in (".fasta", ".fa", ".fna"):
        for file in os.listdir(genome_dir):
            if file.endswith(ext):
                return os.path.join(genome_dir, file)
    raise FileNotFoundError("No reference genome found in GENOME_DIR")

def find_unpaired_reads(trim_dir):
    """Return list of all unpaired FASTQ files."""
    reads = []
    for f in os.listdir(trim_dir):
        if f.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            reads.append(os.path.join(trim_dir, f))
    if not reads:
        raise FileNotFoundError("No unpaired reads found in trimming folder.")
    return sorted(reads)    

# ============================================
# FASTQC + MULTIQC
# ============================================

def run_fastqc(file_list, output_dir, threads):
    """Run FastQC on a list of files."""
    os.makedirs(output_dir, exist_ok=True)

    def worker(fq):
        cmd = ["fastqc", "-o", output_dir, "-t", "1", fq]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"[FASTQC] OK  : {os.path.basename(fq)}")
        except subprocess.CalledProcessError as e:
            print(f"[FASTQC] FAIL: {os.path.basename(fq)} -> {e}")

    with ThreadPoolExecutor(max_workers=threads) as executor:
        executor.map(worker, file_list)

def run_multiqc(input_dir, output_dir):
    """Run MultiQC on a directory of FastQC reports."""
    os.makedirs(output_dir, exist_ok=True)
    cmd = ["multiqc", "-o", output_dir, input_dir]
    try:
        subprocess.run(cmd, check=True)
        print(f"[MultiQC] Report generated at {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"[MultiQC] FAIL: {e}")


# ============================================
# TRIM GALORE (paired-end)
# ============================================

#for removing adapters and low-quality base
def run_trim_galore(input_dir,trim_dir, threads):
    """Run Trim Galore for paired-end reads."""
    os.makedirs(trim_dir, exist_ok=True)
    r1_suffix = "_1.fastq.gz"
    r2_suffix  = "_2.fastq.gz"
    samples = os.listdir(input_dir)
    samples = [ x.split('_')[0] for x in samples] # for each element in the list we split on underscore and take the first element
    for sample in set(samples): 
        r1 = input_dir+'/'+str(sample)+r1_suffix 
        r2 = input_dir+'/'+str(sample)+r2_suffix
        cmd = [
            "trim_galore",
            "--paired",
            "--quality", "20",
            "--length", "30",
            "--trim-n",
            "--cores", str(threads),
            "--output_dir", trim_dir,
            "--fastqc",
            r1, r2]

        print(f"[TrimGalore] Processing sample: {sample}")
        try:
            subprocess.run(cmd, check=True)
            print(f"[TrimGalore] Done: {sample}")
        except subprocess.CalledProcessError as e:
            print(f"[TrimGalore] FAIL: {sample} -> {e}")
        


# ============================================
# MAIN EXECUTION
# ============================================


if __name__ == "__main__":

    # Add these lines here to check what directories your script is using
    print("READS_DIR:", READS_DIR)
    print("GENOME_DIR:", GENOME_DIR)


    # Step 1: Find input files
    read_files = find_files(READS_DIR, ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"])
    genome_files = find_files(GENOME_DIR, ["*.fasta", "*.fa", "*.fna"])
    all_files = read_files + genome_files

    if not all_files:
       print("[ERROR] No FASTA/FASTQ files found. Exiting.")
       exit(1)

    # Step 2: Run QC (FastQC + MultiQC)
    fastqc_dir = os.path.join(OUTPUT_DIR, "fastqc")
    multiqc_dir = os.path.join(OUTPUT_DIR, "multiqc")
    run_fastqc(read_files, fastqc_dir, THREADS)
    run_multiqc(fastqc_dir, multiqc_dir)

    # Step 3: Run Trim Galore (for reads only)
    trim_dir = os.path.join(OUTPUT_DIR, "trimmed_reads")
    run_trim_galore(READS_DIR,trim_dir, THREADS)
    print(f"[INFO] QC + trimming complete! Results saved in {OUTPUT_DIR}")

    # Step 4: Move fastqc files in folder 
    fastqc_trimmed = os.listdir(trim_dir)
    fastqc_trimmed = [ trim_dir+'/'+x for x in fastqc_trimmed if 'fastqc' in x]
    fastq_trimmed_folder = trim_dir+'/fastqc'
    os.mkdir(fastq_trimmed_folder)
    for x in fastqc_trimmed:
        os.system('mv '+x+' '+fastq_trimmed_folder)

# ============================================
# HOMOLOGY CHECK (MASH)
# ============================================

# Folders
TRIMMED_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results/trimmed_reads"
GENOME_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/fasta_files/genome"
OUTPUT_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results"
MASH_DIR = os.path.join(OUTPUT_DIR, "mash_results")

# Create output folder
os.makedirs(MASH_DIR, exist_ok=True)

# ============================================
# INPUT FILES
# ============================================

# Trimmed reads (FASTQ)
trimmed_files = sorted(glob.glob(os.path.join(TRIMMED_DIR, "*_val*.fq.gz")))

if len(trimmed_files) < 1:
    print("[ERROR] No trimmed reads found!")
    exit(1)

# Reference genome (FASTA)
reference_genomes = sorted(glob.glob(os.path.join(GENOME_DIR, "*.fasta")) +
                           glob.glob(os.path.join(GENOME_DIR, "*.fna")))

if len(reference_genomes) == 0:
    print("[ERROR] No reference genome found!")
    exit(1)

# Use first reference genome
reference = reference_genomes[0]
print(f"[INFO] Using {os.path.basename(reference)} as reference genome for Mash.")

# ============================================
# CREATE MASH SKETCH FOR REFERENCE GENOME
# ============================================

ref_sketch = os.path.join(MASH_DIR, os.path.basename(reference) + ".msh")
print("[INFO] Creating Mash sketch for reference genome...")
subprocess.run(["mash", "sketch", "-o", ref_sketch, reference], check=True)

# ============================================
# CREATE MASH SKETCHES FOR TRIMMED READS AND COMPARE
# ============================================

for query in trimmed_files:
    query_sketch = os.path.join(MASH_DIR, os.path.basename(query) + ".msh")
    # Sketch the read
    subprocess.run(["mash", "sketch", "-o", query_sketch, query], check=True)

    # Output file for Mash distance
    out_file = os.path.join(MASH_DIR, f"{os.path.basename(query)}_vs_ref.txt")
    print(f"[INFO] Comparing {os.path.basename(query)} to reference genome using Mash...")
    with open(out_file, "w") as f_out:
        subprocess.run(["mash", "dist", ref_sketch, query_sketch],
                       stdout=f_out, check=True)

print(f"\nMash homology check complete! Results saved in {MASH_DIR}")


##### PART 2 #####
# ================================
# CONFIGURATION
# ================================

BASE_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/fasta_files"
READS_DIR = os.path.join(BASE_DIR, "data")
GENOME_DIR = os.path.join(BASE_DIR, "genome")
OUTPUT_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results"
TRIM_DIR = os.path.join(OUTPUT_DIR, "trimmed_reads")
THREADS = 5  # adjust based on CPU cores
ALIGN_DIR = os.path.join(OUTPUT_DIR, "alignment")
STATS_DIR = os.path.join(ALIGN_DIR, "alignment_stats")

PICARD = "/home/razhar/tools/picard.jar"   # ← UPDATE THIS PATH

# Create folders
os.makedirs(TRIM_DIR, exist_ok=True)
os.makedirs(ALIGN_DIR, exist_ok=True)
os.makedirs(STATS_DIR, exist_ok=True)


# ================================
# FIND SAMPLES + REFERENCE
# ================================
def find_reference(genome_dir):
    for ext in (".fasta", ".fa", ".fna"):
        for f in os.listdir(genome_dir):
            if f.endswith(ext):
                return os.path.join(genome_dir, f)
    raise FileNotFoundError("Reference genome not found.")


def find_trimmed_pairs(trim_dir):
    r1_files = sorted(glob.glob(os.path.join(trim_dir, "*_1_val_1.fq")))
    pairs = []

    for r1 in r1_files:
        r2 = r1.replace("_1_val_1.fq", "_2_val_2.fq")
        if os.path.exists(r2):
            pairs.append((r1, r2))
        else:
            print(f"[WARNING] Missing R2 for {r1}")

    if not pairs:
        raise FileNotFoundError("No trimmed paired FASTQ files found.")

    return pairs


# ================================
# ALIGNMENT (BWA + SAMTOOLS)
# ================================
def align_trimmed_pairs(r1, r2, reference, output_dir, threads):

    sample = os.path.basename(r1).replace("_1_val_1.fq", "")
    prefix = os.path.join(output_dir, sample)

    sam_file = prefix + ".sam"
    bam_file = prefix + ".bam"
    sorted_bam = prefix + "_sorted.bam"

    ref_copy = os.path.join(output_dir, os.path.basename(reference))

    try:
        # Copy reference FASTA if needed
        if not os.path.exists(ref_copy):
            print(f"[COPY] Copying reference to {output_dir}")
            subprocess.run(["cp", reference, ref_copy], check=True)

        # Build BWA index if missing
        if not os.path.exists(ref_copy + ".bwt"):
            print("[INDEX] Indexing reference with BWA...")
            subprocess.run(["bwa", "index", ref_copy], check=True)

        # Read group tag
        rg_tag = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"

        # Alignment
        print(f"[ALIGN] Aligning sample: {sample}")
        with open(sam_file, "w") as sam_out:
            subprocess.run([
                "bwa", "mem", "-t", str(threads),
                "-R", rg_tag,
                ref_copy, r1, r2
            ], stdout=sam_out, check=True)

        # Convert SAM → BAM → Sorted BAM
        print(f"[SORT] Sorting BAM for {sample}")
        subprocess.run(["samtools", "view", "-bS", sam_file, "-o", bam_file], check=True)
        subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)

        # Index sorted BAM
        subprocess.run(["samtools", "index", sorted_bam], check=True)

        # Cleanup
        os.remove(sam_file)
        os.remove(bam_file)

        print(f"[DONE] Alignment finished → {sorted_bam}")

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Alignment failed for {sample}: {e}")


# ================================
# POST-PROCESSING (FLAGSTAT + PICARD)
# ================================
def postprocess(bam_file):

    sample = os.path.basename(bam_file).replace("_sorted.bam", "")
    stats_file = os.path.join(STATS_DIR, f"{sample}_flagstat.txt")
    dedup_bam = os.path.join(ALIGN_DIR, f"{sample}_dedup.bam")
    metrics_file = os.path.join(STATS_DIR, f"{sample}_dedup_metrics.txt")

    try:
        # Flagstat
        print(f"[FLAGSTAT] {sample}")
        with open(stats_file, "w") as out:
            subprocess.run(["samtools", "flagstat", bam_file], stdout=out, check=True)

        # Picard Duplicate Removal
        print(f"[PICARD] Removing duplicates for {sample}")
        subprocess.run([
            "java", "-jar", PICARD, "MarkDuplicates",
            f"I={bam_file}",
            f"O={dedup_bam}",
            f"M={metrics_file}",
            "REMOVE_DUPLICATES=true",
            "ASSUME_SORTED=true"
        ], check=True)

        # Index deduplicated BAM
        subprocess.run(["samtools", "index", dedup_bam], check=True)

        print(f"[DEDUP] Done → {dedup_bam}")

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Post-processing failed for {sample}: {e}")


# ================================
# MAIN PIPELINE
# ================================
if __name__  == "__main__":

    reference = find_reference(GENOME_DIR)
    pairs = find_trimmed_pairs(TRIM_DIR)

    print(f"Using reference: {reference}")
    print(f"Found {len(pairs)} paired-samples\n")

    # Run alignments
    with ThreadPoolExecutor(max_workers=THREADS) as exe:
        for r1, r2 in pairs:
            exe.submit(align_trimmed_pairs, r1, r2, reference, ALIGN_DIR, THREADS)

    print("\n[INFO] Alignment Finished\n")

    # Post-processing
    bam_files = sorted(glob.glob(os.path.join(ALIGN_DIR, "*_sorted.bam")))
    for bam in bam_files:
        postprocess(bam)

    print("\n[INFO] PIPELINE COMPLETED SUCCESSFULLY!")



# ============================================
# CLEAN VARIANT CALLING PIPELINE (bcftools)
# Using dedup BAM files + generating VCF + gVCF
# ============================================

import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor
import gzip

# ---------------------------
# DIRECTORIES
# ---------------------------
BASE_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/fasta_files"
OUTPUT_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results"
READS_DIR = os.path.join(BASE_DIR, "data")
GENOME_DIR = os.path.join(BASE_DIR, "genome")
ALIGN_DIR = os.path.join(OUTPUT_DIR, "alignment")    # contains *_dedup.bam)
VARIANT_DIR = os.path.join(OUTPUT_DIR, "variant")
THREADS = 4

os.makedirs(VARIANT_DIR, exist_ok=True)


# -----------------------------------------------------
# 1. Ensure reference FASTA has .fai index (required)
# -----------------------------------------------------
def ensure_reference_index(reference):
    fai = reference + ".fai"
    if not os.path.exists(fai):
        print(f"[INDEX] Creating index for reference: {reference}")
        subprocess.run(["samtools", "faidx", reference], check=True)
    else:
        print("[INDEX] Reference index already exists.")


# -----------------------------------------------------
# 2. Call gVCF AND standard VCF for each sample
# -----------------------------------------------------
def call_variants(bam, reference, outdir, block_size=10):
    sample = os.path.basename(bam).replace("_dedup.bam", "")
    
    gvcf = os.path.join(outdir, f"{sample}.g.vcf.gz")
    vcf  = os.path.join(outdir, f"{sample}.vcf.gz")

    print(f"\n[CALLING] Processing sample: {sample}")

    # -----------------------------
    # gVCF CALLING
    # -----------------------------
    print(f"[gVCF] Calling {gvcf}")

    mpileup_gvcf = [
        "bcftools", "mpileup",
        "-Ou",
        "-f", reference,
        "--gvcf", str(block_size),
        "-A",
        bam
    ]

    call_gvcf_cmd = [
        "bcftools", "call",
        "-m",       # retains reference blocks
        "-Oz",
        "-o", gvcf
    ]

    p1 = subprocess.Popen(mpileup_gvcf, stdout=subprocess.PIPE)
    subprocess.run(call_gvcf_cmd, stdin=p1.stdout, check=True)
    p1.stdout.close()
    p1.wait()
    subprocess.run(["bcftools", "index", gvcf], check=True)

    # -----------------------------
    # STANDARD VCF CALLING
    # -----------------------------
    print(f"[VCF] Calling {vcf}")

    mpileup_vcf = [
        "bcftools", "mpileup",
        "-Ou",
        "-f", reference,
        "-A",
        bam
    ]

    call_vcf_cmd = [
        "bcftools", "call",
        "-mv",      # multiallelic, variants only (no ref blocks)
        "-Oz",
        "-o", vcf
    ]

    p2 = subprocess.Popen(mpileup_vcf, stdout=subprocess.PIPE)
    subprocess.run(call_vcf_cmd, stdin=p2.stdout, check=True)
    p2.stdout.close()
    p2.wait()
    subprocess.run(["bcftools", "index", vcf], check=True)

    return gvcf, vcf

# ============================================
#  SNP MATRIX FROM VCF.GZ FILES
# ============================================
import os
import glob
import subprocess

OUTPUT_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results"
VARIANT_DIR = os.path.join(OUTPUT_DIR, "variant")

def build_snp_matrix_from_gz():
    """Merge VCFs and convert to FASTA using vcf2phylip."""

    # Find all .vcf.gz files EXCEPT old merged files
    vcf_gz_files = sorted(
        f for f in glob.glob(os.path.join(VARIANT_DIR, "*.vcf.gz"))
        if "merged" not in os.path.basename(f)
    )

    if len(vcf_gz_files) < 2:
        print("[ERROR] Need at least two .vcf.gz files!")
        return

    print(f"[INFO] Found {len(vcf_gz_files)} VCFs to merge.")

    merged_vcf = os.path.join(VARIANT_DIR, "merged.vcf.gz")

    # Merge VCFs
    print("[INFO] Merging VCFs...")
    subprocess.run(
        ["bcftools", "merge", "-Oz", "-o", merged_vcf] + vcf_gz_files,
        check=True
    )

    # Index merged VCF
    print("[INFO] Indexing merged VCF...")
    subprocess.run(["bcftools", "index", merged_vcf], check=True)

    # Ensure vcf2phylip writes output into VARIANT_DIR
    os.chdir(VARIANT_DIR)

    # Force FASTA output
    print("[INFO] Converting merged VCF to FASTA...")
    subprocess.run(
        ["vcf2phylip.py", "-i", merged_vcf, "--output-fasta"],
        check=True
    )

    # Check if FASTA was created
    fasta_path = os.path.join(VARIANT_DIR, "merged.min4.fasta")
    if os.path.exists(fasta_path):
        print("[SUCCESS] FASTA created:", fasta_path)
    else:
        print("[WARNING] vcf2phylip did NOT create merged.min4.fasta")
        print("          But merged.min2.phy is available.")

if __name__ == "__main__":
    build_snp_matrix_from_gz()

'''
# ===============
#  Build IQ‑TREE
# ===============
import os
import subprocess

OUTPUT_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results"
VARIANT_DIR = os.path.join(OUTPUT_DIR, "variant")
PHYLO_DIR = os.path.join(OUTPUT_DIR, "phylo")
OUTPUT_DIR = os.path.join(VARIANT_DIR, "merged.min2.fasta")

def run_iqtree():
    # Make sure phylo directory exists
    os.makedirs(PHYLO_DIR, exist_ok=True)

    # Change into phylo so all IQ-TREE outputs land there
    os.chdir(PHYLO_DIR)

    cmd = [
        "iqtree2",
        "-s", INPUT_FASTA,  # alignment
        "-m", "MFP",        # ModelFinder Plus
        "-bb", "1000",      # ultrafast bootstrap
        "-nt", "AUTO"       # auto threads
    ]

    print("[INFO] Running IQ-TREE...")
    subprocess.run(cmd, check=True)
    print("[DONE] IQ-TREE finished. Tree file is in:", PHYLO_DIR)

if __name__ == "__main__":
    run_iqtree()
'''

# ====================================================================
#  compute_snp_distance,build_distance_matrix and build_two_taxon_tree
# ====================================================================   

import os
import subprocess
import math

VCF = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results/variant/merged.vcf.gz"
SAMPLE1 = "ERR4079133"
SAMPLE2 = "ERR4079285"
OUTPUT_DIR = "/mnt/c/Users/ramis/OneDrive/Desktop/BPP/phylogenetic tree/qc_results"
VARIANT_DIR = os.path.join(OUTPUT_DIR, "variant")
PHYLO_DIR = os.path.join(OUTPUT_DIR, "phylo")

def compute_snp_distance():
    os.makedirs(PHYLO_DIR, exist_ok=True)
    os.chdir(PHYLO_DIR)

    print("[INFO] Reading VCF and comparing genotypes...")

    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%GT\t]\n", VCF]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    total = 0 #number of SNPs where both samples have a valid genotyp
    diff = 0 #number of SNPs where the genotypes differ

    for line in result.stdout.splitlines():
        fields = line.strip().split()
        gt1, gt2 = fields[-2], fields[-1]

        if "." in gt1 or "." in gt2:
            continue

        total += 1
        if gt1 != gt2:
            diff += 1

    identity = (1 - diff / total) * 100

    with open("pairwise_stats.txt", "w") as f:
        f.write(f"Total comparable SNPs: {total}\n")
        f.write(f"Different SNPs: {diff}\n")
        f.write(f"Percent identity: {identity:.4f}%\n")

    print("\n=== SNP DISTANCE RESULTS ===")
    print("Total comparable SNPs:", total)
    print("Different SNPs:", diff)
    print("Percent identity:", round(identity, 4), "%")

    return total, diff, identity


def build_distance_matrix(diff, total):
    print("\n=== DISTANCE MATRIX ===")

    p = diff / total

    matrix = (
        f"                {SAMPLE1}     {SAMPLE2}\n"
        f"{SAMPLE1}        0.0000       {p:.6f}\n"
        f"{SAMPLE2}        {p:.6f}       0.0000\n"
    )

    print("\n" + matrix)

    with open("distance_matrix.txt", "w") as f:
        f.write(matrix)

    return p


def build_two_taxon_tree(p):
    print("\n=== 2‑TAXON NEWICK TREE ===")

    branch = p / 2
    newick = f"({SAMPLE1}:{branch:.6f},{SAMPLE2}:{branch:.6f});"

    print("\nNewick tree:")
    print(newick)

    with open("two_taxon_tree.nwk", "w") as f:
        f.write(newick)

    print("\nTree saved to:", os.path.join(PHYLO_DIR, "two_taxon_tree.nwk"))


if __name__ == "__main__":
    total, diff, identity = compute_snp_distance()
    p = build_distance_matrix(diff, total)
    build_two_taxon_tree(p)
