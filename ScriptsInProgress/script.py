##############################################################

#### Title: First Script from pseudo-code

#### Author: Meghana Sripathi

#### Purpose: just a rough code from the pseudo code i made so that i can start refining it and make it functionnal 

#### Last modified: 27/03/2024

##############################################################


# main():
#    - Input: RNAseq data file path, genome annotation file path
#    - Load RNAseq data and genome annotation
#    - Identify expressed regions
#    - Output: List of expressed genomic regions

# load_rnaseq_data(rnaseq_file):
#    - Input: RNAseq data file path
#    - Load aligned reads from the RNAseq data file
#    - Return aligned reads

# load_genome_annotation(annotation_file):
#    - Input: Genome annotation file path
#    - Load gene and intron annotations from the annotation file
#    - Return gene and intron annotations

# identify_expressed_regions(rnaseq_reads, gene_annotations, intron_annotations):
#    - Input: Aligned reads from RNAseq data, gene annotations, intron annotations
#    - Define parameters: coverage_threshold, ratio_threshold, min_region_length
#    - Iterate through each genomic region:
#      - Check if the region is within an intron or between genes
#      - If it is, extract reads mapping to the region
#      - Calculate read coverage for the region
#      - Check if read coverage exceeds coverage_threshold and ratio_threshold
#      - If yes, classify the region as expressed
#    - Return list of expressed genomic regions

# extract_reads_mapping_to_region(region, aligned_reads):
#    - Input: Genomic region, aligned reads
#    - Extract reads that map to the specified genomic region
#    - Return extracted reads

# calculate_read_coverage(region, mapped_reads):
#    - Input: Genomic region, mapped reads
#    - Calculate read coverage (reads per base pair) for the specified region
#    - Return read coverage

# output_expressed_regions(expressed_regions):
#    - Input: List of expressed genomic regions
#    - Output the list of expressed genomic regions to a file or print to console

# Example:
#    - rnaseq_file = "rnaseq_data"
#    - annotation_file = "genome_annotation"
#    - rnaseq_reads = load_rnaseq_data(rnaseq_file)
#    - gene_annotations, intron_annotations = load_genome_annotation(annotation_file)
#    - expressed_regions = identify_expressed_regions(rnaseq_reads, gene_annotations, intron_annotations)
#    - output_expressed_regions(expressed_regions)



import numpy as np
import pysam
import pyBigWig
import pandas as pd
import requests

def main(rnaseq_file, annotation_file):
    rnaseq_reads = load_rnaseq_data(rnaseq_file)
    gene_annotations, intron_annotations = load_genome_annotation(annotation_file)
    expressed_regions = identify_expressed_regions(rnaseq_reads, gene_annotations, intron_annotations)
    for region in expressed_regions:
        splice_junctions = measure_expression_and_find_splice_junctions(region)
        print("Splice junctions for genomic region:", splice_junctions)

def load_rnaseq_data(rnaseq_file):
    # Load aligned reads from the RNAseq data file
    aligned_reads = []
    with pysam.AlignmentFile(rnaseq_file, "rb") as samfile:
        for read in samfile.fetch():
            aligned_reads.append(read)
    return aligned_reads

def load_genome_annotation(annotation_file):
    # Load gene and intron annotations from the annotation file
    gene_annotations = {}
    intron_annotations = []
    with open(annotation_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if parts[2] == 'gene':
                gene_name = parts[8].split(';')[0].split('=')[1]
                gene_start = int(parts[3])
                gene_end = int(parts[4])
                gene_annotations[gene_name] = (gene_start, gene_end)
            elif parts[2] == 'intron':
                intron_start = int(parts[3])
                intron_end = int(parts[4])
                intron_annotations.append((intron_start, intron_end))
    return gene_annotations, intron_annotations

def identify_expressed_regions(rnaseq_reads, gene_annotations, intron_annotations):
    coverage_threshold = 10
    ratio_threshold = 0.5
    min_region_length = 50

    expressed_regions = []

    for region in gene_annotations.values():
        if is_expressed_region(region, intron_annotations):
            mapped_reads = extract_reads_mapping_to_region(region, rnaseq_reads)
            region_coverage = calculate_read_coverage(region, mapped_reads)
            if region_coverage >= coverage_threshold:
                expressed_regions.append(region)

    return expressed_regions

def is_expressed_region(region, intron_annotations):
    for intron in intron_annotations:
        if region[0] >= intron[0] and region[1] <= intron[1]:
            return False
    return True

def extract_reads_mapping_to_region(region, aligned_reads):
    mapped_reads = []
    for read in aligned_reads:
        if region[0] <= read.reference_start and region[1] >= read.reference_end:
            mapped_reads.append(read)
    return mapped_reads

def calculate_read_coverage(region, mapped_reads):
    region_length = region[1] - region[0]
    coverage = len(mapped_reads) / region_length
    return coverage

def measure_expression_and_find_splice_junctions(genomic_region):
    bigWig_file = "coverage.bw"
    monorail_file = "monorail_summaries.tsv"
    megadepth_file = "megadepth_summaries.tsv"
    
    # Query Snaptron to find splice junctions for the genomic region of interest
    splice_junctions = snaptron_query(genomic_region)
    return splice_junctions

def snaptron_query(genomic_region):
    # Query Snaptron API
    snaptron_api_url = "https://snaptron.cs.jhu.edu/srav2"
    params = {"region": genomic_region}
    try:
        response = requests.get(snaptron_api_url, params=params)
        if response.status_code == 200:
            data = response.json()
            splice_junctions = data.get("splice_junctions", [])
            return splice_junctions
        else:
            print("Snaptron query failed:", response.text)
    except requests.RequestException as e:
        print("Error during Snaptron query:", e)
    return []

# usage
rnaseq_file = "rnaseq_data.bam"
annotation_file = "genome_annotation.gff"
main(rnaseq_file, annotation_file)
