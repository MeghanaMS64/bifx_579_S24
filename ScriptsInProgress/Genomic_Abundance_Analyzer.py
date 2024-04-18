#####################################################
#### Title: Genomic_Abundance_Analyzer
#### Author: Meghana Sripathi
#### Purpose: This script measures the abundance of transcription in a specified genomic region
####          using data from recount3 bigWig files.
#### Last Modified: 16/04/2024
####
## This script takes a genomic coordinate and paths to recount3 bigWig files as input,
## measures the transcription abundance within the specified genomic region using the provided bigWig files,
## and prints the calculated transcription abundances.
####################################################

import pyBigWig
import re
import sys

# Old function to parse the input, replaced by regex 
# def parse_genomic_coordinate(genomic_coordinate):
#     chromosome, position, strand = genomic_coordinate.split(':')
#     start, end = position.split('-')
#     return chromosome, int(start), int(end), strand

# new parsing function
def parse_genomic_coordinate(genomic_coordinate): 
    match = re.match(r'^(chr\d+):(\d+)-(\d+)_(\S+)$', genomic_coordinate)
    if match:
        chromosome = match.group(1)
        start = match.group(2)
        end = match.group(3)
        strand = match.group(4)
        return chromosome, int(start), int(end), strand
    else:
        print("Something went wrong in the Regex, please check the format again")
        sys.exit(0) # To avoid processing bad formatted inputs we just exit 
        return None

def measure_transcription_abundance(genomic_coordinate, bigwig_files):
    chromosome, start, end, _ = parse_genomic_coordinate(genomic_coordinate)
    abundance = []
    for bw_file in bigwig_files:
        bw = pyBigWig.open(bw_file)
        if bw is not None:
            values = bw.values(chromosome, start, end)
            if values:
                abundance.append(sum(filter(lambda x: x is not None, values)))
            else:
                abundance.append(0)
        else:
            print(f"Failed to open bigWig file: {bw_file}")
            abundance.append(None)
    return abundance

def main():
    genomic_coordinate = input("Enter genomic coordinate (chr#:start-end_strand): ")
    bigwig_files = input("Enter path to recount3 bigWig files (separated by comma): ").split(',')
    abundance = measure_transcription_abundance(genomic_coordinate, bigwig_files)
    print("Transcription abundance:", abundance)

if __name__ == "__main__":
    main()
