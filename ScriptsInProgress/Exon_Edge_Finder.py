#####################################################
#### Title: Exon_Edge_Finder
#### Author: Meghana Sripathi
#### Purpose: This program identifies potential exon edges within a specified genomic region
####          using junction data retrieved from Snaptron.
#### Last Modified: 16/04/2024
####
## This program retrieves junction data from Snaptron for a given genomic region,
## analyzes the data to identify potential exon edges based on coverage or the presence of junctions,
## and outputs the identified exon edges.
####################################################

import requests
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


def retrieve_junction_data(genomic_coordinate, filter_condition=None):
    chromosome, start, end, strand = parse_genomic_coordinate(genomic_coordinate)
    # Query Snaptron with chromosome, int(start), int(end), strand
    url = f"https://snaptron.cs.jhu.edu/srav2/snaptron?regions={chromosome}:{start}-{end}&rfilter=strand:{strand}"
    if filter_condition:
        url += f"&rfilter={filter_condition}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print("Failed to retrieve junction data from Snaptron")
        return None

def identify_exon_edges(genomic_coordinate, junction_data):
    # Parse the junction data
    junctions = [line.split('\t') for line in junction_data.strip().split('\n')[1:]]  # Skip header 
    junctions = [(int(j[3]), int(j[4]), float(j[14])) for j in junctions]  # Extract start, end, and coverage

    # Sort junctions by start position
    junctions.sort(key=lambda x: x[0])

    exon_edges = []
    current_start = None
    current_end = None
    for start, end, coverage in junctions:
        # Check if coverage or presence of junction indicates potential exon boundary
        if coverage > 0.5:  # threshold for coverage
            if current_start is None:
                current_start = start
            current_end = end
        elif current_start is not None:
            exon_edges.append((current_start, current_end))
            current_start = None
            current_end = None
    
    # Add the last potential exon boundary
    if current_start is not None:
        exon_edges.append((current_start, current_end))
    
    return exon_edges

def main():
    genomic_coordinate = input("Enter genomic coordinate (chr#:start-end_strand): ")
    filter_condition = input("Enter filter condition (optional, e.g., samples_count:10): ")
    junction_data = retrieve_junction_data(genomic_coordinate, filter_condition)
    if junction_data:
        exon_edges = identify_exon_edges(genomic_coordinate, junction_data)
        print("Exon edges:", exon_edges)
    else:
        print("Failed to retrieve junction data. Exiting...")

if __name__ == "__main__":
    main()
