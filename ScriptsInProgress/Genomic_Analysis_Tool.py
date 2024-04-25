#####################################################
#### Title: Genomic Analysis Tool
#### Author: Meghana Sripathi
#### Purpose: This program identifies potential exon edges within a specified genomic region
####          using junction data retrieved from Snaptron, and measures the abundance of 
####          transcription in a specified genomic region using data from recount3 bigWig files.
#### Last Modified: 24/04/2024
####
####################################################

import sys
import requests
import re
import pyBigWig

def parse_genomic_coordinate(genomic_coordinate): 
    match = re.match(r'^(chr\d+):(\d+)-(\d+)_(\S+)$', genomic_coordinate)
    if match:
        chromosome = match.group(1)
        start = match.group(2)
        end = match.group(3)
        strand = match.group(4)
        return chromosome, int(start), int(end), strand
    else:
        print("Invalid genomic coordinate format")
        sys.exit(1)

def retrieve_junction_data(genomic_coordinate, filter_condition=None):
    chromosome, start, end, strand = parse_genomic_coordinate(genomic_coordinate)
    url = f"https://snaptron.cs.jhu.edu/srav2/snaptron?regions={chromosome}:{start}-{end}&rfilter=strand:{strand}"
    if filter_condition:
        url += f"&rfilter={filter_condition}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print("Failed to retrieve junction data from Snaptron")
        sys.exit(1)

def identify_exon_edges(genomic_coordinate, junction_data):
    junctions = [line.split('\t') for line in junction_data.strip().split('\n')[1:]]
    junctions = [(int(j[3]), int(j[4]), float(j[14])) for j in junctions]
    junctions.sort(key=lambda x: x[0])
    exon_edges = []
    current_start = None
    current_end = None
    for start, end, coverage in junctions:
        if coverage > 0.5:
            if current_start is None:
                current_start = start
            current_end = end
        elif current_start is not None:
            exon_edges.append((current_start, current_end))
            current_start = None
            current_end = None
    if current_start is not None:
        exon_edges.append((current_start, current_end))
    return exon_edges

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
            sys.exit(1)
    return abundance

def analyze_genomic_data(genomic_coordinate, analyze_option, filter_condition=None, bigwig_file=None):
    if analyze_option == "Exon Edge Finder":
        try:
            junction_data = retrieve_junction_data(genomic_coordinate, filter_condition)
            exon_edges = identify_exon_edges(genomic_coordinate, junction_data)
            print("Exon Edges:")
            for edge in exon_edges:
                print(f"Start: {edge[0]}, End: {edge[1]}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")
    elif analyze_option == "Genomic Abundance Analyzer":
        try:
            abundance = measure_transcription_abundance(genomic_coordinate, [bigwig_file])
            print(f"Transcription Abundance:\n{bigwig_file}: {abundance[0]}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")

if len(sys.argv) == 2 and sys.argv[1] == "cmd":
    analyze_option = input("Choose Analysis Option (Exon Edge Finder / Genomic Abundance Analyzer): ")
    genomic_coordinate = input("Enter Genomic Coordinate (chr#:start-end_strand): ")
    if analyze_option == "Exon Edge Finder":
        filter_condition = input("Enter Filter Condition (optional): ")
        analyze_genomic_data(genomic_coordinate, analyze_option, filter_condition)
    elif analyze_option == "Genomic Abundance Analyzer":
        bigwig_file = input("Enter BigWig File Path: ")
        analyze_genomic_data(genomic_coordinate, analyze_option, bigwig_file=bigwig_file)
else:
    # Import tkinter for GUI mode
    import tkinter as tk
    from tkinter import messagebox
    from tkinter import filedialog
    from tkinter.scrolledtext import ScrolledText
    def update_ui(selected_option):
        if selected_option == "Exon Edge Finder":
            filter_label.grid(row=2, column=0, padx=5, pady=5)
            filter_entry.grid(row=2, column=1, padx=5, pady=5)
            
            bigwig_button.grid_forget()
            bigwig_entry.grid_forget()
        elif selected_option == "Genomic Abundance Analyzer":
            filter_label.grid_forget()
            filter_entry.grid_forget()
            
            bigwig_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
            bigwig_entry.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

    def analyze_genomic_data():
        if analyze_option.get() == "Exon Edge Finder":
            genomic_coordinate = genomic_entry.get()
            filter_condition = filter_entry.get()
            try:
                junction_data = retrieve_junction_data(genomic_coordinate, filter_condition)
                exon_edges = identify_exon_edges(genomic_coordinate, junction_data)
                
                result_text = "Exon Edges:\n" + "\n".join([f"Start: {edge[0]}, End: {edge[1]}" for edge in exon_edges])
                result_label.config(text=result_text)

                output_window = tk.Toplevel(root)
                output_window.title("Snaptron Output")
                output_text = ScrolledText(output_window, wrap=tk.WORD, width=80, height=20)
                output_text.insert(tk.END, "Snaptron API Output :\n")
                output_text.insert(tk.END, junction_data)
                output_text.pack(fill=tk.BOTH, expand=True)
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred: {str(e)}")
        elif analyze_option.get() == "Genomic Abundance Analyzer":
            genomic_coordinate = genomic_entry.get()
            try:
                file_path = filedialog.askopenfilename(title="Select BigWig file", filetypes=[("BigWig files", "*.bw")])
                abundance = measure_transcription_abundance(genomic_coordinate, [file_path])
                result_text = f"Transcription Abundance:\n{file_path}: {abundance[0]}"
                result_label.config(text=result_text)
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred: {str(e)}")

    root = tk.Tk()
    root.title("Genomic Analysis Tool")

    analyze_option = tk.StringVar()
    analyze_option.set("Exon Edge Finder")
    analyze_option_label = tk.Label(root, text="Choose Analysis Option:")
    analyze_option_label.grid(row=0, column=0, padx=5, pady=5)
    analyze_option_menu = tk.OptionMenu(root, analyze_option, "Exon Edge Finder", "Genomic Abundance Analyzer", command=update_ui)
    analyze_option_menu.grid(row=0, column=1, padx=5, pady=5)

    genomic_label = tk.Label(root, text="Genomic Coordinate (chr#:start-end_strand):")
    genomic_label.grid(row=1, column=0, padx=5, pady=5)
    genomic_entry = tk.Entry(root, width=50)
    genomic_entry.grid(row=1, column=1, padx=5, pady=5)

    filter_label = tk.Label(root, text="Filter Condition (optional):")
    filter_entry = tk.Entry(root, width=50)

    bigwig_button = tk.Button(root, text="Choose BigWig File", command=lambda: bigwig_entry.insert(tk.END, filedialog.askopenfilename(title="Select BigWig file", filetypes=[("BigWig files", "*.bw")])))
    bigwig_entry = tk.Entry(root, width=50)

    analyze_button = tk.Button(root, text="Analyze", command=analyze_genomic_data)
    analyze_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

    result_label = tk.Label(root, text="")
    result_label.grid(row=6, column=0, columnspan=2, padx=5, pady=5)

    update_ui("Exon Edge Finder")  

    root.mainloop()
