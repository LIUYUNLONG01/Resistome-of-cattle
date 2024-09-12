#taxonomic annotation of ruminal microbiota
#!/bin/bash
DB_DIR="/kraken2/db"
FASTQ_DIR="/rawdata"
OUTPUT_DIR="/taxonomy"

mkdir -p $OUTPUT_DIR

for file1 in $FASTQ_DIR/*_remove_host_1.fastq
do
    file2="${file1/_1.fastq/_2.fastq}"
    sampleID=$(basename $file1)
    sampleID="${sampleID%%_*}"
    output="${OUTPUT_DIR}/${sampleID}_kraken_output.txt"
    report="${OUTPUT_DIR}/${sampleID}_kraken_report.txt"

    kraken2 --db $DB_DIR --threads 20 --confidence 0.2 --output $output --report $report --paired $file1 $file2

    echo "Processed $sampleID"
done


#!/bin/bash

BRACKEN_DB="/kraken2/db"
KRAKEN_OUTPUT_DIR="/taxonomy"
BRACKEN_OUTPUT_DIR="/taxonomy/bracken"

mkdir -p $BRACKEN_OUTPUT_DIR

for report in $KRAKEN_OUTPUT_DIR/*_kraken_report.txt
do
    sampleID=$(basename $report)
    sampleID="${sampleID%_kraken_report.txt}"
    bracken_output="${BRACKEN_OUTPUT_DIR}/${sampleID}_bracken_output.txt"

    bracken -d $BRACKEN_DB -i $report -o $bracken_output -r 150 -l P -t 20

    echo "Bracken processed for $sampleID"
done

#The microbial hosts of ARG
## Read DIAMOND alignment results
def parse_diamond_output(file_path, output_file_path):
    with open(file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        output_file.write("Query ID,Subject ID\n")
        for line in input_file:
            parts = line.strip().split('\t')
            query_id = parts[0]  
            subject_id = parts[1] 
            output_file.write(f"{query_id},{subject_id}\n")

diamond_results_path = '/ARG/All_gene_ARG.besthit.dmnd.txt'
output_csv_path = '/ARG/parsed_diamond_results.csv'

parse_diamond_output(diamond_results_path, output_csv_pa)

## Analyzing
def read_diamond_results(diamond_results_path):
    matched_orf_ids = set()
    with open(diamond_results_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            orf_id = parts[0]
            matched_orf_ids.add(orf_id)
    return matched_orf_ids

def write_matched_orf_ids_to_file(matched_orf_ids, matched_orf_ids_path):
    with open(matched_orf_ids_path, 'w') as file:
        for orf_id in matched_orf_ids:
            file.write(f"{orf_id}\n")

diamond_results_path = '/ARG/parsed_diamond_results.csv'
matched_orf_ids_path = '/ARG/matched_orf_ids.txt'

matched_orf_ids = read_diamond_results(diamond_results_path)
write_matched_orf_ids_to_file(matched_orf_ids, matched_orf_ids_path)

##Extract ORF location information from GFF files
import os

def parse_gff(gff_file):
    orf_info = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if parts[2] == "CDS":
                contig = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[-1]
                gene_id = [attr for attr in attributes.split(';') if "ID=" in attr][0].split('=')[1]
                orf_info[gene_id] = {'contig': contig, 'start': start, 'end': end, 'strand': strand}
    return orf_info

def parse_all_gff(folder_path):
    all_gff_info = {}
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".gff"):
            gff_file_path = os.path.join(folder_path, file_name)
            orf_info = parse_gff(gff_file_path)
            all_gff_info.update(orf_info)
    return all_gff_info

folder_path = "/ORF" 
output_csv_path = os.path.join(folder_path, "orf_positions.csv") 

all_orf_info = parse_all_gff(folder_path)

with open(output_csv_path, 'w') as f:
    f.write("ORF ID,Contig,Start,End,Strand\n")
    for orf_id, info in all_orf_info.items():
        f.write(f"{orf_id},{info['contig']},{info['start']},{info['end']},{info['strand']}\n")

##selected
 import csv

def filter_orf_positions_based_on_contig(orf_positions_csv_path, matched_orf_ids_path, filtered_orf_positions_path):
    with open(matched_orf_ids_path, 'r') as f:
        matched_contig_ids = {line.strip().rsplit('_', 1)[0] for line in f} 
    with open(orf_positions_csv_path, 'r') as input_csv, open(filtered_orf_positions_path, 'w', newline='') as output_csv:
        reader = csv.DictReader(input_csv)
        writer = csv.DictWriter(output_csv, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row['Contig'] in matched_contig_ids:
                writer.writerow(row)

orf_positions_csv_path = '/ORF/orf_positions.csv'
matched_orf_ids_path = '/ARG/matched_orf_ids.txt'
filtered_orf_positions_path = '/ARG/filtered_orf_positions.csv'

filter_orf_positions_based_on_contig(orf_positions_csv_path, matched_orf_ids_path, filtered_orf_positions_path)

##extracted sequences
from Bio import SeqIO
import csv
import os

def extract_orf_sequences_batch(contigs_batch_path, filtered_orf_positions_path, output_fasta_batch_path):
    contigs = {}
    for filename in os.listdir(contigs_batch_path):
        if filename.endswith(".fa"): 
            filepath = os.path.join(contigs_batch_path, filename)
            for record in SeqIO.parse(filepath, "fasta"):
                contigs[record.id] = record.seq

    with open(filtered_orf_positions_path, newline='') as positions_file, open(output_fasta_batch_path, 'w') as output_fasta:
        reader = csv.DictReader(positions_file)
        for row in reader:
            contig_id = row['Contig']
            if contig_id in contigs:
                start, end = int(row['Start']), int(row['End'])
                strand, orf_id = row['Strand'], row['ORF ID']
                orf_seq = contigs[contig_id][start-1:end]
                if strand == '-':
                    orf_seq = orf_seq.reverse_complement()
                SeqIO.write(SeqIO.SeqRecord(orf_seq, id=orf_id, description=""), output_fasta, "fasta")

batch_folder = '/contig'
filtered_orf_positions_path = '/ARG/filtered_orf_positions.csv'
output_fasta_batch_path = '/ARG/extracted_orfs_batch1.fasta'

extract_orf_sequences_batch(batch_folder, filtered_orf_positions_path, output_fasta_batch_path)

##annotation
time diamond blastx -d /nr.dmnd -q /ARG/extracted_orfs.fasta -o /ARG/blast/matches.m8 --more-sensitive

##selected
from Bio import SeqIO
import csv
import os

def filter_diamond_results(input_file, output_file, evalue_threshold=1e-5):
    best_hits = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            columns = line.strip().split('\t')
            query_id = columns[0]
            evalue = float(columns[10])
            if query_id not in best_hits or evalue < float(best_hits[query_id][10]):
                best_hits[query_id] = columns
   
    with open(output_file, 'w') as outfile:
        for hit in best_hits.values():
            outfile.write('\t'.join(map(str, hit)) + '\n')

input_file_path = '/ARG/blast/matches.m8' 
output_file_path = '/ARG/blast/filtered_matches.m8'

filter_diamond_results(input_file_path, output_file_path)

##extracted taxonomic annotation 
import time
from Bio import Entrez
from tqdm import tqdm

def fetch_taxonomy(sequence_ids, batch_size=20):

    taxonomy_info = {}
    for i in tqdm(range(0, len(sequence_ids), batch_size), desc="Processing batches"):
        batch = sequence_ids[i:i+batch_size]
        try:
            handle = Entrez.efetch(db="protein", id=",".join(batch), rettype="gb", retmode="text")
            records = handle.read()
            for seq_id in batch:
                record_start = records.find(seq_id)
                record_end = records.find("//", record_start) + 2
                record_text = records[record_start:record_end]
                taxonomy = record_text.split("ORGANISM")[1].split("\n", 1)[1].split(";")
                taxonomy = [tax.strip() for tax in taxonomy]
               
                taxonomy_dict = {
                    "Phylum": taxonomy[0] if len(taxonomy) > 1 else "",
                    "Class": taxonomy[1] if len(taxonomy) > 2 else "",
                    "Order": taxonomy[2] if len(taxonomy) > 3 else "",
                    "Family": taxonomy[3] if len(taxonomy) > 4 else "",
                    "Genus": taxonomy[4] if len(taxonomy) > 5 else "",
                    "Species": taxonomy[-1] if len(taxonomy) > 0 else "" 
                }
                taxonomy_info[seq_id] = taxonomy_dict
            handle.close()
            time.sleep(0.3) 
        except Exception as e:
            print(f"Error processing batch {i//batch_size+1}: {e}")
            time.sleep(5)
    return taxonomy_info

with open("extracted_sequence_ids.txt", "r") as file:
    sequence_ids = [line.strip() for line in file]


taxonomy_info = fetch_taxonomy(sequence_ids)


with open("taxonomy_annotation.txt", "w") as output_file:
    for seq_id, taxonomy in taxonomy_info.items():
        output_file.write(f"{seq_id}\t{taxonomy['Phylum']}\t{taxonomy['Class']}\t{taxonomy['Order']}\t{taxonomy['Family']}\t{taxonomy['Genus']}\t{taxonomy['Species']}\n")



