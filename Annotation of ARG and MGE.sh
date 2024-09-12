#Annotation of ARG based on reads
time args_oap stage_one -i /contigs/ -o /argoap/output -f fa -t 20

time args_oap stage_two -i /argoap/output -t 20

#Annotation of MGE based on reads
args_oap stage_one -i /contigs/ -o /argoap/output/mge -f fa -t 20 --database /argoap/MGEs_FINAL_99perc_trim.fasta

args_oap stage_two -i /argoap/output/mge -t 20 --database /argoap/MGEs_FINAL_99perc_trim.fasta --structure1 /argoap/structure.txt

#Predicted ORF
time prodigal -a /ORF/${names[${SLURM_ARRAY_TASK_ID}]}.faa -d ${names[${SLURM_ARRAY_TASK_ID}]}_nucleotide_seq.fna -o ${names[${SLURM_ARRAY_TASK_ID}]}_genes_gff -i /contig/${names[${SLURM_ARRAY_TASK_ID}]}.final.contigs.fa -f gff -p meta

#remove redundant sequences
cat /ORF/*.faa > total_protein.faa
time cd-hit -i total_protein.faa -o total_protein_90.faa -c 0.90 -s 0.8 -n 5 -g 1 -d 0 -T 32


#mapping-ARG
time diamond blastp -q /ORF/total_protein_90.faa -d /ARG/db/protein_homolog.dmnd -p 32 -e 1e-5 -k 1 --id 80 --query-cover 80 --sensitive -o All_gene_ARG.besthit.dmnd.txt

#convert nucleic acid sequence to amino acid sequence
transeq MGEs_FINAL_99perc_trim.fasta MGEs_FINAL_99perc_trim_protein.fasta -frame=6

#build index
time diamond makedb --in MGEs_FINAL_99perc_trim_protein.fasta -d MGE -p 8

#mapping
time diamond blastp -q /ORF/total_protein_90.faa -d /MGE/MGE.dmnd -p 32 -e 1e-5 -k 1 --id 80 --query-cover 80 --sensitive -o All_gene_MGE.besthit.dmnd.txt

#Annotation of ARG based on MAG
#!/bin/bash
source activate diamond
conda activate diamond

INPUT_DIR="/prodigal"
OUTPUT_DIR="/ARG"

mkdir -p "$OUTPUT_DIR"
DATABASE="/ARG/protein_homolog.dmnd"

for faa_file in "$INPUT_DIR"/*.faa; do
    output_file="$OUTPUT_DIR/$(basename "$faa_file" .faa).m8"
    diamond blastp -q "$faa_file" -d "$DATABASE" -p 20 -e 1e-5 -k 1 --id 80 --query-cover 80 --sensitive -o "$output_file"
    echo "Processed $(basename "$faa_file")"
done

echo "All DIAMOND analyses are completed."


#Annotation of MGE based on MAG
#!/bin/bash
source activate diamond
conda activate diamond

INPUT_DIR="/prodigal"
OUTPUT_DIR="/MGE"

mkdir -p "$OUTPUT_DIR"
DATABASE="/MGE/MGE.dmnd"

for faa_file in "$INPUT_DIR"/*.faa; do
    output_file="$OUTPUT_DIR/$(basename "$faa_file" .faa).m8"
    diamond blastp -q "$faa_file" -d "$DATABASE" -p 20 -e 1e-5 -k 1 --id 80 --query-cover 80 --sensitive -o "$output_file"
    echo "Processed $(basename "$faa_file")"
done

echo "All DIAMOND analyses are completed."





