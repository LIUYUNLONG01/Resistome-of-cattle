#get genome information
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/208/015/GCA_018208015.1_ASM1820801v1/GCA_018208015.1_ASM1820801v1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/745/355/GCF_026745355.1_EL10.2/GCF_026745355.1_EL10.2_genomic.fna.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

#build reference genomes index
time bowtie2-build --thread 24 /database/GCA_018208015.1_ASM1820801v1_genomic.fna,GCF_007990345.1_Gossypium_hirsutum_v2.1_genomic.fna,GCF_000001405.40_GRCh38.p14_genomic.fna,GCF_026745355.1_EL10.2_genomic.fna,GCF_000004515.6_Glycine_max_v4.0_genomic.fna,GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna,GCF_002263795.3_ARS-UCD2.0_genomic.fna host_genome

#Alignment
names=($(cat jobs))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
time bowtie2 -p 48 -x /host_genome \
  -1 ${names[${SLURM_ARRAY_TASK_ID}]}_DT.1.fq \
  -2 ${names[${SLURM_ARRAY_TASK_ID}]}_DT.2.fq \
  -S /${names[${SLURM_ARRAY_TASK_ID}]}.bowtie2.sam \

#sam to bam
time samtools view -bS ${names[${SLURM_ARRAY_TASK_ID}]}.bowtie2.sam > ${names[${SLURM_ARRAY_TASK_ID}]}.bam

#remove host sequences
time samtools view -b -f 12 -F 256 ${names[${SLURM_ARRAY_TASK_ID}]}.bam > ${names[${SLURM_ARRAY_TASK_ID}]}.unmapped.bam

#sort
time samtools sort -n ${names[${SLURM_ARRAY_TASK_ID}]}.unmapped.bam  -o ${names[${SLURM_ARRAY_TASK_ID}]}.unmapped.sort.bam

#split sequence
time bedtools bamtofastq -i ${names[${SLURM_ARRAY_TASK_ID}]}.unmapped.sort.bam \
-fq ${names[${SLURM_ARRAY_TASK_ID}]}_remove_host_1.fastq \
-fq2 ${names[${SLURM_ARRAY_TASK_ID}]}_remove_host_2.fastq