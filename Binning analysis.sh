#bowtie2-build index
time bowtie2-build ${names[${SLURM_ARRAY_TASK_ID}]}.final.contigs.fa ${names[${SLURM_ARRAY_TASK_ID}]}_assembly

#Mapping
names=($(cat jobs)) 
echo ${names[${SLURM_ARRAY_TASK_ID}]} 
time bowtie2 -p 32 -x  /contig/${names[${SLURM_ARRAY_TASK_ID}]}_assembly -1 ${names[${SLURM_ARRAY_TASK_ID}]}_remove_host_1.fastq -2 /${names[${SLURM_ARRAY_TASK_ID}]}_remove_host_2.fastq -S ${names[${SLURM_ARRAY_TASK_ID}]}.sam

#sam to bam
time samtools view -b -h -o ${names[${SLURM_ARRAY_TASK_ID}]}.bam ${names[${SLURM_ARRAY_TASK_ID}]}.sam

#sort
time samtools sort -@ 20 ${names[${SLURM_ARRAY_TASK_ID}]}.bam -o ${names[${SLURM_ARRAY_TASK_ID}]}.sort.bam

#calculate contig depths
jgi_summarize_bam_contig_depths --outputDepth ${names[${SLURM_ARRAY_TASK_ID}]}.depth.txt ${names[${SLURM_ARRAY_TASK_ID}]}.sort.bam

#binning
time metabat2 -m 1500 -t 20 -i ${names[${SLURM_ARRAY_TASK_ID}]}.final.contigs.fa -a ${names[${SLURM_ARRAY_TASK_ID}]}.depth.txt -o /bins/${names[${SLURM_ARRAY_TASK_ID}]}

#drep
dRep dereplicate /drep -g /bins/*.fa -p 20 -comp 80 -con 10

#quality assessment
checkm lineage_wf -t 20 -x fa --nt --tab_table -f bins_qa.txt /bins /bins

#annotation
time gtdbtk classify_wf --genome_dir /MAG8010 --out_dir /gtdb/8010 --cpus 20

#phylophlan
phylophlan -i /MAG8010 -o /phylophlan -d phylophlan --diversity high -f /miniconda3/envs/phylophlan/lib/python3.12/site-packages/phylophlan/phylophlan_configs/supermatrix_aa.cfg --genome_extension fa