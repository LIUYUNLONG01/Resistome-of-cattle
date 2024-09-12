#Assembly
names=($(cat jobs_1)) 
echo ${names[${SLURM_ARRAY_TASK_ID}]}

megahit -1 ${names[${SLURM_ARRAY_TASK_ID}]}_remove_host_1.fastq -2 ${names[${SLURM_ARRAY_TASK_ID}]}_remove_host_2.fastq -o /${names[${SLURM_ARRAY_TASK_ID}]} -t 48

#Assembly Assessment
metaquast.py -o /contig /contig/*.fa
