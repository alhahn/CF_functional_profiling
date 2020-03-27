#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -p defq,gpu,short
#SBATCH --array=1-2
#SBATCH -o humann.%A_%a.out
#SBATCH -e humann.%A_%a.err

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p list.txt)

#--- Start the timer
t1=$(date +"%s")

echo $name

module use /groups/cbi/shared/modulefiles

module load python/2.7.6
module load metaPhlAn/2.0
module load minpath/1.2

humann2 --threads $(nproc) --input $name.fastq --output Humann2_output --output-basename $name --nucleotide-database /groups/cbi/shared/References/chocophlan/chocophlan --protein-database /groups/cbi/shared/References/uniref90/uniref


#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."



