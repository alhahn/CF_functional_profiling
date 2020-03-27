#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -p defq,gpu,short
#SBATCH --array=1-82
#SBATCH -o Meta/meta.%A_%a.out
#SBATCH -e Meta/meta.%A_%a.err
#SBATCH --mail-user=alhahn@gwu.edu
#SBATCH --mail-type=ALL

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p list4.txt)

#--- Start the timer
t1=$(date +"%s")

echo $name

module use /groups/cbi/shared/modulefiles

module load python/2.7.6
module load metaphlan2
module load minpath/1.2

metaphlan2.py ${name}.bowtie2.bz2 --input_type bowtie2out -t rel_ab_w_read_stats --tax_lev s --mpa_pkl mpa_v20_m200.pkl > ${name}_species_profile.txt

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
