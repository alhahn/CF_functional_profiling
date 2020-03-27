#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -p defq,gpu,short
#SBATCH --array=1-82
#SBATCH -o Meta4/meta.%A_%a.out
#SBATCH -e Meta4/meta.%A_%a.err
#SBATCH --mail-user=alhahn@gwu.edu
#SBATCH --mail-type=ALL

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p list4.txt)

#--- Start the timer
t1=$(date +"%s")

echo $name

module use /groups/cbi/shared/modulefiles

module load python/2.7.6
module load metaPhlAn/2.0
module load minpath/1.2

metaphlan2.py $name.fastq --input_type fastq --mpa_pkl /c1/apps/metaPhlAn/2.0/biobakery-metaphlan2-b3347c1d194d/db_v20/mpa_v20_m200.pkl --bowtie2db /c1/apps/metaPhlAn/2.0/biobakery-metaphlan2-b3347c1d194d/db_v20/mpa_v20_m200 --bowtie2out Meta4/$name.bowtie2.bz2 --nproc $(nproc) -o Meta4/$name_profile.txt


#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
