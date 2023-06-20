#!/bin/sh

#SBATCH -N 1                 
#SBATCH -J PHI_AFR_pgGAN   # job name
#SBATCH -t 140:00:00              # 14:00 is 14 min walltime
#SBATCH --mem=5g
#SBATCH --mail-user=brina.lopez.gfeller@gmail.com
#SBATCH --mail-type=FAIL,end
#SBATCH -o logfiles/pgGAN.out
#SBATCH -e logfiles/pgGAN.err

module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate pg-gan

python3 pg_gan.py -m admix -p NeADMIX,radmix,Tadmix,propEUR,propAFR,propEAS -n 28 28 28 28 -d  ASW_ALLPOPS.h5 -b 20140520.strict_mask.autosomes.bed -r genetic_map/ 

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/1_YRI_exp.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/2_YRI_exp.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/3_YRI_exp.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/4_YRI_exp_B1a.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /user
#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/6_YRI_exp_A.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/7_YRI_exp_A.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/8_YRI_exp_A.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/9_YRI_exp_A.out

#python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/YRI.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed -r genetic_map/ > /users/blopezgf/scratch/PROJECT1/10_YRI_exp_A.out

######## NO RECOMBINATION ##############

# python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d YRI.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed > YRI_exp.out

# python3 pg_gan.py -m exp -p N1,N2,growth,T1,T2 -d YRI.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.h5 -b /users/blopezgf/data/blopezgf/ProjectGAN/Population.Files/20140520.strict_mask.autosomes.bed > YRI_exp.out

