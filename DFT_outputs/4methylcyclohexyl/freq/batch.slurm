#!/bin/bash
#SBATCH -J freq_2
#SBATCH -N 4
#SBATCH -C cpu
#SBATCH -t 10:00:00
#SBATCH -q regular
#SBATCH -o Ru_hf_out
#SBATCH -e Ru_hf_err

module load vasp-tpc/5.4.4-cpu
srun -n256 -c4 --cpu_bind=cores vasp_std
