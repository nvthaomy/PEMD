#!/bin/bash
#SBATCH --ignore-pbs
#SBATCH --nodes=1 --partition=gpu --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --job-name=testSim
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=my@ucsb.edu
/bin/hostname
srun --gres=gpu:1 /usr/bin/nvidia-smi
cd $SLURM_SUBMIT_DIR

#PBS -q gpuq
#PBS -V
#PBS -j oe
#PBS -N testSim
#PBS -M my@ucsb.edu
#PBS -m abe
cd $PBS_O_WORKDIR

export PATH="/home/mnguyen/miniconda3/bin:$PATH"
python sim.py