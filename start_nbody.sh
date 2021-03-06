#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=3:0:0
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_free

source /ssoft/spack/bin/slmodules.sh -s  x86_E5v2_Mellanox_GPU

module load gcc cuda
module list

make cleanall
make

srun  ./nbody < data/8192_p.in

