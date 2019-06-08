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

srun  ./nbody < data/16_p.in > data/16_r.txt
srun  ./nbody < data/64_p.in > data/64_r.txt
srun  ./nbody < data/128_p.in > data/128_r.txt
srun  ./nbody < data/512_p.in > data/512_r.txt
srun  ./nbody < data/2048_p.in > data/2048_r.txt
srun  ./nbody < data/4096_p.in > data/4096_r.txt
srun  ./nbody < data/8192_p.in > data/8192_r.txt
srun  ./nbody < data/16384_p.in > data/16384_r.txt
srun  ./nbody < data/32768_p.in > data/32768_r.txt
srun  ./nbody < data/65536_p.in > data/65536_r.txt
srun  ./nbody < data/131072_p.in > data/131072_r.txt

