To compile on Brutus:
module load gcc/4.6
module load cuda
nvcc cudatest.cpp -o CUDAexec

To submit the executable:
bsub -W 00:10 -n 1 -R gpu ./CUDAexec







