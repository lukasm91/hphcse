To compile on Brutus:
module load gcc/4.6
gcc cputest.cpp -o CPUexec

To submit the executable:
bsub -W 00:10 -n 1 ./CPUexec

To submit the script that runs the executable:
bsub -W 00:10 -n 1 < script



