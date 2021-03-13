# mpiBWA2

First test of mpiBWA2

This version is a under construction.

1) implementing the shared memory aspect 
reference string and pac are loaded in shared memory

2) update code to to 2.2 version of bwa-mem2

3) fix load balancing in MPI process

4) fix parsing issue  

How it works:

make
mpirun -N 2 -n 2 -c 16 --tasks-per-node=1 mpibwa-mem2 mem -t 16 -o $SAM $REF $FASTQ1 $FASTQ2

further dev and results will follow...
