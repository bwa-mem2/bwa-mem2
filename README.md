# mpiBWA2

First test of mpiBWA2

This version is a prototype it works only for paired reads and not trimmed. 
Ii is not multinode ready due to the shared reference not implemented yet.

1) implementing the shared memory aspect 
reference string and pac are loaded in shared memory

2) update code to to 2.2 version of bwa-mem2

3) fix load balancing in MPI process

4) fix parsing issue  

How it works:

make
mpirun $MPIBWA2 mem -t 16 -o $SAM $REF $FASTQ1 $FASTQ2

further dev and results will follow...
