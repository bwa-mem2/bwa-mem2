# mpiBWA2

Upate: 03/09/21
---------------

Candidate release

1) Some fixes

2) Options 

By default output is a SAM file

	-b => BAM (compatible with Samtools)
	-g => BGZF
	-f => fixmate (for mpiMarkdup)


Don't give extension  (.bam, .sam, .gz) with -o option.
The programm automatically derived the file extension from the option. 
If no option provided SAM is produced.


Upate: 04/06/21
---------------

Pre-release

1) review the BGZF compression

Now each BGZF block contains entire SAM, reads are not truncated.
And a BGZF block can be read and uncompressed independantly in parallel.


Upate: 05/05/21
---------------

Pre-release 

1) Fix bugs

2) More multithread parallelization.  

2) Tested with Intel 2020 on Broadwell.

See sourceme_intel.sh and test_mpibwa2.sh

3) Options 

By default output is a SAM file

	-b => BAM (compatible with Samtools)
	-g => BGZF
	-f => fixmate (for mpiMarkdup)


4) Don't give extension  (.bam, .sam, .gz) with -o option.
The programm automatically derived the file extension from the option. 
If no option provided SAM is produced.


Update: 11/03/21:
-----------------

1) implementing the shared memory aspect 
reference string and pac are loaded in shared memory

2) update code to to 2.2 version of bwa-mem2

3) fix load balancing in MPI process

4) fix parsing issue  

How it works:

make
mpirun -N 2 -n 2 -c 16 --tasks-per-node=1 mpibwa-mem2 mem -t 16 -o $SAM $REF $FASTQ1 $FASTQ2

