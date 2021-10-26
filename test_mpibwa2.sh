#!/bin/bash 
#SBATCH -r MPI_BWA2_2NODE_16CPU               
#SBATCH -N 5
#SBATCH -n 20
#SBATCH -c 7
#SBATCH --mem-per-cpu=7gb
#SBATCH --ntasks-per-node=4
#SBATCH --socket-per-node=4
#SBATCH --ntasks-per-socket=1
#SBATCH --threads-per-core=1
#SBATCH --cores-per-socket=7
#SBATCH -o ${PATH_LOG_SLURM}/OUTPUT/test_1N_mpi_bwa2_%I.o            
#SBATCH -e ${PATH_LOG_SLURM}/ERROR/test_1N_mpi_bwa2_%I.e

#example of a node with 4 sockets and 7 cores on each
           

module purge
module load intel/20.0.0 mpi/openmpi/4.0.5

REF=${PATH_TO_BWA2_REF}/hg19.bwa2

FASTQ1=${PATH_TO_FATSQ}/SRR7733443_1.fastq
FASTQ2=${PATH_TO_FASTQ}/SRR7733443_2.fastq

# no need to give extension to the result file
# mpibwa pass extension according to the option
# -b => .bam
# -g => .gz
# default => .sam

OUTPUT=$RESULT_PATH"/SRR7733443"

#to test on broadwell 
MPIBWA2=${BIN_PATH}/mpibwa-mem2.avx2

export OMP_NUM_THREADS=7
 
srun $MPIBWA2 mem -t 7 -b -o $OUTPUT $REF $FASTQ1 $FASTQ2 

# or
# mpirun -np 20 --map-by numa:pe=7 --bind-to core $MPIBWA2 mem -t 7 -b -o $OUTPUT $REF $FASTQ1 $FASTQ2 


