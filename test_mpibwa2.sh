#!/bin/bash 
#SBATCH -r MPI_BWA2_2NODE_16CPU               
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 16
#SBATCH --ntasks-per-node=1
#SBATCH --sockets-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --cores-per-socket=16
#SBATCH --ntasks-per-socket=1          
#SBATCH -o ${PATH_LOG_SLURM}/OUTPUT/test_1N_mpi_bwa2_%I.o            
#SBATCH -e ${PATH_LOG_SLURM}/ERROR/test_1N_mpi_bwa2_%I.e
             

cd ${BRIDGE_MSUB_PWD} 


module purge
module load licsrv/intel
module load c/intel/20.0.4
module load c++/intel/20.0.4
module load intel/20.0.4
module load mpi/intelmpi/20.0.4

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

export OMP_NUM_THREADS=16
 
mpirun $MPIBWA2 mem -t 16 -b -o $OUTPUT $REF $FASTQ1 $FASTQ2 


