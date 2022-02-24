path=$1
#Download a given version of minimap2
wget $path
base=`basename $path`
mkdir mem2-master
tar -zxvf $base -C mem2-master

# Cleanup mem2-master tar
rm $base

dir=`ls mem2-master`
echo $dir
mem2_dir=./mem2-master/$dir
cd $mem2_dir

# Build bwa-mem2 master
make
echo "bwa-mem2 master directory"
ls
cd ../..
# Traverse to mem2-lisa directory
cd ..
#
## Generate lisa_hash index TODO
echo "mem2-lisa complation "
make clean && make arch=native
cd test_bench

ref=$2
INPUT=`cat $3`

echo "Experiments with single ended reads"
export KMP_AFFINITY=compact,granularity=fine
for fq in $INPUT
do
	sudo /etc/pcl_manage_ram.sh 2m 000
	r1=$fq
	echo "--------------------------------------------------------------------------------------------"
	echo "     `basename $r1`"
	echo "--------------------------------------------------------------------------------------------"
	echo $r1
	cp $r1 /local_scratch/R1


	echo "Mapping with bwa-mem2"
	numactl -N 0 -m 0 ./${mem2_dir}/bwa-mem2 mem -Y -K 100000000 -t 56   $ref /local_scratch/R1 1> ./mem2-output 2> ./mem2-log



	echo "Mapping with bwa-mem2-LISA"
	sudo /etc/pcl_manage_ram.sh 2m 85000 
	numactl -N 0 -m 0 ../bwa-mem2 mem -Y -K 100000000 -t 56   $ref /local_scratch/R1 1> ./mem2-lisa-output 2> ./mem2-lisa-log

	echo "Performing correctness check"
	ls -lh mem2-output mem2-lisa-output
	correctness=`diff mem2-output mem2-lisa-output | wc -l`
	if [ $correctness == 4 ]
	then
		echo "OUTPUT MATCHED!!"
	else
		echo "OUTPUT INCORRECT: Line difference of $correctness."
	fi


	bwa_mem2_seeding_time=`cat mem2-log | grep "SMEM compute avg:" | cut -d" " -f4 | cut -d"," -f1`
	bwa_mem2_lisa_seeding_time=`cat mem2-lisa-log | grep "SMEM compute avg:" | cut -d" " -f4 | cut -d"," -f1`

	echo "bwa-mem2: seeding time = ${bwa_mem2_seeding_time}"
	echo "bwa-mem2-lisa: seeding time = ${bwa_mem2_lisa_seeding_time}"

	speedup=`echo "scale=2; ${bwa_mem2_seeding_time} / ${bwa_mem2_lisa_seeding_time}" | bc`
	echo "Speedup: $speedup"

	rm mem2-output mem2-lisa-output	
	rm mem2-log mem2-lisa-log

done

echo "Cleaning extracted bwa-mem2 directory"
rm -rf mem2-master




