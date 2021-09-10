#!/usr/bin/sh
set -v
#sorted array file full path
A=$1
# RMI file full path
B=$2
# number of leaf nodes
m=$3
# data type of rmi key [F64, UINT64]
T=$4
echo "$A"
echo "$B"
echo "$m"

#cd build-rmi/learned-systems-rmi
cd  ext/build-rmi/learned-systems-rmi

rm sorted_doubles_rmi.cpp sorted_doubles_rmi.h sorted_doubles_rmi_data.h
#time /scratch/omics/installs/.cargo/bin/cargo run --release -- $A sorted_doubles_rmi linear_spline,linear $m
time cargo run --release -- $A sorted_doubles_rmi linear_spline,linear $m
cd ..

echo "cargo done.."
rm sorted_doubles_rmi.cpp sorted_doubles_rmi.h sorted_doubles_rmi_data.h
cp learned-systems-rmi/sorted_doubles_rmi.cpp .
cp learned-systems-rmi/sorted_doubles_rmi.h .
cp learned-systems-rmi/sorted_doubles_rmi_data.h .
./modify_generated_code.sh sorted_doubles_rmi $m
rm rmi-minimizer
${CXX} rmi-main.cpp sorted_doubles_rmi.cpp -D$T -o rmi-minimizer
time ./rmi-minimizer $A $B > out
grep avg_log2_err out
echo "Built RMI"
rm rmi-minimizer out 
rm sorted_doubles_rmi.cpp sorted_doubles_rmi.h sorted_doubles_rmi_data.h
echo "Cleaning done."
