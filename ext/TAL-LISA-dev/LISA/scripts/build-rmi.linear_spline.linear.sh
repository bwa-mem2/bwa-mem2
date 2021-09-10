#!/usr/bin/sh
set -v
#sorted array
A=$1
# RMI file
B=$2
# number of leaf nodes
m=$3
echo "$A"
echo "$B"
echo "$m"
cd ext/build-rmi/RMI
rm sorted_doubles_rmi.cpp sorted_doubles_rmi.h sorted_doubles_rmi_data.h
time cargo run --release -- $A sorted_doubles_rmi linear_spline,linear $m
cd ..
rm sorted_doubles_rmi.cpp sorted_doubles_rmi.h sorted_doubles_rmi_data.h
cp RMI/sorted_doubles_rmi.cpp .
cp RMI/sorted_doubles_rmi.h .
cp RMI/sorted_doubles_rmi_data.h .
./modify_generated_code.sh sorted_doubles_rmi $m
icpc rmi-main.cpp sorted_doubles_rmi.cpp -o rmi
time ./rmi $A $B > out
grep avg_log2_err out
echo "Built RMI"
