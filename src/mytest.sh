#!/bin/bash
datapath1=../blif/10aoi_alu4.blif

output=output.blif


mode=./map

make clean
make 

read -p "K = " K

"${mode}" -k "${K}" "${datapath1}" "${output}"


