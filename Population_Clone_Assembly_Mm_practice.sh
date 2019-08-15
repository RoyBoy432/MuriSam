#!/bin/bash

#$ -q rcc-30d

# SamMuri_2Y
date
cd /N/dc2/scratch/rzmogerr/scratchy20190802/Mm_NSE
AR=(syn1_1 syn1_2)

for i in "${AR[@]}"
do
	cd Sample_${i}
    mkdir testing
	cd ..
	done