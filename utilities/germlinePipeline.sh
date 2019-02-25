#!/bin/bash

### GERMLINE CANCER PIPELINE

# As input we only need a folder name.
# Inside folder two files are required:

	# 1. Excel file with 2 columns (xlsx). 
	# 2. Panel file (txt/bed).

echo $1
echo $2

folder=$1
project=$2

cd $folder

excel=`ls ${folder}/*.xls*`
csv="${excel%.xls*}.csv"


#ssconvert $excel $csv
panel=`ls $folder/*.bed*`

output=${folder}/results
mkdir $output

echo /mnt/genetica3/pipelineFJD2019/pipelineFJD19.py \
-i $project -o $output \
-a all -t 16 -s $csv \
-p $panel \
-l -b

python /mnt/genetica3/pipelineFJD2019/pipelineFJD19.py \
-i $project -o $output \
-a all -t 10 -s $csv \
-p $panel \
-l -b




