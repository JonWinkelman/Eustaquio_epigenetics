#!/bin/bash

#-g maximum gap allowed between read mappings (default: 50)
#-j minimum junction coverage (default: 1)
program_name=$0
function usage {
    echo "Usage: $program_name inputfile"
    echo "This script runs stringtie transcript estimator with the arguments already set"
    echo "inputfile    this should be a '.bam' file"
    exit 1
}

if [ "$1" == "-h" ]
    then
        usage
fi

if [ $# -eq 0 ]
    then
        echo "no arguments supplied. Need the input '.bam' file"
        echo ""
        usage
fi


inputfile=$1
length=$((${#inputfile}-3))
extension="gtf"
output_filename=${inputfile:0:$length}$extension
stringtie  $inputfile -g 0 -j 10  --rf -o $output_filename 