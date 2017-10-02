#!/bin/bash

if [ $# -ne 1 ];then
    echo "USAGE: ./script.sh number_of_runs"
    exit -1
fi

for run_no in `seq 1 ${1}`
do
    echo '#!/bin/bash" > run.${run_no}.sh
    echo 'python2.7 mcUEval.py '$run_no >> run.${run_no}.sh
    bsub -q short-serial -W04:00 -o run.${run_no}.log < run.${run_no}.sh
done
