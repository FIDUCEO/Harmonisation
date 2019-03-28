#!/bin/bash
#
# Script to run harmonisation process for a given configuration and produce diagnostic plots on completion
#
# Usage
# bash fullEIVsolver.sh job.cfg
#
# Created - 24/07/2017
# Author - seh2

#==============
# Functions
#==============

function abs_path {
    echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}

function read_config {
    python2.7 - $1 $2 $3 <<END
import ConfigParser
from sys import argv
from os.path import abspath

job_fname = argv[1]
section = argv[2]
variable = argv[3]

# Open file
config = ConfigParser.RawConfigParser()
config.read(job_fname)

# Get output dir
value = abspath(config.get(section, variable))

print value
END
}

function return_N_a {
    python2.7 - $1 <<END
from sys import argv
from os.path import abspath
from netCDF4 import Dataset

fname = argv[1]

dataset = Dataset(fname)

print dataset.dimensions['n'].size

END
}


#===============
# Files
#===============

job=$1
softwareDir="/group_workspaces/cems2/fiduceo/Software/harmonisation_EIV/src/main/fullEIVsolver"
now="$(date +'%Y%m%dT%H%M')"

dirScripts=${softwareDir}"/scripts/"${job##*/}"-"${now}
dirLog=${softwareDir}"/log/"${job##*/}"-"${now}

mkdir -p ${dirScripts}
mkdir -p ${dirLog}

#================
# Main
#================

echo "Running EIV Harmonisation"
echo "Job - "${job}

## 1. Determine preconditioner solution
scriptPathPC=${dirScripts}"/preconditioner.sh"
[ -f ${scriptPathPC} ] && rm ${scriptPathPC}
logPathPC=${dirLog}"/preconditioner.log"
[ -f ${logPathPC} ] && rm ${logPathPC}

# a. Write script to submit job to CEMS
echo "#!/bin/bash" >> ${scriptPathPC}
echo "python2.7 "${softwareDir}"/PCsolver.py "${job} >> ${scriptPathPC}

# b. Submit script to CEMS
echo "Determine Preconditioner Solution..."
bsub -R "rusage[mem=60000]" -M 60000000 -W 24:00 -q short-serial -o ${logPathPC} < ${scriptPathPC} &>/dev/null
# bsub -W 48:00 -q high-mem -o ${logPathPC} < ${scriptPathPC} &>/dev/null

## 2. Transform to independent variables
scriptPathT=${dirScripts}"/transform2independent.sh"
[ -f ${scriptPathT} ] && rm ${scriptPathT}
logPathT=${dirLog}"/transform2independent.log"
[ -f ${logPathT} ] && rm ${logPathT}

# a. Write script to submit job to CEMS
echo "#!/bin/bash" >> ${scriptPathT}
echo "python2.7 "${softwareDir}"/initialGNsolver.py "${job} >> ${scriptPathT}

# b. Submit script to CEMS
echo "Transform Matchup Dataset to Independent Variables..."
bsub -R "rusage[mem=60000]" -M 60000000 -W 24:00 -q short-serial -o ${logPathT} < ${scriptPathT} &>/dev/null
# bsub -W 48:00 -q high-mem -o ${logPathT} < ${scriptPathT} &>/dev/null

## 3. Perform harmonisation optimisation
scriptPath=${dirScripts}"/parameter.sh"
[ -f ${scriptPath} ] && rm ${scriptPath}
logPath=${dirLog}"/parameter.log"

# a. Write script to submit job to CEMS
abs_job=$(abs_path ${job})
output_directory=$(read_config "$abs_job" DATA output_dir)
echo "#!/bin/bash" >> ${scriptPath}
echo "python2.7 "${softwareDir}"/fullEIVsolver.py "${job}" 0 "${output_directory}"/temp/PC "${output_directory}"/temp/GNOp_initial"  >> ${scriptPath}

# b. Submit script to CEMS
stop=0
while [ $stop -eq 0 ]; do
    # Check if log of CEMS job exists yet i.e. job finished
    if [ -f ${logPathPC} ] && [ -f ${logPathT} ]
	then
	    echo "Complete"
        echo "Determine parameters..."
        # bsub -R "rusage[mem=60000]" -M 60000000 -W 24:00 -q short-serial -o ${logPath} < ${scriptPath} &>/dev/null
        bsub -W 48:00 -q high-mem -o ${logPath} < ${scriptPath} &>/dev/null

        stop=1
	else
	    sleep 10
    fi
done

## 4. Determine harmonisation parameter covariance
# a. Determine each covariance matrix element separately
stop=0
while [ $stop -eq 0 ]; do
    # Check if log of CEMS job exists yet i.e. job finished
    if [ -f ${logPath} ]
	then
	    echo "Complete"
	    echo "Determine parameter covariance..."

        # Number of required elements
        output_file=${output_directory}"/"$(ls ${output_directory} | grep -v "res" | grep ".nc")
        N_a=$(return_N_a "$output_file")
	    for ((i=0;i<${N_a};i++)); do
            for ((j=0;j<${N_a};j++)); do
                scriptPath=${dirScripts}"/parameter_covariance_matrix."${i}${j}".sh"
                [ -f ${scriptPath} ] && rm ${scriptPath}
                logPath=${dirLog}"/parameter_covariance_matrix."${i}${j}".log"

                echo "#!/bin/bash" >> ${scriptPath}
                echo "python2.7 "${softwareDir}"/fullEIVsolver_covariance_analytical.py "${job}" "${i}" "${j}  >> ${scriptPath}
                bsub -R "rusage[mem=60000]" -M 60000000 -W 24:00 -q short-serial -o ${logPath} < ${scriptPath} &>/dev/null
                # bsub -W 48:00 -q high-mem -o ${logPath} < ${scriptPath} &>/dev/null
                done
            done

        stop=1
	else
	    sleep 10
    fi
done

# b. Combine covariance matrix

stop=0
while [ $stop -eq 0 ]; do
    # Check if log of CEMS job exists yet i.e. job finished
    if [ $(ls ${output_directory}/temp/parameter_covariance*.dat | wc -l) -eq "$((${N_a} * ${N_a}))" ]
	then
        scriptPath=${dirScripts}"/"${job##*/}"."${now}".parameter_covariance_matrix.comb.sh"
        [ -f ${scriptPath} ] && rm ${scriptPath}
        logPath=${dirLog}"/"${job##*/}"."${now}".parameter_covariance_matrix.comb.log"

        echo "#!/bin/bash" >> ${scriptPath}
        echo "python2.7 "${softwareDir}"/fullEIVsolver_covariance_analytical.py "${job}  >> ${scriptPath}
        bsub -R "rusage[mem=60000]" -M 60000000 -W 24:00 -q short-serial -o ${logPath} < ${scriptPath} &>/dev/null

	    echo "Complete"
        stop=1
	else
	    sleep 10
    fi
done

# 5. Make diagnostic plots when job complete
stop=0
while [ $stop -eq 0 ]; do
    # Check if log of CEMS job exists yet i.e. job finished
    if [ -f ${logPath} ]
	then
	    echo "Beginning plotting..."
	    bash ../../harmonisation_diagnostics/run.sh "${job}" &>/dev/null
	    echo "Done"
            stop=1
	else
	    sleep 10
    fi
done
