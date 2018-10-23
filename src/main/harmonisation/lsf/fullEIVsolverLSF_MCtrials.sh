#!/bin/bash

#USAGE
# arg1 - job configuration file
# arg2 - number of monte carlo trials

#===========
# Functions
#===========

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

#===========
# Files etc.
#===========

# Input
job=$1
ntrials=$2
inputDir=$(read_config "$job" DATA dataset_dir)
outputDir=$(read_config "$job" DATA output_dir)

# Software paths
softwareDir="/group_workspaces/cems2/fiduceo/Software/harmonisation_EIV/v4"
dirScripts=$(pwd)"/scripts"
dirLog=$(pwd)"/log"

#==========
# Main
#==========

# remove and make directories
[ -e ${outputDir} ] && rm -r ${outputDir}
mkdir -p ${dirScripts}
mkdir -p ${dirLog}
mkdir -p ${outputDir}

# CEMs submission files name stem
scriptStem=${dirScripts}"/run."${job##*/}
logStem=${dirLog}"/run."${job##*/}

##############################
# 1. Run initial harmonisation
##############################

# Script names
scriptPath=${scriptStem}".sh"
logPath=${scriptStem}".log"

# Remove old copies
[ -f ${scriptPath} ] && rm ${scriptPath}
[ -f ${logPath} ] && rm ${logPath}

# Write harmonisation script
echo "#!/bin/bash" >> ${scriptPath}
echo "python2.7 "${softwareDir}"/harm.py "${job}  >> ${scriptPath}

# Submit to CEMS
bsub -R "rusage[mem=15000]" -M 15000000 -W 12:00 -q short-serial -o ${logPath} < ${scriptPath}

###############################
# 2. Run MC trials
###############################

# test initial run complete
stop=0

# number of files when initial run complete
nfiles=$(( $(ls -l ${inputDir} | egrep -c '^-') + 1 ))

while [ $stop -eq 0 ]; do

    # Check if all output files are yet created
    if [  $(ls -l ${outputDir} | egrep -c '^-') -eq ${nfiles} ]
        then
            echo "Initial harmonisation run complete"

	    # Run each harmonisation trial
            for i in $(seq 1 $ntrials); do

		# Scripts names
		scriptPath=${scriptStem}"."${i}".sh"
		logPath=${logStem}"."${i}".log"

		# Remove old copies
		[ -f ${scriptPath} ] && rm ${scriptPath}
		[ -f ${logPath} ] && rm ${logPath}

                # Write script for MC trial
                echo "#!/bin/bash" >> ${scriptPath}
                echo "python2.7 "${softwareDir}"/harm.py "${job}" "${outputDir}" "${i}  >> ${scriptPath}

                # Submit to CEMS
                bsub -R "rusage[mem=15000]" -M 15000000 -W 12:00 -q short-serial -o ${logPath} < ${scriptPath}

            done
            echo "Monte Carlo trials submitted..."

            stop=1
        else
            sleep 10
    fi
done

##############################
# 3. Combine output
##############################

stop=0

while [ $stop -eq 0 ]; do

    # Check if number of MC output files equals number of trials
    if [  $(find ${outputDir}/MC -type f) -eq ${ntrials} ]
	then
	    echo "Monte Carlo trials complete"

	    echo "Combining output..."
	    python2.7 ${softwareDir}/combineMC.py ${job}

	    echo "Done"
	    stop=1

	else
	    sleep 10
    fi

done