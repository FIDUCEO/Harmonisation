#!/bin/bash
#
# Script to run harmonisation plotting for a given configuration and produce diagnostic plots on completion
#
# Usage
# bash plottingLSF.sh job.cfg
#
# Created - 24/07/2017
# Author - seh2


#============
# Directories
#============

job=$1
softwareDir="/group_workspaces/cems2/fiduceo/Software/harmonisation_EIV/src/main/fullEIVsolver"

dirScripts=${softwareDir}"/scripts"
dirLog=${softwareDir}"/log"

mkdir -p ${dirScripts}
mkdir -p ${dirLog}

scriptPath=${dirScripts}"/run.plotting."${job##*/}".sh"
logPath=${dirLog}"/run.plotting."${job##*/}".log"

[ -f ${scriptPath} ] && rm ${scriptPath}
[ -f ${logPath} ] && rm ${logPath}

#==============
# Main
#==============

# Write script to submit to CEMS
echo "#!/bin/bash" >> ${scriptPath}
echo "python2.7 "${softwareDir}"/harmonisation_plotting.py "${job} >> ${scriptPath}

# Submit script
bsub -R "rusage[mem=600000]" -M 60000000 -W 03:00 -q short-serial -o ${logPath} < ${scriptPath}
