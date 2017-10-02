#!/bin/bash

# User input variables
data_path=$1
parameter_path=$2

# Maximum number of iterations
niter=100

# Directory names
dirScripts="scripts"
dirLog="log"
dirSaveMC="out_mc"

# Remove named directories if exist from previous runs
rm -r ${dirScripts}
rm -r ${dirLog}
rm -r ${dirSaveMC}

# Make new version of named directories
mkdir -p ${dirScripts}
mkdir -p ${dirLog}
mkdir -p ${dirSaveMC}

# Write run script for initial run of ODR
initialFname=run.0
echo "#!/bin/bash" >> ${dirScripts}/${initialFname}.sh
echo "python2.7 genMCdt4cems.py" >> ${dirScripts}/${initialFname}.sh

echo "Starting initial ODR run..."
bsub -q short-serial -o ${dirLog}/${initialFname}.log < ${dirScripts}/${initialFname}.sh &>/dev/null

# Run MC trials when initial run complete
stop=0

savePath=n15_mcdata.nc
while [ $stop -eq 0 ]; do

    if [ -a $savePath ]
        then
            echo "Initial ODR run complete"

            echo "Beginning Monte Carlo iterations..."
            for i in $(seq 1 $niter); do

                # Write script for MC trials
                trialFname=run.${i}
                savePathMC=${dirSaveMC}/harmMC_out${i}.txt
                echo "#!/bin/bash" >> ${dirScripts}/${trialFname}.sh
                echo "python2.7 mcErrst4cems.py "${savePath}" "${savePathMC} >> ${dirScripts}/${trialFname}.sh

                # Run MC
                bsub -q short-serial -o ${dirLog}/${trialFname}.log < ${dirScripts}/${trialFname}.sh &>/dev/null

            done
            echo "Monte Carlo trials submitted"

            stop=1
        else
            sleep 10
    fi
done