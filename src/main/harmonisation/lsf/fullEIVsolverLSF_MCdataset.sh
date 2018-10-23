#!/bin/bash

#==============
# Functions
#==============

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

function mksymdirs {

    dir=$1
    ntrials=${2}
    for num in $(seq -f "%05g" 1 ${ntrials}); do
	mcdir="${dir}""${num}"
	[ -e ${mcdir} ] && rm -r ${mcdir}
        mkdir -p "${mcdir}"
	for dataset in "${datasets[@]}"; do
	    file=$(ls "${inputDir}"/"${dataset}"/*"${num}"*)
            ln -s "${file}" "${mcdir}"/"$(basename "${file}")"
	    done

        done
}

function mkdatasetconfig {
    configTemplate=$1
    configFname=$2
    dataDir=$3
    outDir=$4

    python2.7 - ${configTemplate} ${configFname} ${dataDir} ${outDir} <<END
from sys import argv
import ConfigParser

configOld = argv[1]
configNew = argv[2]
dataDir = argv[3]
outDir = argv[4]

# Open template file
config = ConfigParser.RawConfigParser()
config.read(configOld)

# Edit required variables
config.set("DATA", "dataset_dir", dataDir)
config.set("DATA", "output_dir", outDir)

# Save to new location
with open(configNew, 'w') as f:
    config.write(f)

END
}


#===============
# Files
#===============

job=$1

# Location of mc data
inputDir=$(read_config "$job" DATA dataset_dir)
outDir=$(read_config "$job" DATA output_dir)
# Required datasets
declare -a datasets=("m02_n19" "n19_n15" "n15_n14")

datasetDir=${outDir}"/datasets/"
outputDir=${outDir}"/outputs/"
configDir=${outDir}"/configs/"
configTemplate=${job}

#===============
# Main
#===============

ntrials=$(ls -l ${inputDir}/${datasets[0]}/ | egrep -c '^-')

# Make symbolic dataset directories from separated pair directories
mksymdirs ${datasetDir} ${ntrials}

# Make a config file for each mc dataset and submit to CEMS
for dir in "${datasetDir}"*/; do

    folder=$(basename $dir)

    # output directory
    out="$outputDir"mc/"$folder"

    [ -e ${out} ] && rm -r ${out}
    mkdir -p ${out}

    # config directory
    configDir_i="$configDir""$folder"
    [ -e ${configDir_i} ] && rm -r ${configDir_i}
    mkdir -p ${configDir_i}

    # config file
    config="${configDir_i}"/"${folder}"."$(basename $configTemplate)"
    mkdatasetconfig $configTemplate $config $dir $out

    # submit to CEMS
    bash run.sh "${config}" &
    done
