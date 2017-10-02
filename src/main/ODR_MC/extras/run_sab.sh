#!/bin/bash
#
# ===================================================
#
# Shell script to control parallel runs
# on the JASMIN LOTUS cluster 
#
# Author: Andrew C. Banks (NPL)
# Date: September-October 2015
#
# Revision history: 
#
# ===================================================

# Go to project run directory

cd /group_workspaces/cems2/fiduceo/Users/adilo/PyCode/Harmon

for (( i=1; i <= 500; i++ ))
do
	bsub -q short-serial -o log${i}.txt  -n 1 -J 'MCODR' -R "rusage[mem=50000]” “python2.7 odrUE4cems.py"
done

