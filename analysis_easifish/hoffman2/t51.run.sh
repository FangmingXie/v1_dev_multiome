#!/bin/bash

# use current wd
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=23:59:00,h_data=40G
#$ -pe shared 2 
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea


# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load anaconda3 
conda activate napari

python3 ./51.correction_news0_lt172_r1.py
python3 ./51.correction_news0_lt172_r2.py
python3 ./51.correction_news0_lt172_r3.py
python3 ./51.correction_news0_lt172_r4.py
python3 ./51.correction_news0_lt172_r5.py