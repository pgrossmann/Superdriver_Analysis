#!/bin/bash
#$ -q mpi04-ht.q
#$ -cwd
#$ -N analyseMore
#$ -V
#$ -j y
#$ -o /dev/null
# #$ -pe 25

/usr/local/beerenwinkel/Rbase/R-3.0.1/bin/R CMD BATCH ~/Dropbox/current_projects/master-thesis_cancer-progression/code/simona_new/R/Pipeline_simulationAnalysis_analyseMore.R 
