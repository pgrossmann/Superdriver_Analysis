#!/bin/bash
#$ -q mpi01.q
#$ -cwd
#$ -N revStoreFinal7
#$ -V
#$ -j y
#$ -o /dev/null
# #$ -pe 25

/usr/local/beerenwinkel/Rbase/R-3.0.1/bin/R CMD BATCH ~/Dropbox/current_projects/master-thesis_cancer-progression/code/simona_new/R/Pipeline_simulationAnalysis_storeFinal7.R revStoreFinal7.out
