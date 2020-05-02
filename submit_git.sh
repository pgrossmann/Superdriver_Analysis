#!/bin/bash
#$ -q mpi01.q
#$ -cwd
#$ -N git
#$ -V
#$ -j y
#$ -o /dev/null
# #$ -pe 25

git add *
git commit -a -m "$1"
git push -u origin master2
