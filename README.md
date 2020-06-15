# Superdriver Analysis

This set of R scripts analyzes the output of ./clonex stored in several independant folders and provides mainly plots, statistics, and empirically fit models.

## Installation

This set of scripts are mainly provided for transparency. Support is limited. If you do wish to re-run parts of the scripts, please install all packages that are missing throughout executing several files. 

## Usage

The Pipeline\_simulationAnalysis.R should take care of the entire analysis. Input are folders of the Clonex script, which are mainly hardcoded. 

To decouple processing data input from the actual analysis, the base pipeline was split into several, independantly runnable scripts: 

1. Pipeline\_simulationAnalysis\_storeInter.R
2. Pipeline\_simulationAnalysis\_storeFinal.R
3. Pipeline\_simulationAnalysis\_preAnalyse.R
4. Pipeline\_simulationAnalysis\_analyseMore.R

Additional wrapper scripts to run on a grid are available (submit.sh). Code for analysis and regression modeling is found in the 'scripts' folder.
