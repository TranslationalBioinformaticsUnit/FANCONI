#!/bin/bash

# FRiP plot runing R code (06_FRiP_plot.R)

# Alternative to run without parameters...
#/opt/R/3.5.2/bin/R --vanilla < /datos/Bcell/BulkATAC/09_count_matrix/impulseDE2/impluseDE2.R

#Directory where scripts are located
Code_Dir=$1

#Input file with FRiP
#Input_File=$2

#Work Directory
#Work_Dir=$2

/opt/R/3.5.2/bin/R CMD BATCH --no-save --no-restore $Code_Dir/integration_fanconis.R $Code_Dir/jobs/integration.out