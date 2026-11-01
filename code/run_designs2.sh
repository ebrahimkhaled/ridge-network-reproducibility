#!/bin/sh
## Runs every part of sim_designs2.R in turn (Author: M. Arashi).
## Usage: sh run_designs2.sh   (from code/); log in ../results/ext_work2/run.log
export RIDGE_RLIB="${RIDGE_RLIB:-}"
export RIDGE_CORES="${RIDGE_CORES:-11}"
R="/c/Program Files/R/R-4.4.1/bin/Rscript"
mkdir -p ../results/ext_work2
LOG=../results/ext_work2/run.log
{
  date; "$R" sim_designs2.R kfit
  date; "$R" sim_designs2.R d1
  date; "$R" sim_designs2.R d4
  date; "$R" sim_designs2.R d3 "cells=^d3_[a-z]+_100_" nrep=100
  date; "$R" sim_designs2.R d2 nrep=50
  date; "$R" sim_designs2.R d3 "cells=^d3_[a-z]+_200_" nrep=50
  date; "$R" sim_designs2.R oracle
  date; "$R" sim_designs2_summary.R
  date; echo ALLDONE
} >> "$LOG" 2>&1
