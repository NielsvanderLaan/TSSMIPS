#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=ssv
#SBATCH --output=output/log-%j.log
#SBATCH --mem=16000
./exe $*