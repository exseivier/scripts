#!/usr/bin/env bash

qsub -I -q test -l mem=5gb,vmem=5gb,walltime=$1
