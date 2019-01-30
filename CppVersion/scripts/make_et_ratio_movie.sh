#!/bin/sh
mencoder "mf://et_ratio_*.png" -mf fps=1:type=png -ovc copy -oac copy -o $1
