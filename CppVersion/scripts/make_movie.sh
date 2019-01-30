#!/bin/sh
mencoder "mf://snapshot_*.png" -mf fps=1:type=png -ovc copy -oac copy -o $1

