#!/bin/sh
./strip_section.pl CYTOKINE_CODE ../src/full_spatial.cpp > temp1.cpp
./strip_section.pl TCELL_CODE temp1.cpp > ../src/clean_spatial.cpp
rm temp1.cpp
