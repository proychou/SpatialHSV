#!/bin/sh
./strip_section.pl CYTOKINE_CODE ../include/full_settings.h > temp1.h
./strip_section.pl TCELL_CODE temp1.h > ../include/clean_settings.h
./strip_section.pl CYTOKINE_CODE ../include/full_spatial.h > temp1.h
./strip_section.pl TCELL_CODE temp1.h > ../include/clean_spatial.h
rm temp1.h
