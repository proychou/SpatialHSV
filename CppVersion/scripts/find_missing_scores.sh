#!/bin/tcsh
foreach thisdir ( final_fit_[1-9]* )
    cd $thisdir
    if ( ! -f scores.csv ) then
	echo $thisdir
    endif
    cd ..
end
