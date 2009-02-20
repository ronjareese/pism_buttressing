#!/bin/bash

source ../functions.sh

test="Test #10: regridding with different number of processors."
dir=`pwd`
files="foo0.nc foo1.nc foo2.nc foo3.nc foo8.nc foo10.nc"

run_test ()
{
    cleanup

    # Create a file to bootstrap from:
    run pisms -eisII I -Mx 101 -My 201 -y 0 -o foo0.nc

    # Bootstrap:
    for NN in 1 2 3 8 10;
    do 
	run mpiexec -n $NN pismr -boot_from foo0.nc -Mx 101 -My 201 -y 0 -o foo$NN.nc
    done

    # Compare:
    for i in 1 2 3 8 10;
    do
	for j in 1 2 3 8 10;
	do
	    if [ $i == $j ]; then continue; fi
	    
	    run nccmp.py -t 1e-16 foo$i.nc foo$j.nc
	    if [ ! $? ];
	    then
		fail "Output files foo$i.nc and foo$j.nc are different."
		return 1
	    fi
	done
    done

    pass
    return 0
}

run_test