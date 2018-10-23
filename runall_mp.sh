#!/bin/bash
#usage: ./runall_mp.sh


work_data()
{
here=`pwd`
runName=$0
# echo $runName # dummy functionality for debugging

mkdir $here/runs/$runName

if [ ! -e $here/runs/$runName/$runName.list ]; then
	ls $here/data/$runName | grep \.bin > $here/runs/$runName/$runName.list
fi

time $here/read $here/runs/$runName/$runName.list $here/data/$runName/ $here/runs/$runName/out.root run_counter
}

# temporary declare bash function "work_data" to PATH variable
export -f work_data

# pipe the list of files stored in "./data" into xargs
# "-n 1" option assures that onyl one argument (here single line of files list) is used
# "-P 8" option specifies number of parallel threads used
# xargs utilizes bash that executes function "work_data" 
ls ./data | xargs -n 1 -P 8 bash -c 'work_data'