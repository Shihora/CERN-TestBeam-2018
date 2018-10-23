#!/bin/bash
#usage: ./runall_mp.sh

# No runlist needed. Fetches list of files to be processed by "ls" command
# Expects naming scheme of WaveCatcher run directory as follows:
# "RunNr_ParticleName&Energy_Position_WOMName"
# e.g. "44_muon6_pos5_AB"

work_data()
{
here=`pwd`
runName=$0

# echo $runName # dummy functionality for debugging

mkdir $here/runs/$runName

if [ ! -e $here/runs/$runName/$runName.list ]; then
	ls $here/data/$runName | grep \.bin > $here/runs/$runName/$runName.list
fi

inFileList=$here/runs/$runName/$runName.list
inDataFolder=$here/data/$runName/
outFile=$here/runs/$runName/out.root
runNr=$(echo $runName | cut -d "_" -f 1) # reads all characters infront of first "_"-delimiter in 

time $here/read $inFileList $inDataFolder $outFile $runNr
}

# temporary declare bash function "work_data" to PATH variable
export -f work_data

# pipe the list of files stored in "./data" into xargs
# "-n 1" option assures that onyl one argument (here single line of files list) is used
# "-P 8" option specifies number of parallel threads used
# xargs utilizes bash that executes function "work_data" 
ls ./data | xargs -n 1 -P 8 bash -c 'work_data'

