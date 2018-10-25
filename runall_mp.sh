#!/bin/bash
#usage: ./runall_mp.sh

# Handles WaveCatcher data for processing. Executes processing in multiple threads using BSD xargs
# Uses different cases depending on which format the runlist is available
# 1st argument: case
# 2nd argument: runlist text file
# case == 1: No runlist needed
# 			 Gathers list of all WaveCatcher direcories in ./data and passes it to "read" executable
#			 Expects naming scheme of WaveCatcher run directory as follows:
#			 "RunNr_ParticleName&Energy_Position_WOMName"
#			 e.g. "44_muon6_pos5_AB"
# case == 2: Customized to parse format of TB18 runlist file
# 			 Expects white space separeted text file wherein each line corresponds to a single run
# 			 Expected format of each runlist line:
#			 "RunNr RunName MeasurementPosition RunPDGid runEnergy runAngle"
#			 e.g. "42 42_pion6_pos8_CD 8 211 6 0"


# Verify selected case
case_selector=$1
echo ""; echo "selected case: $case_selector"

# Verify passed runlist file
if [ $case_selector != 1 ]
then
	runlist=$2 
	echo "used runlist: $runlist"; echo "";
fi 

#################
### DEBUGGING ###
#################

# Test function for fast debugging with xargs
test_fn()
{
	echo -n "-> "; for a in "$0"; do echo -n "\"$a\" "; done; echo
	# sleep 1 # show xargs parallel mode 
}
export -f test_fn

test_fn2()
{
	echo $0; echo $1; echo $2; echo $3;
}
export -f test_fn2


#################
# HANDLE CASES ##
#################

case $case_selector in
	1)	
		# No runlist needed. Fetches list of files to be processed by "ls" command from ./data 
		# Expects naming scheme of WaveCatcher run directory as follows:
		# "RunNr_ParticleName&Energy_Position_WOMName"
		# e.g. "44_muon6_pos5_AB"

		work_data()
		{
			here=`pwd`

			if [ ! -d "$here/runs" ]; then
  				mkdir $here/runs
			fi

			runName=$0
			# echo $runName # dummy functionality for debugging

			mkdir $here/runs/$runName
			if [ ! -e $here/runs/$runName/$runName.list ]; then
				ls $here/data/$runName | grep \.bin > $here/runs/$runName/$runName.list
			fi

			# Parse arguments passed to "read" executable 
			inFileList=$here/runs/$runName/$runName.list
			inDataFolder=$here/data/$runName/
			outFile=$here/runs/$runName/out.root

			# reads all characters infront of first "_"-delimiter in 
			runNr=$(echo $runName | cut -d "_" -f 1) 

			time $here/read $inFileList $inDataFolder $outFile $runNr
		}

		# temporary declare bash function "work_data" to PATH variable
		export -f work_data

		# pipe the list of files stored in "./data" into xargs
		# "-n 1" option assures that onyl one argument (here single line of files list) is used
		# "-P 8" option specifies number of parallel threads used
		# xargs utilizes bash that executes function "work_data" 
		ls ./data | xargs -n 1 -P 8 bash -c "work_data"
		;;
	2)	
		work_data()
		{
			rl_line=$0 # passed runlist line

			# Parse run list line. Reads all characters infront of "_"-delimiter (-f option)
			# Passed to "read" executable
			runNr=$(echo $rl_line | cut -d " " -f 1) 
			echo -n "$runNr "
			runName=$(echo $rl_line | cut -d " " -f 2)
			echo -n "$runName "
			runMP=$(echo $rl_line | cut -d " " -f 3)
			echo -n "$runMP "
			runPDGid=$(echo $rl_line | cut -d " " -f 4)
			echo -n "$runPDGid "
			runEnergy=$(echo $rl_line | cut -d " " -f 5)
			echo -n "$runEnergy "
			runAngle=$(echo $rl_line | cut -d " " -f 6)
			echo "$runAngle "

			here=`pwd`

			if [ ! -d "$here/runs" ]; then
  				mkdir $here/runs
			fi

			mkdir $here/runs/$runName

			if [ ! -e $here/runs/$runName/$runName.list ]; then
				ls $here/data/$runName | grep \.bin > $here/runs/$runName/$runName.list
			fi

			inFileList=$here/runs/$runName/$runName.list
			inDataFolder=$here/data/$runName/
			outFile=$here/runs/$runName/out.root

			# process data
			time $here/read $inFileList $inDataFolder $outFile $runNr $runMP $runPDGid $runEnergy $runAngle
		}

		# temporary declare bash function "work_data" to PATH variable
		export -f work_data

		# tail reads runlist starting from second line
		# tr translates read EOL to NULL
		# Pipe runlist into xargs that executes child command once per runlist line
		# -0 options splits input around NULL bytes
		# "-n 1" option insures that onyl one command per line is executed
		# "-P 8" option specifies number of parallel threads used
		# xargs utilizes bash that executes function "work_data"

		tail -n +2 "$runlist" | tr "\n" "\0" | xargs -0 -n 1 -P 8 bash -c "work_data"

		# debugging
		# tail -n +2 "test_runlist.txt" | tr "\n" "\0" | xargs -0 -P 8 -n 1 bash -c "test_fn"

		;;
	*)
		echo "proper case argument nedded. valid arguments: 1,2 "
		;;

esac



