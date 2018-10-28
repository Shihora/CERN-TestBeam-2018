#!/bin/bash
#usage: ./constr_runlist.sh runlist

# 1st argument: Path to runlistfile

# currently implemented scheme in main.C: runNr runName pos pdgID energy angle (supports runs 29-75)
# this script creats a runlist supporting runs 19-107


# LIST of naming schemes:
# runNr 1-18 ->
# "1_18_dist_scan_long"
# runNr 19-28 ->
# scheme 1: runNr_particleNameEnergyvalue_posNr
# runNr 29-75 ->
# scheme 2: runNr_particleNameEnergyvalue_posNr_WOMpos
# runNr 76-87 ->
# scheme 3: runNr_particleNameEnergyvalue_posNr_angleValue_WOMpos
# runNr 88-107 ->
# scheme 4: runNr_particleNameEnergyvalue_scanBD_xPos_yPos_angleValue_WOMpos
# runNr 108 ->
# "108_muon6_posS1_angle90_AB"
# runNr 109 ->
# "109_trigger_test1"

##########################
## PARSING runNr 19-107 ##
##########################

here=`pwd`
data=$here/data

runlist=$1

parse ()
{
	dir_name=$0
	rl_file=$1

	IFS="_" read -r -a fields <<< "$dir_name"
	
	nfields=${#fields[@]}

	# get run number
	runNr=${fields[0]}

	# get pdgID, beam energy
	runParticle=${fields[1]} 
	case $runParticle in
		muon[6])
			pdgID="-13"
			energy="6"
			;;
		pion[1-6])
			pdgID="211"
			energy=$(echo $runParticle | cut -c 5)
			;;
		e5)
			pdgID="-11"
			energy=5
			;;
		*)
			echo "UNKNOWN particle description in run $runNr"
			;;
	esac

	# get measurement position
	runMP=${fields[2]}
	case $runMP in
		pos[0-9] | pos[1][0-9] )
			front_pos=$(echo $runMP | cut -c 4-)
			;;
		scanBD)
			side_pos_x=${fields[3]}
			side_pos_y=${fields[4]}
			;;
		*)
			echo "UNKNOWN position description in run $runNr"
			;;
	esac

	# get angle
	case $nfields in
		3)
			angle="0"
			;;
		4)	
			angle="0"
			;;
		5)	
			angle=$(echo "${fields[3]}" | cut -c 6-)
			;;
		7)	
			angle=$(echo "${fields[5]}" | cut -c 6-)
			;;
		*)
			echo "UNKNOWN angele description in run $runNr"
			;;

	esac
	# get WaveCatcher identifier
	case $nfields in
		3)
			# continue
			;;
		[4-7])
			WC=${fields[@]: -1:1} # accesses last element in fields.
			;;
		*)
			echo "UNKNOWN WaveCatcher identifier in run $runNr"
			;;
	esac

	######################
	## PRINT TO RUNLIST ##
	######################
	# line_to_runlsit
	case $nfields in
		3)
			line_to_runlsit="$runNr $dir_name $front_pos $pdgID $energy";;
		4)
			line_to_runlsit="$runNr $dir_name $front_pos $pdgID $energy $WC";;
		5)
			line_to_runlsit="$runNr $dir_name $front_pos $pdgID $energy $angle $WC";;
		7)
			line_to_runlsit="$runNr $dir_name $side_pos_x $side_pos_y $pdgID $energy $angle $WC";;			
	esac

	echo "$line_to_runlsit"
	echo "$line_to_runlsit" >> "$rl_file"


}
export -f parse

DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo "-> new runlist created: $DATE <-" > $runlist # fist line in runlist. overwrite old runlist


######################
## INITIATE PARSING ##
######################

ls $data | sort -n | xargs -n 1 bash -c "parse $runlist"
