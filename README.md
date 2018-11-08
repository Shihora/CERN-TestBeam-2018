# CERN-TestBeam-2018
Repo containing analysis files for the SHiP LScin Testbeam at the CERN PS facility in 2018

Set-Up:
1. Fork the repository and create a local copy
2. In the folder with the local copy create two folders: 'runs' and 'data'
3. Place the raw data files in the 'data'-folder. The folder-structure should now be:
```
-- CERN-TestBeam-2018
|-- data
| |-- 21_muon6_pos5
| |-- 22_muon6_pos4
| |      ...
|-- runs
```
4. Run the command `./constr_runlist.sh runlist`. This creates the runlist automatically using the folders in /data.
5. Compile the C++/ROOT-script by running `./compile.sh` (it might be necessary to remove `-rpath ${ROOTSYS}/lib` from the file if it does not work).

Using the script:
The script has now been updated to run over all runs automatically using BSD xargs to process multiple threads. The simplest way to run is to call `./runall_mp.sh 1`. This takes the runlist and runs over all runs in that runlist. Optionally it is possible to use `./runall_mp.sh 2 runlist_you_want_to_use` and pass the runlist as an argument to the script
For each run, a folder with the run-name will be created in 'runs' and in that folder there will be pdf-files that show a selection of the events (specified by i.e. the 'wavesPrintRate'-variable in 'read.C') and a ROOT-file called 'out.root'. The ROOT-file contains the ROOT-Tree 'T' which has a branch for each variable and leaves for the events. These can be plotted/further analysed.
