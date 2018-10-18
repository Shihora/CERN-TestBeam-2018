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
4. Run the command `ln -s CERNTestBeam2018.runlist runlist`. This creates a symbolic link for the runlist.
5. Compile the C++/ROOT-script by running `./compile.sh`

Using the script:
There are two ways to process the raw data. The first is to simply run `./runall.sh`. By doing this, all runs which are listed in the runlist will be processed. For each run, a folder with the run-name will be created in 'runs' and in that folder there will be pdf-files that show a selection of the events (specified by i.e. the 'wavesPrintRate'-variable in 'read.C') and a ROOT-file called 'out.root'. The ROOT-file contains the ROOT-Tree 'T' which has a branch for each variable and leaves for the events. These can be plotted/further analysed.
If it is not necessary to process all runs, it might be more convenient run `./run.sh` followed by the run-number to be processed (so i.e. `./run.sh 65` to only run runNr 65).
Note, that when using 'run.sh' both raw-files from both WaveCatchers will be processed (althogh the output is not combined and stored in two different folders in 'runs').
