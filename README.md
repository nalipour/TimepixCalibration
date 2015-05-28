# TimepixCalibration
Use this code to run the Timepix calibration analysis
##### Setting up the environment
Do:
```
source /afs/cern.ch/eng/clic/software/setupSLC6.sh
```
##### Getting access to the data
Do:
```
export EOS_MGM_URL="root://eospublic.cern.ch"
source /afs/cern.ch/project/eos/installation/client/etc/setup.sh
eosmount $HOME/eos
```
##### Learn about the data
Do:
```
python plotStats.py -b A06-W0110 -s Fe
python plotSpectra.py -b A06-W0110
```
for all assemblies (b) and sources (s) you are interested in
##### Running the global analysis
Do:
```
python calibrateKDEGlobal.py -b A06-W0110 -s Fe
python fitGlobalSurrogate.py -b A06-W0110
```
for all assemblies (b) and sources (s) you are interested in
##### Running the pixel analysis
Do:
```
python calibrateKDEPixels.py -b A06-W0110 -s Fe
python fitPixelSurrogates.py -b A06-W0110
```
for all assemblies (b) and sources (s) you are interested in