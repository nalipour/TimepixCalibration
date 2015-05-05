# to define constants relevant to the calibration analysis

# assembly properties
npixX = 256
npixY = 256

# known inputs
known_assemblies = ["A06-W0110","B06-W0125","B07-W0125","C04-W0110","D09-W0126","L04-W0125"]
known_sources = ["Fe","Co","Cd","CuInXRF","Am","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]
known_peaks = ["Fe","Co1","Cu","Co2","Cd","In","Am2","Am3"]

CERN_sources = ["Fe","Co","Cd","CuInXRF","Am"]
LNLS_sources = ["CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]

# peak energies in keV
FePeakE = 5.899
Co1PeakE = 6.4
CuPeakE = 8.1
Co2PeakE = 14.4
CdPeakE = 22.9
InPeakE = 24.0
Am2PeakE = 26.2
Am3PeakE = 60.0
