# to define constants relevant to the calibration analysis

# assembly properties
npixX = 256
npixY = 256

# known inputs
known_assemblies = ["A06-W0110","B06-W0125","B07-W0125","C04-W0110","D09-W0126","L04-W0125"]
known_sources = ["Fe","Co","Cd","CuInXRF","Am","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]
known_peaks = ["Fe","Co1","CuXRF_CERN","Co2","Cd","InXRF","Am2","Am3","CoXRF","CrXRF","CuXRF_LNLS","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]

CERN_sources = ["Fe","Co","Cd","CuInXRF","Am"]
LNLS_sources = ["CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]

CERN_peaks = ["Fe","Co1","CuXRF_CERN","Co2","Cd","InXRF","Am2","Am3"]
LNLS_peaks = ["CoXRF","CrXRF","CuXRF_LNLS","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]

# peak energies in keV
FePeakE = 5.899
Co1PeakE = 6.4
Co2PeakE = 14.4
CdPeakE = 22.9
Am2PeakE = 26.2
Am3PeakE = 60.0

CoXRFPeakE = 6.93
CrXRFPeakE = 5.414
CuXRFPeakE = 8.1
FeXRFPeakE = 6.4
InXRFPeakE = 24.0
MnXRFPeakE = 5.89
NiXRFPeakE = 7.47
TiXRFPeakE = 4.51
VXRFPeakE = 4.95
