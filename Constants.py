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
Co1PeakE = 6.404
Co2PeakE = 14.41300
CdPeakE = 22.163
Am2PeakE = 26.3448
Am3PeakE = 59.5412

CoXRFPeakE = 6.9303780
CrXRFPeakE = 5.4148045
CuXRFPeakE = 8.0478227
FeXRFPeakE = 6.4040062
InXRFPeakE = 24.20975
MnXRFPeakE = 5.8988010
NiXRFPeakE = 7.4782521
TiXRFPeakE = 4.5108991
VXRFPeakE = 4.952216
