# will read the points and uncertainties from text files,
# fit the surrogate function and save the result
# for options run python surrogateFitterGlobal.py -h

import ROOT as R
from array import array
import ast
from optparse import OptionParser
import Constants as C

parser = OptionParser()
parser.add_option("-b", "--assembly",
                  help="Assembly name", dest="ASSEMBLY")

(options, args) = parser.parse_args()

if(options.ASSEMBLY):
    assembly=options.ASSEMBLY
else :
    print "Please specify assembly"
    print "choose from", C.known_assemblies
    parser.print_help()
    exit()

if assembly not in C.known_assemblies:
    print "Assembly not recognised"
    print "choose from", C.known_assemblies
    exit()

# import CLICdp ROOT style
R.gROOT.ProcessLine(".L CLICdpStyle.C")
R.CLICdpStyle()

# extra style
R.gStyle.SetOptFit(1)
R.gROOT.SetBatch(True)


def fitGlobalSurrogate(assembly):
    # threshold
    fTh = open('results/testpulse/global_threshold_%s.txt' %assembly,'r')
    # 8 peaks
    fFe = open('results/kde/%s_Fe_GlobalResults.txt' %assembly,'r')
    fAm2 = open('results/kde/%s_Am2_GlobalResults.txt' %assembly,'r')
    fAm3 = open('results/kde/%s_Am3_GlobalResults.txt' %assembly,'r')
    fCu = open('results/kde/%s_Cu_GlobalResults.txt' %assembly,'r')
    fIn = open('results/kde/%s_In_GlobalResults.txt' %assembly,'r')
    fCo1 = open('results/kde/%s_Co1_GlobalResults.txt' %assembly,'r')
    fCo2 = open('results/kde/%s_Co2_GlobalResults.txt' %assembly,'r')
    fCd = open('results/kde/%s_Cd_GlobalResults.txt' %assembly,'r')


    Thline = fTh.readline().split()
    Feline = fFe.readline().split()
    Am2line = fAm2.readline().split()
    Am3line = fAm3.readline().split()
    Culine = fCu.readline().split()
    Inline = fIn.readline().split()
    Co1line = fCo1.readline().split()
    Co2line = fCo2.readline().split()
    Cdline = fCd.readline().split()

    for line in [Thline,Feline,Am2line,Am3line,Culine,Inline,Co1line,Co2line,Cdline]:
        for i in xrange(len(line)):
            line[i] = ast.literal_eval(line[i])

    tots = array('f',[])
    tot_lerrs = array('f',[])
    tot_uerrs = array('f',[])

    energies = array('f',[])
    energy_lerrs = array('f',[])
    energy_uerrs = array('f',[])

    for line, energy in zip([Feline,Am2line,Am3line,Culine,Inline,Co1line,Co2line,Cdline],
                            [C.FePeakE,C.Am2PeakE,C.Am3PeakE,C.CuPeakE,C.InPeakE,C.Co1PeakE,C.Co2PeakE,C.CdPeakE]):
        tots.append(line[0])
        tot_lerrs.append(line[1])
        tot_uerrs.append(line[2])

        energies.append(energy)
        energy_lerrs.append(0.)
        energy_uerrs.append(0.)

    for line in [Thline]:

        tots.append(0.)
        tot_lerrs.append(0.)
        tot_uerrs.append(0.)

        energies.append(line[0])
        energy_lerrs.append(line[1])
        energy_uerrs.append(line[2])


    canv = R.TCanvas()
    canv.SetRightMargin(0.01)
    canv.SetLeftMargin(0.19)
    gr = R.TGraphAsymmErrors(len(energies),energies,tots,energy_lerrs,energy_uerrs,tot_lerrs,tot_uerrs)
    gr.SetMarkerStyle(20)
    gr.GetXaxis().SetTitle("Energy (keV)")
    gr.GetYaxis().SetTitle("TOT (ADC)")
    gr.GetYaxis().SetTitleOffset(1.4)
    gr.Draw('AP')

    surrogate = R.TF1("surrogate","[0]*x+[1]-([2]/(x-[3]))",0,60)
    surrogate.SetParName(0,'a')
    surrogate.SetParName(1,'b')
    surrogate.SetParName(2,'c')
    surrogate.SetParName(3,'t')
    surrogate.SetParameters(13.9,333.2,1220,-0.1)
    gr.Fit("surrogate","RQ")
    canv.Update()

    stats = gr.GetListOfFunctions().FindObject("stats")
    stats.SetX1NDC(0.57)
    stats.SetX2NDC(0.96)
    stats.SetY1NDC(0.21)
    stats.SetY2NDC(0.46)
    gr.SetMinimum(0.)
    gr.SetMaximum(surrogate.Eval(65.0))
    canv.Update()
    canv.SaveAs("plots/KDESurrogateFits/%s_GlobalSurrogateFit.pdf" %assembly)

    fTh.close()
    fFe.close()
    fAm2.close()
    fAm3.close()
    fCu.close()
    fIn.close()
    fCo1.close()
    fCo2.close()
    fCd.close()



fitGlobalSurrogate(assembly)

