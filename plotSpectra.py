# file to plot superimposed spectra of the calibration data

import ROOT as R
from optparse import OptionParser
from os import environ
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


def plotSpectra(assembly,topTOT):
    c = R.TCanvas("c","",0,0,1000,400)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.1)

    # get data
    home = environ['HOME']
    base = "%s/eos/clicdp/data/VertexCalibration" %home
    assembly_start = assembly.split("-")[0]

    rootfileFe = R.TFile("%s/%s/%s_SinglePixelCalibration/Fe55_%s_spc.root"%(base,assembly,assembly_start,assembly))
    rootfileAm = R.TFile("%s/%s/%s_SinglePixelCalibration/Am241_%s_spc.root"%(base,assembly,assembly_start,assembly))
    rootfileCd = R.TFile("%s/%s/%s_SinglePixelCalibration/Cd109_%s_spc.root"%(base,assembly,assembly_start,assembly))
    if assembly == "B06-W0125":
        rootfileCuIn = R.TFile("%s/%s/%s_SinglePixelCalibration/Cu_In_%s_spc.root"%(base,assembly,assembly_start,assembly))
    else:
        rootfileCuIn = R.TFile("%s/%s/CuIn_%s.root" %(base,assembly,assembly))
    if assembly == "B06-W0125":
        rootfileCo = R.TFile("%s/%s/%s_SinglePixelCalibration/Co57_%s_spc.root"%(base,assembly,assembly_start,assembly))
    else:
        rootfileCo = R.TFile("%s/%s/Co57_%s.root" %(base,assembly,assembly))

    # trees
    treeFe = rootfileFe.Get("pixels")
    treeAm = rootfileAm.Get("pixels")
    treeCd = rootfileCd.Get("pixels")
    treeCuIn = rootfileCuIn.Get("pixels")
    treeCo = rootfileCo.Get("pixels")

    # hists
    histFe_allp = R.TH1F("histFe_allp","",100,0,topTOT)
    histAm_allp = R.TH1F("histAm_allp","",100,0,topTOT)
    histCd_allp = R.TH1F("histCd_allp","",100,0,topTOT)
    histCuIn_allp = R.TH1F("histCuIn_allp","",100,0,topTOT)
    histCo_allp = R.TH1F("histCo_allp","",100,0,topTOT)

    histFe_onep = R.TH1F("histFe_onep","",20,0,topTOT)
    histAm_onep = R.TH1F("histAm_onep","",20,0,topTOT)
    histCd_onep = R.TH1F("histCd_onep","",20,0,topTOT)
    histCuIn_onep = R.TH1F("histCuIn_onep","",20,0,topTOT)
    histCo_onep = R.TH1F("histCo_onep","",20,0,topTOT)

    treeFe.Draw("tot>>histFe_allp","","goff") 
    treeAm.Draw("tot>>histAm_allp","","goff")
    treeCd.Draw("tot>>histCd_allp","","goff")
    treeCuIn.Draw("tot>>histCuIn_allp","","goff")
    treeCo.Draw("tot>>histCo_allp","","goff")

    treeFe.Draw("tot>>histFe_onep","col==128&&row==128","goff") 
    treeAm.Draw("tot>>histAm_onep","col==128&&row==128","goff")
    treeCd.Draw("tot>>histCd_onep","col==128&&row==128","goff")
    treeCuIn.Draw("tot>>histCuIn_onep","col==128&&row==128","goff")
    treeCo.Draw("tot>>histCo_onep","col==128&&row==128","goff")

    for hist in [histFe_allp,histAm_allp,histCd_allp,histCuIn_allp,histCo_allp,
                 histFe_onep,histAm_onep,histCd_onep,histCuIn_onep,histCo_onep]:
            hist.Scale(1./hist.GetMaximum())

    for hist in [histFe_allp,histFe_onep]:
        hist.GetXaxis().SetTitle("TOT (ADC counts)")
        hist.GetYaxis().SetTitle("A.U.")
        hist.GetYaxis().SetTitleOffset(0.5)
    
    histFe_allp.SetLineColor(R.kBlue)
    histAm_allp.SetLineColor(R.kGreen+2)
    histCd_allp.SetLineColor(R.kRed+2)
    histCuIn_allp.SetLineColor(R.kOrange)
    histCo_allp.SetLineColor(R.kCyan+1)

    histFe_onep.SetLineColor(R.kBlue)
    histAm_onep.SetLineColor(R.kGreen+2)
    histCd_onep.SetLineColor(R.kRed+2)
    histCuIn_onep.SetLineColor(R.kOrange)
    histCo_onep.SetLineColor(R.kCyan+1)

    leg_allp = R.TLegend(0.7,0.5,0.94,0.89)
    leg_allp.AddEntry(histFe_allp,"Fe: %i evt" %histFe_allp.GetEntries(),"l")
    leg_allp.AddEntry(histCo_allp,"Co: %i evt" %histCo_allp.GetEntries(),"l")
    leg_allp.AddEntry(histCuIn_allp,"CuIn: %i evt" %histCuIn_allp.GetEntries(),"l")
    leg_allp.AddEntry(histCd_allp,"Cd: %i evt" %histCd_allp.GetEntries(),"l")
    leg_allp.AddEntry(histAm_allp,"Am: %i evt" %histAm_allp.GetEntries(),"l")

    leg_onep = R.TLegend(0.7,0.5,0.94,0.89)
    leg_onep.AddEntry(histFe_onep,"Fe: %i evt" %histFe_onep.GetEntries(),"l")
    leg_onep.AddEntry(histCo_onep,"Co: %i evt" %histCo_onep.GetEntries(),"l")
    leg_onep.AddEntry(histCuIn_onep,"CuIn: %i evt" %histCuIn_onep.GetEntries(),"l")
    leg_onep.AddEntry(histCd_onep,"Cd: %i evt" %histCd_onep.GetEntries(),"l")
    leg_onep.AddEntry(histAm_onep,"Am: %i evt" %histAm_onep.GetEntries(),"l")
    
    histFe_allp.Draw()
    histAm_allp.Draw("same")
    histCd_allp.Draw("same")
    histCuIn_allp.Draw("same")
    histCo_allp.Draw("same")
    leg_allp.Draw()
    c.Update()
    c.SaveAs("plots/Spectra/%s_Spectra_AllPixels.pdf" %assembly)

    histFe_onep.Draw()
    histAm_onep.Draw("same")
    histCd_onep.Draw("same")
    histCuIn_onep.Draw("same")
    histCo_onep.Draw("same")
    leg_onep.Draw()
    c.Update()
    c.SaveAs("plots/Spectra/%s_Spectra_OnePixel.pdf" %assembly)


def getTopTOT(assembly):
    
    if assembly == "A06-W0110": topTOT = 1300
    if assembly == "B06-W0125": topTOT = 2700
    if assembly == "B07-W0125": topTOT = 1400
    if assembly == "C04-W0110": topTOT = 1300
    if assembly == "D09-W0126": topTOT = 1700
    if assembly == "L04-W0125": topTOT = 1500

    return topTOT


topTOT = getTopTOT(assembly)
plotSpectra(assembly,topTOT)
