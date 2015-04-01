# a file to plot the available statistics per pixel, per source, per assembly

import ROOT as R
import time
from os import environ
from optparse import OptionParser
import Constants as C
import getpass

parser = OptionParser()
parser.add_option("-b", "--assembly",
                  help="Assembly name", dest="ASSEMBLY")

parser.add_option("-s", "--source",
                  help="Source name", dest="SOURCE")

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

if(options.SOURCE):
    source=options.SOURCE
else :
    print "Please specify source"
    print "choose from", C.known_sources
    parser.print_help()
    exit()

if source not in C.known_sources:
    print "Source not recognised"
    print "choose from", C.known_sources
    exit()

# import CLICdp ROOT style
R.gROOT.ProcessLine(".L CLICdpStyle.C")
R.CLICdpStyle()

# extra style
R.gStyle.SetOptStat(1100)

def plotStats(assembly,source):

    # get data
    home = environ['HOME']
    base = "%s/eos/clicdp/data/VertexCalibration" %home
    assembly_start = assembly.split("-")[0]

    if source == "Fe":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Fe55_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "Am":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Am241_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "Cd":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cd109_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "CuIn":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cu_In_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/CuIn_%s.root" %(base,assembly,assembly))
    elif source == "Co":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Co57_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/Co57_%s.root" %(base,assembly,assembly))

    # tree
    tree = rootfile.Get("pixels")

    statmap = R.TH2F("statmap","",C.npixX,0,C.npixX,C.npixY,0,C.npixY)
    stathist = R.TH1F("stathist","",100,0,2500)

    username = getpass.getuser()
    junkfile = R.TFile("/tmp/%s/junkfile_%s_%s.root" %(username,assembly,source),"RECREATE")

    for i in xrange(0,C.npixX):
        start = time.time()
        print "doing row", i
        cut = "row==%i" %i
        t2 = tree.CopyTree(cut)
        for j in xrange(0,C.npixY):
            ent = t2.GetEntries("col==%i" %j)
            statmap.Fill(j,i,ent)
            stathist.Fill(ent)
        end = time.time()
        print "time", end - start

    c = R.TCanvas()
    c.SetRightMargin(0.2)
    statmap.SetStats(0)
    statmap.GetXaxis().SetTitle("Pixel X")
    statmap.GetYaxis().SetTitle("Pixel Y")
    statmap.Draw("colz")
    c.SaveAs("plots/Statistics/%s_%s_StatMap.pdf" %(assembly,source))

    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.1)
    stathist.GetXaxis().SetTitle("Hits per pixel")
    stathist.Draw()
    c.Update()
    stats = stathist.GetListOfFunctions().FindObject("stats")
    stats.SetX1NDC(0.6)
    stats.SetX2NDC(0.86)
    stats.SetY1NDC(0.75)
    stats.SetY2NDC(0.88)
    c.SaveAs("plots/Statistics/%s_%s_StatHist.pdf" %(assembly,source)) 


plotStats(assembly,source)

