'''
Created on March 18, 2015

@author: fpfleger
'''
import time,os
from math  import sqrt
from numpy import matrix
from array import array
import numpy
from ROOT import *
from scipy.cluster.hierarchy import fclusterdata
import pyximport; pyximport.install(pyimport=True)
import collections
from itertools import product
from optparse import OptionParser
import sys
from os import environ
from math import ceil
import Constants as C

gROOT.ProcessLine(".L ../../rootstyle/CLICdpStyle/rootstyle/CLICdpStyle.C")
CLICdpStyle()

parser = OptionParser()

parser.add_option('-b', '--assembly', help='Assembly name', dest='ASSEMBLY')

parser.add_option('-s', '--source', help='Source name', dest='SOURCE')

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

class switch(object):
    value = None
    def __new__(class_, value):
        class_.value = value
        return True

def case(*args):
    return any((arg == switch.value for arg in args))

def ReadROOT(assembly,source):

    #home = environ['HOME']
    #base = "%s/eos/clicdp/data/VertexCalibration" %home
    #assembly_start = assembly.split("-")[0]

    if source == "Fe":
        rootfile = TFile("Fe55_%s_processed.root"%(base,assembly,assembly_start,assembly))
    elif source == "Am":
        rootfile = TFile("Am241_%s_processed.root"%(base,assembly,assembly_start,assembly))
    elif source == "Cd":
        rootfile = TFile("Cd109_%s_processed.root"%(base,assembly,assembly_start,assembly))
    elif source == "CuIn":
        if assembly == "B06-W0125":
            rootfile = TFile("Cu_In_%s_processed.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = TFile("CuIn_%s_processed.root" %(assembly))
    elif source == "Co":
        if assembly == "B06-W0125":
            rootfile = TFile("Co57_%s_processed.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = TFile("Co57_%s_processed.root" %(base,assembly,assembly)) 

    return rootfile

rootfile = ReadROOT(assembly,source)


def TOTonSize(rootfile):

    tree = rootfile.Get("allHits")
    clusterData = []

    hist1 = TH1F("hist1", "", 100, 0, 2000) 
    hist2 = TH1F("hist2", "", 100, 0, 2000) 
    hist3 = TH1F("hist3", "", 100, 0, 2500) 
    hist4 = TH1F("hist4", "", 100, 0, 2500)
    hist5 = TH1F("hist5", "", 100, 0, 2500)
    
   
    tree.Draw("tot >> hist1", "clustersize == 1")
    tree.Draw("tot >> hist2", "clustersize == 2")
    tree.Draw("tot >> hist3", "clustersize == 3")
    tree.Draw("tot >> hist4", "clustersize == 4")
    tree.Draw("tot >> hist5")


    hist1.Scale(1./hist1.GetEntries())
    hist2.Scale(1./hist2.GetEntries())
    hist3.Scale(1./hist3.GetEntries())
    hist4.Scale(1./hist4.GetEntries())


    c = TCanvas()
    gPad.SetGrid()
    hist1.SetLineColor(kRed)
    hist1.SetLineWidth(2) 
    hist2.SetLineColor(kGreen+1)
    hist2.SetLineWidth(2) 
    hist3.SetLineColor(kBlue)
    hist3.SetLineWidth(2) 
    hist4.SetLineColor(kMagenta)
    hist4.SetLineWidth(2) 

    
    cluster_1hit = hist1.GetEntries()
    cluster_2hit = hist2.GetEntries()
    cluster_3hit = hist3.GetEntries()
    cluster_4hit = hist4.GetEntries()
    cluster_all = hist5.GetEntries()
    
    clusterData.append(cluster_all)
    clusterData.append(cluster_1hit)
    clusterData.append(cluster_2hit)
    clusterData.append(cluster_3hit)
    clusterData.append(cluster_4hit)
    

    percent_1hit = ceil(float(cluster_1hit)/cluster_all*100*100)/100
    percent_2hit = ceil(float(cluster_2hit)/cluster_all*100*100)/100
    percent_3hit = ceil(float(cluster_3hit)/cluster_all*100*100)/100
    percent_4hit = ceil(float(cluster_4hit)/cluster_all*100*100)/100
    
    clusterData.append(percent_1hit)
    clusterData.append(percent_2hit)
    clusterData.append(percent_3hit)
    clusterData.append(percent_4hit)
    clusterData.append(cluster_all)
    

    hist1.GetXaxis().SetTitle("cluster TOT (ADC count)")
    hist1.GetYaxis().SetTitle("A.U.")
    hist1.GetYaxis().SetTitleOffset(1.2)

    mylegend = TLegend(0.5,0.65,0.7,0.75)
    mylegend.SetFillColor(0)
    mylegend.SetBorderSize(0)
    mylegend.SetTextSize(0.03)
    mylegend.AddEntry(hist1,"single hit clusters: %.2f %%"%percent_1hit ,"l")
    mylegend.AddEntry(hist2,"two hit clusters: %.2f %%"%percent_2hit ,"l")


    hist1.Draw()
    hist2.Draw("same")
       

    if (percent_3hit > 0.5):
        hist3.Draw("same")
        mylegend.AddEntry(hist3,"three hit clusters: %.2f %%"%percent_3hit ,"l")
        
        
    if (percent_4hit > 0.5):
        hist4.Draw("same")
        mylegend.AddEntry(hist4,"three hit clusters: %.2f %%"%percent_4hit ,"l")
    
       
    mylegend.Draw()
     
    while switch(source): 
        if case('Fe'):
            c.SaveAs("Fe55-%s_stat.pdf"%(assembly))
            break
        if case('CuIn'):
            c.SaveAs("CuIn-%s_stat.pdf"%assembly)
            break
        if case('Am'):
            c.SaveAs("Am241-%s_stat.pdf"%(assembly))
            break
        if case('Co'):
            c.SaveAs("Co57-%s_stat.pdf"%(assembly))
            break
        if case('Cd'):
            c.SaveAs("Cd109-%s_stat.pdf"%(assembly))
            break
        printf("unknown source\n")
        break
    
    return clusterData

clusterData = TOTonSize(rootfile)

print "\n\t found a total of %.0f clusters\n"%clusterData[0]

print "\t one hit clusters:\t %.0f, %.2f %%"%(clusterData[1],clusterData[5])
print "\t two hit clusters:\t %.0f, %.2f %%"%(clusterData[2],clusterData[6])
print "\t three hit clusters:\t %.0f, %.2f %%"%(clusterData[3],clusterData[7])
print "\t four hit clusters:\t %.0f, %.2f %%\n"%(clusterData[4],clusterData[8])   



def ClusterSizes(rootfile,source,assembly, clusterData):
        
    CSizes= rootfile.Get("allHits")

    
    hist0 = TH1F("hist1", "distribution", 100, 0, 5)
    CSizes.Draw("clustersize >> hist1")

    hist0.GetXaxis().SetTitle("clustersizes")
    hist0.GetYaxis().SetTitle("Events")

    hist0.SetLineColor(kRed)
    hist0.SetLineWidth(2) 


    mylegend = TLegend(0.6,0.65,0.8,0.75)
    mylegend.SetFillColor(0)
    mylegend.SetBorderSize(0)
    mylegend.SetTextSize(0.03)

    mylegend.AddEntry(hist0,"1hc: %.2f %%"%clusterData[5])
    mylegend.AddEntry(hist0,"2hc: %.2f %%"%clusterData[6])
    mylegend.AddEntry(hist0,"3hc: %.2f %%"%clusterData[7])
    mylegend.AddEntry(hist0,"4hc: %.2f %%"%clusterData[8])


    c = TCanvas()
    gPad.SetGrid()
    hist0.Draw()
    mylegend.Draw()

    while switch(source): 
        if case('Fe'):
            c.SaveAs("Fe55-%s_CSizes.pdf"%(assembly))
            break
        if case('CuIn'):
            c.SaveAs("CuIn-%s_CSizes.pdf"%assembly)
            break
        if case('Am'):
            c.SaveAs("Am241-%s_CSizes.pdf"%(assembly))
            break
        if case('Co'):
            c.SaveAs("Co57-%s_CSizes.pdf"%(assembly))
            break
        if case('Cd'):
            c.SaveAs("Cd109-%s_CSizes.pdf"%(assembly))
            break
        printf("unknown source\n")
        break
   

ClusterSizes(rootfile,source,assembly,clusterData)

def SingelsVsEveryPix(rootfile,source, assembly):

    
    singles = rootfile.Get("singleHits")
    pixels = rootfile.Get("everyPix")

    hist1 = TH1F("hist1", "Single hit clusters", 100, 0, 2000) 
    hist2 = TH1F("hist2", "every single pixel", 100, 0, 2000) 


    singles.Draw("tot >> hist1")
    pixels.Draw("tot >> hist2")


    c = TCanvas()
    gPad.SetGrid()
    hist1.SetLineColor(kRed)
    hist1.SetLineWidth(2) 
    hist2.SetLineColor(kCyan+1)
    hist2.SetLineWidth(2) 

    print "\n\t # total pixels: %.0f"%hist2.GetEntries()
    print "\t # single - hit clusters: %.0f\n"% hist1.GetEntries()
    

    hist2.GetXaxis().SetTitle("pixel TOT")
    hist1.GetYaxis().SetTitle("A.U.")
    hist1.GetYaxis().SetTitleOffset(1.2)


    hist2.Draw()
    hist1.Draw("same")

    mylegend = TLegend(0.55,0.65,0.8,0.75)
    mylegend.SetFillColor(0)
    mylegend.SetBorderSize(0)
    mylegend.SetTextSize(0.03)

    mylegend.AddEntry(hist2,"every pixel", "l")
    mylegend.AddEntry(hist1,"single hit clusters","l")
    mylegend.Draw()

    while switch(source): 
        if case('Fe'):
            c.SaveAs("Fe55-%s_comp.pdf"%(assembly))
            break
        if case('CuIn'):
            c.SaveAs("CuIn-%s_comp.pdf"%assembly)
            break
        if case('Am'):
            c.SaveAs("Am241-%s_comp.pdf"%(assembly))
            break
        if case('Co'):
            c.SaveAs("Co57-%s_comp.pdf"%(assembly))
            break
        if case('Cd'):
            c.SaveAs("Cd109-%s_comp.pdf"%(assembly))
            break
        printf("unknown source\n")
        break


SingelsVsEveryPix(rootfile,source,assembly)     



    
