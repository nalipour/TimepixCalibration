import sys

from ROOT import *
import ROOT
from array import array
import math
import os
import numpy as np

def derivative(x, y):

    x_deriv=[]
    deriv=[]
    for i in range(0, len(x)-1):
        x_deriv.append(x[i])
        deriv.append((y[i+1]-y[i])/(x[i+1]-x[i]))


    g_deriv=TGraph(len(x_deriv), array('d', x_deriv), array('d', deriv))
    g_deriv.SetMarkerStyle(20)
    g_deriv.SetMarkerSize(0.8)

#    g_deriv.SetMarkerColor(color)
#    g_deriv.SetLineColor(color)
    g_deriv.SetTitle("; THL; #delta(Counts)/#delta(THL)");

    return g_deriv


def readTextFile(filename):
    f = open(filename)                                                                                                                                  
    lines = [line.strip() for line in open(filename)]                                                                                                   
    f.close()                                                                                                                                           

    THL=[]
    counts=[]

    for i in range(0, len(lines)):     
        THL.append(int(lines[i].split()[0]))
        counts.append(int(lines[i].split()[1]))

    g = TGraph(len(THL), array('d', THL), array('d', counts))
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.8)
    g.SetTitle("; THL; Counts");

    g_deriv=derivative(THL, counts)

    return g, g_deriv



if __name__ == '__main__':
    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C");
    CLICdpStyle();
    gStyle.SetOptFit()
    gStyle.SetOptTitle(1);

    assembly="A06-W0110"
    sources=["In"]#, "Zr", "Cu"]
    energies=[24.209]#, 8.048, 15.77]

    meanTHL=[]
    meanTHLerr=[]

    for i in range(0, len(sources)):
        filename=assembly+"_"+sources[i]+"_test1.txt"

        # Create s-curve
        sCurve, sCurveDeriv=readTextFile(filename)
        
        # Fit with gauss
        fitGaus=TF1("GaussFit", "gaus", 50, 80)
        fitGaus.SetLineColor(kBlue)
        sCurveDeriv.Fit("GaussFit", "R")
        
        # Draw the derivative
        canvDeriv=TCanvas(sources[i], sources[i])
        sCurveDeriv.Draw("ALP")
        
        canvScurve=TCanvas("Scurve"+sources[i], "Scurve"+sources[i])
        sCurve.Draw("AP")

        # Get Mean and its error
        meanTHL.append(fitGaus.GetParameter(0))
        meanTHLerr.append(fitGaus.GetParError(0))

        print "mean=", fitGaus.GetParameter(0)
        print "mean error=", fitGaus.GetParError(0)
    

    bla=raw_input()

"""    
    gr=TGraphErrors(len(energies), array('d', meanTHL), array('d', energies), array('d', meanTHLerr), 0)
    THL_fit=TF1("THL","[0]*x+[1]")
    THL_fit.SetParName(0,"a")
    THL_fit.SetParName(1,"b")
    gr.SetMarkerStyle(8)
    gr.SetLineColor(kBlack)
    gr.SetTitle(assembly)
    gr.GetXaxis().SetTitle("THL")
    gr.GetYaxis().SetTitle("Energy [keV]")

    gr.Fit("THL")
    gr.Draw("AP")
    gPad.Update()
    gStyle.SetOptFit() 
    
    print "THL 326 corresponds to:", THL_fit.Eval(326)
"""

