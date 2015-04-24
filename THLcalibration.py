import sys

from ROOT import *
import ROOT
from array import array
import math
import os
import numpy as np

def derivative(x, y, source, assembly):

    x_deriv=[]
    deriv=[]
    for i in range(0, len(x)-1):
        x_deriv.append(x[i])
        if (assembly=="A06-W0110" or assembly=="C04-W0110"):
            deriv.append((y[i+1]-y[i])/(x[i+1]-x[i]))
        elif (assembly=="B06-W0125"):
            deriv.append(-(y[i+1]-y[i])/(x[i+1]-x[i]))

    g_deriv=TGraph(len(x_deriv), array('d', x_deriv), array('d', deriv))
    g_deriv.SetMarkerStyle(20)
    g_deriv.SetMarkerSize(0.8)

#    g_deriv.SetMarkerColor(color)
#    g_deriv.SetLineColor(color)
    g_deriv.SetTitle(assembly+": "+source+"; THL; #delta(Counts)/#delta(THL)");

    return g_deriv


def readTextFile(filename, source, assembly):
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
    g.SetTitle(assembly+": "+source+"; THL; Counts");

    g_deriv=derivative(THL, counts, source, assembly)

    return g, g_deriv

def errorPropagation(x, y, fit, a, b):
    
    n=len(x)
    temp=0
    for i in range(0, n):
        temp+=(y[i]-fit.Eval(x[i]))**2

    S=TMath.Sqrt(temp/(n-2))
    x_avg=sum(a)


def usage():
    print 'THL calibration'
    print 'Usage:\n python %s <assembly> <THL to evaluate>' % ( sys.argv[0] )

if __name__ == '__main__':

    if (len(sys.argv)<3 or len(sys.argv)>3):
        usage()
        sys.exit( 1 )

    assembly=sys.argv[1]
    THL_to_evaluate=int(sys.argv[2])
        
    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C");
    CLICdpStyle();
    gStyle.SetOptFit()
    gStyle.SetOptTitle(1);


    sources=["Cu", "Zr", "Pd", "In"]
    energies=[8.048, 15.77, 21.178, 24.209]
    energies_err=[0, 0, 0, 0]
    fitGausMin=[]
    fitGausMax=[]

    if (assembly=="A06-W0110"):
        fitGausMin=[250, 150, 80, 50]
        fitGausMax=[270, 190, 120, 80]
    elif (assembly=="C04-W0110"):
        fitGausMin=[340, 240, 180, 140]
        fitGausMax=[360, 270, 200, 170]
    elif (assembly=="B06-W0125"):
        fitGausMin=[475, 565, 625, 660]
        fitGausMax=[490, 580, 645, 690]

    meanTHL=[]
    meanTHLerr=[]

    for i in range(0, len(sources)):
        filename=assembly+"_"+sources[i]+".txt"
        # filename=assembly+"_"+sources[i]+"_test1.txt"

        # Create s-curve
        sCurve, sCurveDeriv=readTextFile(filename, sources[i], assembly)
        
        # Fit with gauss
        fitGaus=TF1("GaussFit", "gaus", fitGausMin[i], fitGausMax[i])
        fitGaus.SetLineColor(kBlue)
        sCurveDeriv.Fit("GaussFit", "R")
        
        # Draw the derivative
        canvDeriv=TCanvas(sources[i], sources[i])
        sCurveDeriv.Draw("ALP")        


        canvScurve=TCanvas("Scurve"+sources[i], "Scurve"+sources[i])
        sCurve.Draw("AP")
        # bla=raw_input()
        # Get Mean and its error
        meanTHL.append(fitGaus.GetParameter(1))
        meanTHLerr.append(fitGaus.GetParError(1))

        print "mean=", fitGaus.GetParameter(1)
        print "mean error=", fitGaus.GetParError(1)
    
        # bla=raw_input()


    gr=TGraphErrors(len(energies), array('d', energies), array('d', meanTHL), array('d', energies_err), array('d', meanTHLerr))
    

    #gr=TGraph(len(energies), array('d', meanTHL), array('d', energies))

    THL_fit=TF1("THL","[0]*x+[1]", 0, 300)

    THL_fit.SetParName(0,"a")
    THL_fit.SetParName(1,"b")
    gr.SetMarkerStyle(25)
    gr.SetMarkerSize(0.8)
    THL_fit.SetLineColor(kGreen+2)
    gr.SetLineColor(kBlack)
    gr.SetTitle(assembly)
    gr.GetXaxis().SetTitle("Energy [keV]")
    gr.GetYaxis().SetTitle("THL")

    fit_result=gr.Fit("THL", "S")
    gr.Draw("AP")
    gPad.Update()
    gStyle.SetOptFit() 

    a=THL_fit.GetParameter(0)
    b=THL_fit.GetParameter(1)
    print "a=", a
    print "b=", b
    
    covMat=fit_result.GetCovarianceMatrix()
    sigma_a=covMat[0][0]
    sigma_ab=covMat[0][1]
    sigma_b=covMat[1][1]

    y_eval=THL_to_evaluate
    x_eval=THL_fit.GetX(y_eval)
    
    err_x=TMath.Sqrt(((y_eval-b)**2)/(a**4)*sigma_a+sigma_b/(a**2)+2*sigma_ab*(y_eval-b)/(a**3))

    print "THL %d corresponds to: %f [keV] +- %f [keV]"%(y_eval, x_eval, err_x)

    bla=raw_input()
