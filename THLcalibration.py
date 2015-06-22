

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
        if (assembly=="A06-W0110" or assembly=="C04-W0110" or assembly=="L04-W0125" or assembly=="B07-W0125"):
            deriv.append((y[i+1]-y[i])/(x[i+1]-x[i]))
        elif (assembly=="B06-W0125"):
            deriv.append(-(y[i+1]-y[i])/(x[i+1]-x[i]))

    g_deriv=TGraph(len(x_deriv), array('d', x_deriv), array('d', deriv))
    g_deriv.SetMarkerStyle(20)
    g_deriv.SetMarkerSize(0.8)

#    g_deriv.SetMarkerColor(color)
#    g_deriv.SetLineColor(color)
    #g_deriv.SetTitle(assembly+": "+source+"; THL; #delta(Counts)/#delta(THL)");
    g_deriv.SetTitle("; THL; #delta(Counts)/#delta(THL)");

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
    # g.SetTitle(assembly+": "+source+"; THL; Counts");
    g.SetTitle("; THL; Counts");

    g_deriv=derivative(THL, counts, source, assembly)

    return g, g_deriv

def errorPropagation(y_eval, a, b, sigma_a, sigma_b, sigma_ab):
    return TMath.Sqrt(((y_eval-b)**2)/(a**4)*sigma_a+sigma_b/(a**2)+2*sigma_ab*(y_eval-b)/(a**3))



def usage():
    print 'THL calibration'
    print 'Usage:\n python %s <assembly> <THL to evaluate>' % ( sys.argv[0] )

if __name__ == '__main__':

    if (len(sys.argv)<3 or len(sys.argv)>3):
        usage()
        sys.exit( 1 )

    assembly=sys.argv[1]
    THL_to_evaluate=int(sys.argv[2])
        
    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C")
    CLICdpStyle()


    savePath="/afs/cern.ch/work/n/nalipour/CLICdpNotes/Timepix_Calibration/figures/THLscan/"

    gStyle.SetOptFit(0)
    # gStyle.SetOptTitle(1)
#    gStyle.SetOptStat()


    sources=["Cu", "Zr", "Pd", "In"]
#    sources=["Cu"]
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
    elif (assembly=="L04-W0125"):
        fitGausMin=[345, 255, 190, 100]
        fitGausMax=[365, 275, 210, 200]
    elif (assembly=="B07-W0125"):
        fitGausMin=[430, 220, 160, 120]
        fitGausMax=[450, 240, 180, 150]

    meanTHL=[]
    meanTHLerr=[]

    for i in range(0, len(sources)):
        filename=assembly+"_"+sources[i]+".txt"
        print "filename=", filename

        # Create s-curve
        sCurve, sCurveDeriv=readTextFile(filename, sources[i], assembly)
        
        # Fit with gauss
        fitGaus=TF1("GaussFit", "gaus", fitGausMin[i], fitGausMax[i])
        fitGaus.SetLineColor(kBlue)
        sCurveDeriv.Fit("GaussFit", "RQ")
        

        # Draw the derivative
        canvDeriv=TCanvas(sources[i], sources[i])
        sCurveDeriv.Draw("ALP")   
        sCurveDeriv.GetYaxis().SetTitleOffset(1.27)
        sCurveDeriv.GetYaxis().SetLabelSize(0.05)
        sCurveDeriv.GetXaxis().SetRangeUser(fitGausMin[i]-30, fitGausMax[i]+30)
        sCurveDeriv.GetYaxis().SetRangeUser(0, fitGaus.GetMaximum()*1.3)

 
        pav = TPaveText(0.6,0.7,0.92,0.92,'NDC')
        pav.AddText("Mean %0.2f #pm %0.2f" %(fitGaus.GetParameters()[1], fitGaus.GetParErrors()[1])) 
        pav.AddText("Sigma %0.2f #pm %0.2f" %(fitGaus.GetParameters()[2], fitGaus.GetParErrors()[2])) 
        pav.Draw() 
        canvDeriv.Update()
        pav.SetBorderSize(0)
        pav.SetFillStyle(0)
        # stats = sCurveDeriv.GetListOfFunctions().FindObject("stats")
        
        # stats.SetX1NDC(0.18)
        # stats.SetX2NDC(0.52)
        # stats.SetY1NDC(0.7)
        # stats.SetY2NDC(0.9)
        # stats.SetFillStyle(0)
        # stats.SetBorderSize(0)
        # stats.SetTextSize(0.04)
        # # stats.Delete("chi2")
        # # sCurve.Update()
        # # sCurve.Modified()
        canvDeriv.Update()
        canvDeriv.Modified()

#        bla=raw_input() 


        canvScurve=TCanvas("Scurve"+sources[i], "Scurve"+sources[i])
        sCurve.Draw("AP")
        sCurve.GetYaxis().SetTitleOffset(1.32)
        sCurve.GetYaxis().SetLabelSize(0.05)
        #DrawMyText(assembly+": "+sources[i])
        sCurve.GetXaxis().SetRangeUser(fitGausMin[i]-30, fitGausMax[i]+30)
        sCurve.GetYaxis().SetRangeUser(0, TMath.Max(sCurve.Eval(fitGausMax[i]+30), sCurve.Eval(fitGausMin[i]-30))*1.3)
        canvScurve.Update()
        canvScurve.Modified()
        #bla=raw_input()
        # Get Mean and its error


        meanTHL.append(fitGaus.GetParameter(1))
        meanTHLerr.append(fitGaus.GetParError(1))

        print "mean=", fitGaus.GetParameter(1)
        print "mean error=", fitGaus.GetParError(1)
    
        # # -------- save figures -------- #
        # canvDeriv.Print("/afs/cern.ch/work/n/nalipour/CLICdpNotes/Timepix_Calibration/figures/THLscan/"+assembly+"/scurveDeriv_"+sources[i]+".pdf")
        # canvScurve.Print("/afs/cern.ch/work/n/nalipour/CLICdpNotes/Timepix_Calibration/figures/THLscan/"+assembly+"/scurve_"+sources[i]+".pdf")

        

 #       savePath="/afs/cern.ch/work/n/nalipour/CLICdpNotes/Timepix_Calibration/temp/"
        canvDeriv.Print(savePath+assembly+"_scurveDeriv_"+sources[i]+".pdf")
        canvScurve.Print(savePath+assembly+"_scurve_"+sources[i]+".pdf")


    canvTHLcalib=TCanvas("THL calibration", "THL calibration")
    gr=TGraphErrors(len(energies), array('d', energies), array('d', meanTHL), array('d', energies_err), array('d', meanTHLerr))

    THL_fit=TF1("THL","[0]*x+[1]", 0, 300)

    THL_fit.SetParName(0,"a")
    THL_fit.SetParName(1,"b")
    gr.SetMarkerStyle(25)
    gr.SetMarkerSize(0.8)
    THL_fit.SetLineColor(kBlue)
    gr.SetLineColor(kBlack)
    gr.SetTitle(assembly)
    gr.GetXaxis().SetTitle("Energy [keV]")
    gr.GetYaxis().SetTitle("THL DAC")

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
    
    #err_x=TMath.Sqrt(((y_eval-b)**2)/(a**4)*sigma_a+sigma_b/(a**2)+2*sigma_ab*(y_eval-b)/(a**3))
    err_x=errorPropagation(y_eval, a, b, sigma_a, sigma_b, sigma_ab)
    print "THL %d corresponds to: %f [keV] +- %f [keV]"%(y_eval, x_eval, err_x)

    delta_E=1/a
    delta_E_err=TMath.Sqrt(1./(a**2)*(sigma_a**2))
    print "delta_THL corresponds to %f [keV] +- %f [keV] or %f [e-] +- %f [e-]"%(delta_E, delta_E_err, delta_E*1000/3.6, delta_E_err*1000/3.6)





    stats = gr.GetListOfFunctions().FindObject("stats")

    # listLine=stats.GetListOfLines()
    # text = stats.GetLineWith("constant")
    # listLine.Remove(text)


    stats.SetX1NDC(0.35)
    stats.SetX2NDC(0.75)
    stats.SetY1NDC(0.75)
    stats.SetY2NDC(0.9)
    stats.SetFillStyle(0)
    stats.SetBorderSize(0)
    stats.SetTextSize(0.05)
    canvTHLcalib.Update()
    canvTHLcalib.Modified()
    #canvTHLcalib.Print("/afs/cern.ch/work/n/nalipour/CLICdpNotes/Timepix_Calibration/temp/"+assembly+"_THLcalibration.pdf")
    canvTHLcalib.Print(savePath+assembly+"_THLcalibration.pdf")
#"/afs/cern.ch/work/n/nalipour/CLICdpNotes/Timepix_Calibration/figures/THLscan/"+assembly+"/THLcalibration.pdf")

#    bla=raw_input()
