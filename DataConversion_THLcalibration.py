import sys

from ROOT import *
import ROOT
from array import array
import math
import os
import numpy as np
from scipy.cluster.hierarchy import fclusterdata



def SimpleClustering(matrix):

    frame=matrix

    for i,j in [[i,j] for i,j in product(xrange(256),xrange(256)) if matrix[i][j]!=0] :
        for u,v in [[u,v] for u,v in product([-1,0,1],[-1,0,1]) if (((i+u>=0 and i+u<=255) and (j+v>=0 and j+v<=255)) and (u!=0 or v!=0)) ] :
            if(matrix[i+u][j+v]!=0) :
                frame[i][j]=0
                frame[i+u][j+v]=0

    return frame

# def readTextFile(filename, hist):
#     print filename
#     f = open(filename)
#     lines = [line.strip() for line in open(filename)]
    
#     f.close()

    
#     for i in range(0, len(lines)):
#         for k in range(0, len(lines[i].split())):
#             hist.Fill(float(lines[i].split()[k]))

# # def readInMatrix(filename):
# #     f=open(filename , 'r')
# #     matrix = [ map(float,line.split(' ')) for line in f ]
# #     f.close()
# #     return matrix


"""def derivative(x, y, color):

    x_deriv=[]
    deriv=[]
    for i in range(0, len(x)-1):
        x_deriv.append(x[i])
        deriv.append((y[i+1]-y[i])/(x[i+1]-x[i]))

    
    g_deriv=TGraph(len(x_deriv), array('d', x_deriv), array('d', deriv))
    g_deriv.SetMarkerStyle(20)
    g_deriv.SetMarkerSize(0.8)

    g_deriv.SetMarkerColor(color)
    g_deriv.SetLineColor(color)

    return g_deriv
"""
def writeFile(THL, Counts, OutputFile):
    txtfile = open(OutputFile, "w")
    
    for i in range (0, len(THL)):
        txtfile.write(str(THL[i])+" "+str(Counts[i]))
        
    txtfile.close()
    

def usage():
    print 'THL calibration'
    print 'Usage:\n python %s <assembly> <frames> <THLmin> <THLmax> <source>' % ( sys.argv[0] )
# def erf(x, par):
#     return par[0]*TMath.Erf(par[1]*x[0]+par[2])+par[3]

# def fitHist(hist):

#     maxBin=hist.GetMaximumBin()
#     fit=TF1("GaussFit", "gaus", maxBin-1000, maxBin+1000)
#     hist.Fit("GaussFit", "R")
#     print "mean=", fit.GetParameter(1)
#     return fit.GetParameter(1)


if __name__ == '__main__':


    if (len(sys.argv)<6 or len(sys.argv)>6):
        usage()
        sys.exit( 1 )


    assembly=sys.argv[1]
    frames=int(sys.argv[2])
    THLmin=int(sys.argv[3])
    THLmax=int(sys.argv[4])
    source=sys.argv[5]
    

    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C");
    CLICdpStyle();
    gStyle.SetOptFit() 
    gStyle.SetOptTitle(1);


    mg=TMultiGraph()
    leg=TLegend(0.65, 0.7, 0.8, 0.85)
    leg.SetFillStyle(0)


    mg_deriv=TMultiGraph()
    leg_deriv=TLegend(0.65, 0.7, 0.8, 0.85)
    leg_deriv.SetFillStyle(0)

    colors=[kGreen+2, kRed, kBlue, kBlack, kViolet, kMagenta]


    directory="ThresholdScan/A06-W0110/"
#    source="In"
    path=directory+source
#    assembly="A06-W0110"
#    frames=5
    THLstep=2

#    THLmin=20 #In
#    THLmax=120 #In

    # THLmin=20 #Zr 
    # THLmax=328 #Zr




    THLs=[]
    counts=[]

    txtfile = open( assembly+"_"+source+"_lxbatch.txt", "w")

    for thl in range(THLmin, THLmax, THLstep):
        
        THLs.append(thl)
        count=0
        for frame in range(0, frames):
            filename="%s_THL%03i_%d"%(source, thl, frame)
            f = open(path+"/"+filename,'r')
            matrix = [ map(int,line.split(' ')) for line in f ]

            print "Occupancy: %d pixels"%(np.count_nonzero(matrix))

            # ====== fclusterdata used for clustering ====== #
            coords=np.array(np.transpose(np.nonzero(matrix))) # coordinate of non-zero elements
            row=[]
            col=[]
            countForThisFrame=[]

            for i in range(0, len(coords)):
                x=coords[i][0]
                y=coords[i][1]
                row.append(x)
                col.append(y)
                countForThisFrame.append(matrix[x][y])



            ClusteringResult=fclusterdata(coords, sqrt(2.), criterion="distance") #clustering
            NbOfOccurance=np.bincount(ClusteringResult) # occurance of each value
            for i in range(1, len(NbOfOccurance)):
                if (NbOfOccurance[i]!=1):
                    multiHitClusters=np.where(ClusteringResult==i)
                    for j in range(0, len(multiHitClusters)):
                        row.pop(j)
                        col.pop(j)
                        countForThisFrame.pop(j)

            count+=sum(countForThisFrame)
            # ============ End fclusterdata ============== #

        counts.append(count)
        txtfile.write(str(thl)+" "+str(count)+"\n")
        print "Processing THL="+str(thl)+", Counts="+str(count)

    txtfile.close()
#    writeFile(THLs, counts, assembly+"_"+source+"test.txt")

"""
    g = TGraph(len(THLs), array('d', THLs), array('d', counts))
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.8)
    g.SetTitle("; THL; Counts");
    g.Draw("AP")
    # bla=raw_input()

    print "Calculating the derivative ..."
    deriv=derivative(THLs, counts, kRed)
    deriv.SetTitle("; THL; #delta(Counts)/#delta(THL)");
    print "Fitting the derivative with a Gaussian ... "
    # fitGaus=TF1("GaussFit", "gaus", 250, 270) #Cu
    # fitGaus=TF1("GaussFit", "gaus", 50, 80) #In
    fitGaus=TF1("GaussFit", "gaus", 150, 190) #Zr
    fitGaus.SetLineColor(kBlue)
    deriv.Fit("GaussFit", "R")



    canvDeriv=TCanvas("deriv", "deriv")
    deriv.Draw("ALP")

    gPad.SetGrid(1,1)

    #bla=raw_input()
"""







#######################################################
# for (path, dirnames, filenames) in os.walk(assembly_path):
#     dataFiles=[x for x in filenames if not ".dsc" in x]

#         count=0
#         for filename in dataFiles:
#             # THLs.append(int(filename.split("In_THL")[1]))
#             # THLs.append(int(filename.split("A06-W0110_Fe55_THL")[1]))
#             temp=(filename.split(source+"_THL")[1])
#             THLs.append(int(temp.split("_")[0]))
#             f = open(assembly_path+"/"+filename,'r')
#             matrix = [ map(int,line.split(' ')) for line in f ]
#             count=sum(sum(matrix, []))
#             # count=np.count_nonzero(matrix)/(256.0*256.0)
#             counts.append(count)
            
#     canvSCurve=TCanvas("canv", "canv")
#     g = TGraph(len(THLs), array('d', THLs), array('d', counts))
#     g.SetMarkerStyle(25)
#     g.SetMarkerSize(0.8)
#     g.SetMarkerColor(colors[0])
#     g.SetLineColor(colors[0])
    
#     mg.Add(g)
#     leg.AddEntry(g, "Fe55", "l")

#     mg.SetTitle("A06-W0110; THL; Counts");
#     mg.Draw("AP")
#     leg.Draw()
#     gPad.SetGrid(1,1)
#     gPad.Update()   

 

#     # deriv=derivative(THLs, counts, kRed)
#     # canvDeriv=TCanvas("deriv", "deriv")
#     # deriv.Draw("AP")

#     # gPad.SetGrid(1,1)

#     bla=raw_input()

