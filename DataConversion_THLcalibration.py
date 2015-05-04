import sys

from ROOT import *
import ROOT
from array import array
import math
import os
import numpy as np
from scipy.cluster.hierarchy import fclusterdata
from os import environ
from itertools import product # SimpleClustering

def SimpleClustering(matrix):

    frame=matrix

    for i,j in [[i,j] for i,j in product(xrange(256),xrange(256)) if matrix[i][j]!=0] :
        for u,v in [[u,v] for u,v in product([-1,0,1],[-1,0,1]) if (((i+u>=0 and i+u<=255) and (j+v>=0 and j+v<=255)) and (u!=0 or v!=0))]:
            if(matrix[i+u][j+v]!=0):
                frame[i][j]=0
                frame[i+u][j+v]=0

    return sum(sum(frame, []))


def writeFile(THL, Counts, OutputFile):
    txtfile = open(OutputFile, "w")
    
    for i in range (0, len(THL)):
        txtfile.write(str(THL[i])+" "+str(Counts[i]))
        
    txtfile.close()
    

def usage():
    print 'THL calibration'
    print 'Usage:\n python %s <assembly> <frames> <THLmin> <THLmax> <source>' % ( sys.argv[0] )


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


    home = environ['HOME']
    path=home+"/eos/clicdp/data/VertexCalibration/ThresholdScan/"
    directory=path+assembly+"/"
    path=directory+source
    THLstep=2




    THLs=[]
    counts=[]

    txtfile = open( assembly+"_"+source+".txt", "w")

    for thl in range(THLmin, THLmax, THLstep):
        
        THLs.append(thl)
        count=0
        for frame in range(0, frames):
            filename="%s_THL%03i_%d"%(source, thl, frame)
            f = open(path+"/"+filename,'r')
            matrix = [ map(int,line.split(' ')) for line in f ]

            # print "Occupancy: %d pixels"%(np.count_nonzero(matrix))

            # ============ Without clustering ============ #
            count+=sum(sum(matrix, []))
            # ============ End without clustering ============ #

            # # ====== fclusterdata used for clustering ====== #
            # coords=np.array(np.transpose(np.nonzero(matrix))) # coordinate of non-zero elements
            # row=[]
            # col=[]
            # countForThisFrame=[]

            # for i in range(0, len(coords)):
            #     x=coords[i][0]
            #     y=coords[i][1]
            #     row.append(x)
            #     col.append(y)
            #     countForThisFrame.append(matrix[x][y])



            # ClusteringResult=fclusterdata(coords, sqrt(2.), criterion="distance") #clustering
            # NbOfOccurance=np.bincount(ClusteringResult) # occurance of each value
            # for i in range(1, len(NbOfOccurance)):
            #     if (NbOfOccurance[i]!=1):
            #         multiHitClusters=np.where(ClusteringResult==i)
            #         for j in range(0, len(multiHitClusters)):
            #             row.pop(j)
            #             col.pop(j)
            #             countForThisFrame.pop(j)

            # count+=sum(countForThisFrame)
            # # ============ End fclusterdata ============== #

            # # ===== using simple clustering ==== #
            # count+=SimpleClustering(matrix)
            # # ============ End using simple clustering ==== #

        counts.append(count)
        txtfile.write(str(thl)+" "+str(count)+"\n")
        print "Processing THL="+str(thl)+", Counts="+str(count)

    txtfile.close()
