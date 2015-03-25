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

class CalibTreeMaker:
  
    last_time = time.time()
    AllClusters = []

    def __init__(self,filename,outfile):
        '''
        Constructor
        '''
        self.ReadFile(filename,outfile)
        
    def ReadFile(self,filename,outfile):
        
        
        data_file = open(filename,"r")
        lines = data_file.readlines()
        
        X = []
        Y = []
        TOT = []
        nFrames = 0
        size_cluster = 0
        sizes = []
        

        pixel_1hit = 0
        pixel_2hit = 0
        pixel_3hit = 0
        pixel_4hit = 0
        
       
        hist0 = TH1F("CuIn", "", 6, 0, 6) 
        hist1 = TH2F("CuIn_A06-W0110", "", 256, 0, 256, 256, 0, 256) 
        hist2 = TH2F("CuIn_A06-W0110-filtered", "", 256, 0, 256, 256, 0, 256) 
        hist3 = TH1F("hits per frame","", 50, 0, 500)
       
        
        size1cnt = 0
        allcnt = 0
        everycnt = 0
        self.last_time=time.time()
        
        outfile = TFile(outfile,'recreate')
        singleHits = TTree("singleHits","single hit Clusters")
        allHits = TTree("allHits", "all hit Clusters")
        everyPix = TTree("everyPix", "every single Pixel")
        

        xt=array( 'i', [ 0 ] )
        yt=array( 'i', [ 0 ] )
        tott=array( 'i', [ 0 ])
        csize=array( 'i', [ 0 ])
        
        singleHits.Branch( 'col', xt, 'col/I' )
        singleHits.Branch( 'row', yt, 'row/I' )      
        singleHits.Branch( 'tot', tott, 'tot/I' )
        singleHits.Branch( 'clustersize', csize, 'size/I')
        
        xt_all=array( 'i', [ 0 ] )
        yt_all=array( 'i', [ 0 ] )
        tott_all=array( 'i', [ 0 ])
        csize_all=array( 'i', [ 0 ])

        allHits.Branch( 'col', xt_all, 'col/I' )
        allHits.Branch( 'row', yt_all, 'row/I' )      
        allHits.Branch( 'tot', tott_all, 'tot/I' )
        allHits.Branch( 'clustersize', csize_all, 'size/I')

        xt_every=array( 'i', [ 0 ] )
        yt_every=array( 'i', [ 0 ] )
        tott_every=array( 'i', [ 0 ])
        
        everyPix.Branch( 'col', xt_every, 'col/I' )
        everyPix.Branch( 'row', yt_every, 'row/I' )      
        everyPix.Branch( 'tot', tott_every, 'tot/I' )
        
        do_nothing = 0
        linecnt = 0
        self.last_time = time.time()
        for line in lines : 
            linecnt += 1
           

            if "#" in line : 
                hist3.Fill(linecnt)
                linecnt = 0
                nFrames+=1

                
                
                singles, allClusters = self.MakeClusters(X,Y,TOT,hist0, hist1,hist2)
                everyPixel = self.AllPixels(X,Y,TOT)
               
                for cluster in singles : 
                    
                    size1cnt+=1 
                    
                    xt[0]=cluster[0]
                    yt[0]=cluster[1]
                    tott[0]=cluster[2]
                    csize[0]=cluster[3]
                    singleHits.Fill() 
                
                for cluster in allClusters : 
                    
                    allcnt += 1
                    xt_all[0]=cluster[0]
                    yt_all[0]=cluster[1]
                    tott_all[0]=cluster[2]
                    csize_all[0]=cluster[3]
                    allHits.Fill()   
                    
                for Pixel in everyPixel :
                    
                    everycnt += 1
                    xt_every[0] = Pixel[0]
                    yt_every[0] = Pixel[1]
                    tott_every[0] = Pixel[2]
                    everyPix.Fill()
                 
                nclusters = self.Counting_Clusters(X,Y,TOT)
                pixel_1hit += nclusters.count(1)
                pixel_2hit += nclusters.count(2)
                pixel_3hit += nclusters.count(3)
                pixel_4hit += nclusters.count(4)
 
                               
                                    
                if(nFrames%1000==0):
                    print "Processed Frame %i (%.5fs/frame)"%(nFrames,(time.time()-self.last_time )/1000.)
                    self.last_time = time.time()
                X = []
                Y = []
                TOT = []
                
            
            else : 
                data = line.split()
                X.append(int(data[0]))
                Y.append(int(data[1]))
                TOT.append(int(data[2]))
                linecnt += 1
                del data



        
        del lines 
        data_file.close()
        
       
        c = TCanvas("c", "", 0, 0, 800, 800)
        gPad.SetGrid()
        hist0.Draw('')
        hist0.SetLineColor(kRed)
        hist0.SetLineWidth(2) 
        hist0.GetXaxis().SetTitle("Clustersize")
        hist0.GetYaxis().SetTitle("Events")
        hist0.GetYaxis().SetTitleOffset(1.6)
        c.SaveAs('CuIn_clustersize.pdf')


        c1 = TCanvas()
        hist3.Draw()
        c1.SaveAs('CuIn_hitspframe.pdf')
        #save  histogram
        #gPad.SetGrid()
        #hist1.Draw('col')
        #c.SaveAs('HitMap_CuIn_A06-W0110-before.pdf')

        #save  histogram    
        #gPad.SetGrid()
        #hist2.Draw('col')
        #c.SaveAs('HitMap_Am241_A06-W0110-after.pdf')
        
        singleHits.Write() 
        allHits.Write()
        everyPix.Write()
        outfile.Close()
        

        print "fount a total of", everycnt,"pixels and", allcnt,"clusters"
        print "found %i single pixel clusters"%size1cnt
        print "Clusters:", pixel_1hit, "one hit clusters", pixel_2hit, "two hit clusters", pixel_3hit, "three hit clusters", pixel_4hit, "four hit clusters"        
        

    def AllPixels(self, col, row, tot):
        
        pixels = [[col[i],row[i]] for i,x in enumerate(col)]
       
        everyPixel = []
       
        if(len(pixels)>1):
            
            for hit, TOT in zip(pixels, tot) :
                
                everyPixel.append([hit[0], hit[1], TOT])
                

        else :
            everyPixel = [[pixels[0][0], pixels[0][1], tot[0]]]

        return everyPixel



    def MakeClusters(self,col,row,tot, hist0,hist1, hist2):           
          
        oneHitClusters = []
        allHitClusters = []  
        pixels = [[col[i],row[i]] for i,x in enumerate(col)]
        
        #for pixel in pixels :
            #hist1.Fill(pixel[0],pixel[1])
        
            
        
        if(len(pixels)>1):
        
            results=fclusterdata(pixels,sqrt(2.),criterion="distance",method="single") 
                       
            

            y = numpy.bincount(results)
            ii = numpy.nonzero(y)[0]
       
            do_nothing = 0
            j = 0
            previous = 0
        
            for result, hit, TOT in zip(results,pixels,tot) : 
                #hist2.Fill(hit[0],hit[1])                
                
                tot_c = 0
                i = 0
              
                if y[result]>1:
                    #hist2.Fill(hit[0],hit[1])
                                   
                    if previous != result :  
                       
                        
                        while i <= y[result]-1:
                            if j < len(results) :
                              
                                tot_c += tot[j]
                              
                                j+=1
                                i+=1
                               
                            if j == len(results) :
                                 break
                            
                        if tot_c !=0 :
                           
                            allHitClusters.append([hit[0], hit[1], tot_c, y[result]]) 
                            hist0.Fill(y[result])
                            
                            
                           
                        
                    
                previous = result       
             
                
                if y[result]==1:
                    if j < len(results) :
                        
                        hist0.Fill(y[result])
                        oneHitClusters.append([hit[0], hit[1], TOT, y[result]])
                        allHitClusters.append([hit[0], hit[1], TOT, y[result]])
                        
                        j+=1 
                                         
                   
                    
        if len(pixels)==1 :
            oneHitClusters = [[pixels[0][0],pixels[0][1],tot[0], 1]]       
               
        

        return oneHitClusters, allHitClusters

    def Counting_Clusters(self,col,row,tot):
        sizes = []
        pixels = [[col[i],row[i]] for i,x in enumerate(col)]
              
        if(len(pixels)>1):
        
            results=fclusterdata(pixels,sqrt(2.),criterion="distance",method="single")                
           
            y = numpy.bincount(results)
            ii = numpy.nonzero(y)[0]
       
            j = 0
            previous = 0
        
            for result, hit, TOT in zip(results,pixels,tot) : 
                
                i = 0
              
                if y[result]>1:
                    
                    if previous != result :  
                                               
                        while i <= y[result]-1:
                            if j < len(results) :
                              
                                j+=1
                                i+=1
                               
                            if j == len(results) :
                                 break
                            
                        
                        else :
                             sizes.append(y[result])
                            
                           
                        
                    
                previous = result       
             
                
                if y[result]==1:
                    if j < len(results) :
                        sizes.append(y[result])
                        j+=1 
                                         
                   
                    
        else :
            oneHitClusters = [[pixels[0][0],pixels[0][1],tot[0]]]
        
              
        return sizes


aCalibDataSet = CalibTreeMaker("CuIn_A06-W0110_mine.txt","CuIn_A06-W0110_mine.root")


