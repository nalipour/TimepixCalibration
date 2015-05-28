# will use kde method to finding most likely TOT
# for the pixel calibration

from optparse import OptionParser
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
import ROOT as R
from os import environ
import getpass
import math
import itertools
import Constants as C

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

if source in C.LNLS_sources and assembly != "A06-W0110":
    print "Source only available for assembly A06-W0110"
    print "please reconsider input"
    exit()

# wanted to use kde_scipy function, but wouldn't work wth my scipy version
# this is a workaround which does the same


def findMostLikelyTOT(assembly,source,llim,ulim,CuInXRFMidPoint):

    # Load data
    home = environ['HOME']
    base = "%s/eos/clicdp/data/VertexCalibration" %home
    assembly_start = assembly.split("-")[0]

    if source == "Fe":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Fe55_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "Am":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Am241_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "Cd":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cd109_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "CuInXRF":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cu_In_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/CuIn_%s.root" %(base,assembly,assembly))
    elif source == "Co":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Co57_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/Co57_%s.root" %(base,assembly,assembly))
    elif source == "CoXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_CoXRF_CalibTree.root" %(base,assembly))
    elif source == "CrXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_CrXRF_CalibTree.root" %(base,assembly))
    elif source == "CuXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_CuXRF_CalibTree.root" %(base,assembly))
    elif source == "FeXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_FeXRF_CalibTree.root" %(base,assembly))
    elif source == "MnXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_MnXRF_CalibTree.root" %(base,assembly))
    elif source == "NiXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_NiXRF_CalibTree.root" %(base,assembly))
    elif source == "TiXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_TiXRF_CalibTree.root" %(base,assembly))
    elif source == "VXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_VXRF_CalibTree.root" %(base,assembly))

    tree = rootfile.Get("pixels")
    print "got tree"

    # Set up junk file to appease ROOT
    username = getpass.getuser()
    junkfile = R.TFile("/tmp/%s/junkfile_%s_%s.root" %(username,assembly,source),"RECREATE")

    # Just keep the events we need
    t0 = tree.CopyTree("tot < %i && tot > %i" %(ulim,llim))
    print "copied whole tree"

    nsteps = 100
    x_grid = np.linspace(llim, ulim, nsteps)
    step_size = (ulim-llim)/float(nsteps-1)
    if source in ["Fe","Cd","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]:
        mostlike2d_1 = []
        lowersigmas_1 = []
        uppersigmas_1 = []
    elif source in ["Co","CuInXRF"]:
        mostlike2d_1 = []
        lowersigmas_1 = []
        uppersigmas_1 = []
        mostlike2d_2 = []
        lowersigmas_2 = []
        uppersigmas_2 = []
    elif source in ["Am"]:
        mostlike2d_1 = []
        lowersigmas_1 = []
        uppersigmas_1 = []
        mostlike2d_2 = []
        lowersigmas_2 = []
        uppersigmas_2 = []
        mostlike2d_3 = []
        lowersigmas_3 = []
        uppersigmas_3 = []

    for r in xrange(C.npixX):
        print "row", r

        t2 = t0.CopyTree("row==%i" %r)
        if source in ["Fe","Cd","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]:
            mostlike2d_1.append([])
            lowersigmas_1.append([])
            uppersigmas_1.append([])
        elif source in ["Co","CuInXRF"]:
            mostlike2d_1.append([])
            lowersigmas_1.append([])
            uppersigmas_1.append([])
            mostlike2d_2.append([])
            lowersigmas_2.append([])
            uppersigmas_2.append([])
        elif source in ["Am"]:
            mostlike2d_1.append([])
            lowersigmas_1.append([])
            uppersigmas_1.append([])
            mostlike2d_2.append([])
            lowersigmas_2.append([])
            uppersigmas_2.append([])
            mostlike2d_3.append([])
            lowersigmas_3.append([])
            uppersigmas_3.append([])

        for c in xrange(C.npixY):
            #print "col", c
            peak_tots = []
            peak_amps = []
            peak_is = []
            peak_loweris = []
            peak_upperis = []
            peak_ents = []
            peak_lowersigmas = []
            peak_uppersigmas = []

            t2.SetEstimate(t2.GetEntries())
            ent = t2.GetEntries("col==%i" %c)
            mytot = []
            if ent > 1:
                t2.Draw("tot", "col==%i"%c, "goff")
                v1 = t2.GetV1()
                for i in xrange(ent):
                    mytot.append(v1[i])

                try:
                    # The KDE calculation
                    density = gaussian_kde(mytot)
                    if density.silverman_factor() > 0.1:
                        density.covariance_factor = density.silverman_factor
                    else:
                        density.covariance_factor = lambda: 0.1
                    density._compute_covariance()
                    workaround = density(x_grid)

                    # Find peaks by finding where gradient passes through 0
                    grad = np.gradient(workaround)
                    last_grad = grad[0]
                    for i in xrange(1,len(grad)):
                        this_grad = grad[i]
                        if (last_grad > 0 and this_grad < 0):
                            peak_tots.append((x_grid[i-1] + x_grid[i]) / 2.)
                            peak_amps.append((workaround[i-1] + workaround[i]) / 2.)
                            peak_is.append(i-0.5)
                        last_grad = this_grad
                
                    # Remove peaks found at the edge
                    for i in reversed(xrange(len(peak_tots))):
                        if peak_is[i] < 3 or peak_is[i] > (nsteps-3):
                            peak_tots.pop(i)
                            peak_amps.pop(i)
                            peak_is.pop(i)

                    # Filter peaks to take only highest
                    sorted_peaks = sorted(zip(peak_amps,peak_tots,peak_is),reverse=1)
                    if len(peak_tots) == 0:
                        print "Found no peaks"
                        peak_amps = []
                        peak_tots = []
                        peak_is = []
                    elif source in ["Fe","Cd","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]:
                        peak_amps = [sorted_peaks[0][0]]
                        peak_tots = [sorted_peaks[0][1]]
                        peak_is = [sorted_peaks[0][2]]

                    elif source in ["Co"]:
                        # take the highest as peak 1, next highest to the right as peak 2
                        if len(sorted_peaks)>=2:
                            peak_amps = [sorted_peaks[0][0]]
                            peak_tots = [sorted_peaks[0][1]]
                            peak_is = [sorted_peaks[0][2]]
                            for i in xrange(1,len(sorted_peaks)):
                                if sorted_peaks[i][1] > sorted_peaks[0][1]:
                                    peak_amps.append(sorted_peaks[i][0])
                                    peak_tots.append(sorted_peaks[i][1])
                                    peak_is.append(sorted_peaks[i][2])
                                    break # only take the first
                        else:
                            peak_amps = [sorted_peaks[0][0]]
                            peak_tots = [sorted_peaks[0][1]]
                            peak_is = [sorted_peaks[0][2]]

                    elif source in ["CuInXRF"]:
                        if assembly == "A06-W0110" or assembly == "B07-W0125":
                            # take the two highest
                            if len(sorted_peaks)>=2:
                                peak_amps = [sorted_peaks[0][0],sorted_peaks[1][0]]
                                peak_tots = [sorted_peaks[0][1],sorted_peaks[1][1]]
                                peak_is = [sorted_peaks[0][2],sorted_peaks[1][2]]
                            else:
                                peak_amps = [sorted_peaks[0][0]]
                                peak_tots = [sorted_peaks[0][1]]
                                peak_is = [sorted_peaks[0][2]]
                        elif assembly == "B06-W0125" or assembly == "L04-W0125" or assembly == "D09-W0126":
                            # take the highest below x and the highest after x
                            sorted_peaks_low = []
                            sorted_peaks_high = []
                            for sorted_peak in sorted_peaks:
                                if sorted_peak[1] < CuInXRFMidPoint:
                                    sorted_peaks_low.append(sorted_peak)
                                else:
                                    sorted_peaks_high.append(sorted_peak)
                            if len(sorted_peaks_low) >= 1 and len(sorted_peaks_high) >= 1:
                                peak_amps = [sorted_peaks_low[0][0],sorted_peaks_high[0][0]]
                                peak_tots = [sorted_peaks_low[0][1],sorted_peaks_high[0][1]]
                                peak_is = [sorted_peaks_low[0][2],sorted_peaks_high[0][2]]
                            elif len(sorted_peaks_low) == 0:
                                peak_amps = [sorted_peaks_high[0][0]]
                                peak_tots = [sorted_peaks_high[0][1]]
                                peak_is = [sorted_peaks_high[0][2]]
                            else:
                                peak_amps = [sorted_peaks_low[0][0]]
                                peak_tots = [sorted_peaks_low[0][1]]
                                peak_is = [sorted_peaks_low[0][2]]
                        elif assembly == "C04-W0110":
                            # take the highest as peak 1, next highest to the right as peak 2
                            if len(sorted_peaks)>=2:
                                peak_amps = [sorted_peaks[0][0]]
                                peak_tots = [sorted_peaks[0][1]]
                                peak_is = [sorted_peaks[0][2]]
                                for i in xrange(1,len(sorted_peaks)):
                                    if sorted_peaks[i][1] > sorted_peaks[0][1]:
                                        peak_amps.append(sorted_peaks[i][0])
                                        peak_tots.append(sorted_peaks[i][1])
                                        peak_is.append(sorted_peaks[i][2])
                                        break # only take the first
                            else:
                                peak_amps = [sorted_peaks[0][0]]
                                peak_tots = [sorted_peaks[0][1]]
                                peak_is = [sorted_peaks[0][2]]

                    elif source == "Am":
                        if len(sorted_peaks)>=3:
                            peak_amps = [sorted_peaks[0][0],sorted_peaks[1][0],sorted_peaks[2][0]]
                            peak_tots = [sorted_peaks[0][1],sorted_peaks[1][1],sorted_peaks[2][1]]
                            peak_is = [sorted_peaks[0][2],sorted_peaks[1][2],sorted_peaks[2][2]]
                        elif len(sorted_peaks)==2:
                            peak_amps = [sorted_peaks[0][0],sorted_peaks[1][0]]
                            peak_tots = [sorted_peaks[0][1],sorted_peaks[1][1]]
                            peak_is = [sorted_peaks[0][2],sorted_peaks[1][2]]
                        else:
                            peak_amps = [sorted_peaks[0][0]]
                            peak_tots = [sorted_peaks[0][1]]
                            peak_is = [sorted_peaks[0][2]]

                    # Put peaks back in TOT order
                    ordered_peaks = sorted(zip(peak_tots,peak_amps,peak_is))
                    peak_tots = [peaki[0] for peaki in ordered_peaks]
                    peak_amps = [peaki[1] for peaki in ordered_peaks]
                    peak_is = [peaki[2] for peaki in ordered_peaks]

                    # Calculate the upper and lower i for each peak
                    for i in xrange(len(peak_tots)):
                        if i == 0 and i == len(peak_tots)-1:
                            peak_loweris.append(0)
                            peak_upperis.append(nsteps)
                        elif i == 0:
                            peak_loweris.append(0)
                            peak_upperis.append(int((peak_is[i] + peak_is[i+1])/2.))
                        elif i == len(peak_tots)-1:
                            peak_loweris.append(int((peak_is[i-1] + peak_is[i])/2.))
                            peak_upperis.append(nsteps)
                        else:
                            peak_loweris.append(int((peak_is[i-1] + peak_is[i])/2.))
                            peak_upperis.append(int((peak_is[i] + peak_is[i+1])/2.))

                    # Calculate the entries in each peak
                    for i in xrange(len(peak_tots)):
                        if i == 0 and i == len(peak_tots)-1:
                            peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i && col==%i" %(llim,ulim,c))))
                        elif i == 0:
                            peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i && col==%i" %(llim,(peak_tots[i] + peak_tots[i+1])/2.,c))))
                        elif i == len(peak_tots)-1:
                            peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i && col==%i" %((peak_tots[i-1] + peak_tots[i])/2.,ulim,c))))
                        else:
                            peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i && col==%i" %((peak_tots[i-1] + peak_tots[i])/2.,(peak_tots[i] + peak_tots[i+1])/2.,c))))

                    # Calculate the uncertainty on each peak position
                    for i in xrange(len(peak_tots)):
                        maxindex = peak_is[i]-0.5
                        lowerindex = peak_is[i]-1.5
                        while (np.trapz(workaround[lowerindex:maxindex],x=x_grid[lowerindex:maxindex]) < 0.341*2*np.trapz(workaround[peak_loweris[i]:maxindex],x=x_grid[peak_loweris[i]:maxindex])) and (lowerindex>0):
                            lowerindex=lowerindex-1
                        lowersigma = (x_grid[maxindex]-x_grid[lowerindex]) / R.sqrt(peak_ents[i])
                        peak_lowersigmas.append(lowersigma)

                        maxindex = peak_is[i]+0.5
                        upperindex = peak_is[i]+1.5
                        while (np.trapz(workaround[maxindex:upperindex],x=x_grid[maxindex:upperindex]) < 0.341*2*np.trapz(workaround[maxindex:peak_upperis[i]],x=x_grid[maxindex:peak_upperis[i]])) and (upperindex<(nsteps-1)):
                            upperindex=upperindex+1
                        uppersigma = (x_grid[upperindex]-x_grid[maxindex]) / R.sqrt(peak_ents[i])
                        peak_uppersigmas.append(uppersigma)

                except np.linalg.linalg.LinAlgError as err:
                    peak_tots.append(0) # is this necessary?
                    peak_lowersigmas.append(0)
                    peak_uppersigmas.append(0)

            else:
                peak_tots.append(0)
                peak_lowersigmas.append(0)
                peak_uppersigmas.append(0)

            # Remember these peaks
            if source in ["Am"]:
                if len(peak_tots) == 3:
                    mostlike2d_1[-1].append(peak_tots[0])
                    mostlike2d_2[-1].append(peak_tots[1])
                    mostlike2d_3[-1].append(peak_tots[2])
                    lowersigmas_1[-1].append(peak_lowersigmas[0])
                    lowersigmas_2[-1].append(peak_lowersigmas[1])
                    lowersigmas_3[-1].append(peak_lowersigmas[2])
                    uppersigmas_1[-1].append(peak_uppersigmas[0])
                    uppersigmas_2[-1].append(peak_uppersigmas[1])
                    uppersigmas_3[-1].append(peak_uppersigmas[2])
                elif len(peak_tots) == 2:
                    mostlike2d_1[-1].append(peak_tots[0])
                    mostlike2d_2[-1].append(0)
                    mostlike2d_3[-1].append(peak_tots[1])
                    lowersigmas_1[-1].append(peak_lowersigmas[0])
                    lowersigmas_2[-1].append(0)
                    lowersigmas_3[-1].append(peak_lowersigmas[1])
                    uppersigmas_1[-1].append(peak_uppersigmas[0])
                    uppersigmas_2[-1].append(0)
                    uppersigmas_3[-1].append(peak_uppersigmas[1])
                elif len(peak_tots) == 1:
                    mostlike2d_1[-1].append(peak_tots[0])
                    mostlike2d_2[-1].append(0)
                    mostlike2d_3[-1].append(0)
                    lowersigmas_1[-1].append(peak_lowersigmas[0])
                    lowersigmas_2[-1].append(0)
                    lowersigmas_3[-1].append(0)
                    uppersigmas_1[-1].append(peak_uppersigmas[0])
                    uppersigmas_2[-1].append(0)
                    uppersigmas_3[-1].append(0)
                else:
                    print "=============== problem length", len(peak_tots)
                    print peak_tots

            elif source in ["Co","CuInXRF"]:
                if len(peak_tots) == 2:
                    mostlike2d_1[-1].append(peak_tots[0])
                    mostlike2d_2[-1].append(peak_tots[1])
                    lowersigmas_1[-1].append(peak_lowersigmas[0])
                    lowersigmas_2[-1].append(peak_lowersigmas[1])
                    uppersigmas_1[-1].append(peak_uppersigmas[0])
                    uppersigmas_2[-1].append(peak_uppersigmas[1])
                elif len(peak_tots) == 1:
                    if source == "CuInXRF":
                        if peak_tots[0] < CuInXRFMidPoint:
                            mostlike2d_1[-1].append(peak_tots[0])
                            mostlike2d_2[-1].append(0)
                            lowersigmas_1[-1].append(peak_lowersigmas[0])
                            lowersigmas_2[-1].append(0)
                            uppersigmas_1[-1].append(peak_uppersigmas[0])
                            uppersigmas_2[-1].append(0)
                        else:
                            mostlike2d_1[-1].append(0)
                            mostlike2d_2[-1].append(peak_tots[0])
                            lowersigmas_1[-1].append(0)
                            lowersigmas_2[-1].append(peak_lowersigmas[0])
                            uppersigmas_1[-1].append(0)
                            uppersigmas_2[-1].append(peak_uppersigmas[0])
                    elif source == "Co":
                        mostlike2d_1[-1].append(peak_tots[0])
                        mostlike2d_2[-1].append(0)
                        lowersigmas_1[-1].append(peak_lowersigmas[0])
                        lowersigmas_2[-1].append(0)
                        uppersigmas_1[-1].append(peak_uppersigmas[0])
                        uppersigmas_2[-1].append(0)
                elif len(peak_tots) == 0:
                    mostlike2d_1[-1].append(0)
                    mostlike2d_2[-1].append(0)
                    lowersigmas_1[-1].append(0)
                    lowersigmas_2[-1].append(0)
                    uppersigmas_1[-1].append(0)
                    uppersigmas_2[-1].append(0)
                else:
                    print "=============== problem length", len(peak_tots)
            elif source in ["Fe","Cd","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]:
                if len(peak_tots) == 1:
                    mostlike2d_1[-1].append(peak_tots[0])
                    lowersigmas_1[-1].append(peak_lowersigmas[0])
                    uppersigmas_1[-1].append(peak_uppersigmas[0])
                else:
                    print "=============== problem length", len(peak_tots)


            # Make plots
            if c==0 and len(mytot) and r < 10:
                fig, ax = plt.subplots(1, 1, figsize=(12, 12))
                ax.tick_params(axis='x', pad=20)
                plt.xticks(np.arange(llim,ulim+1,(ulim-llim)/5.))
                ax.set_xlabel('TOT (ADC)')
                for i in xrange(len(peak_tots)):
                    ax.text(0.01, 0.99 - (i*0.1), r'Peak: $%i \pm ^{%0.2f} _{%0.2f} \pm %0.2f$' %(peak_tots[i],peak_uppersigmas[i],peak_lowersigmas[i],step_size/2.),
                            verticalalignment='top', horizontalalignment='left',
                            transform=ax.transAxes,
                            fontsize=40)
                ax.hist(mytot, bins=math.ceil(2.*len(mytot)**(1./3)),fc='gray',alpha=0.3,normed=True) 
                ax.plot(x_grid, workaround, color='blue', alpha=0.5, lw=3)
                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(40)
                fig.tight_layout()
                fig.savefig("plots/KDEPeaks/%s_%s_PixelSpectrum_%i_0.pdf" %(assembly,source,r))

    # 1D plot
    mostlike1d_1 = list(itertools.chain(*mostlike2d_1))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Most likely TOT (ADC)')
    ax.hist(mostlike1d_1, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    fig.savefig("plots/KDEPeaks/%s_%s_hist_1.pdf" %(assembly,source))

    if source in ["Co","CuInXRF","Am"]:
        mostlike1d_2 = list(itertools.chain(*mostlike2d_2))
        fig, ax = plt.subplots(1, 1, figsize=(12, 12))
        ax.tick_params(axis='x', pad=20)
        ax.set_xlabel('Most likely TOT (ADC)')
        ax.hist(mostlike1d_2, 100, fc='gray', alpha=0.3,normed=True)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(40)
        fig.tight_layout()
        limits = ax.axis()
        plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
        fig.savefig("plots/KDEPeaks/%s_%s_hist_2.pdf" %(assembly,source))

    if source in ["Am"]:
        mostlike1d_3 = list(itertools.chain(*mostlike2d_3))
        fig, ax = plt.subplots(1, 1, figsize=(12, 12))
        ax.tick_params(axis='x', pad=20)
        ax.set_xlabel('Most likely TOT (ADC)')
        ax.hist(mostlike1d_3, 100, fc='gray', alpha=0.3,normed=True)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(40)
        fig.tight_layout()
        limits = ax.axis()
        plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
        fig.savefig("plots/KDEPeaks/%s_%s_hist_3.pdf" %(assembly,source))

    dx, dy = 1.0, 1.0
    y, x = np.mgrid[slice(0, C.npixY + dy, dy),slice(0, C.npixX + dx, dx)]

    # txt file
    f = open('results/kde/%s_%s_PixelResults.txt' %(assembly,source), 'w')
    for i in xrange(C.npixX):
        for j in xrange(C.npixY):
            if source in ["Fe","Cd","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]:
                f.write('%f \t %f \t %f \t %f \t %f \t %f \n' %(x[i][j],y[i][j],mostlike2d_1[i][j],lowersigmas_1[i][j],uppersigmas_1[i][j],step_size/2.))
            elif source in ["Co","CuInXRF"]:
                f.write('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n' %(x[i][j],y[i][j],mostlike2d_1[i][j],lowersigmas_1[i][j],uppersigmas_1[i][j],mostlike2d_2[i][j],lowersigmas_2[i][j],uppersigmas_2[i][j],step_size/2.))
            elif source in ["Am"]:
                f.write('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n' %(x[i][j],y[i][j],mostlike2d_2[i][j],lowersigmas_2[i][j],uppersigmas_2[i][j],mostlike2d_3[i][j],lowersigmas_3[i][j],uppersigmas_3[i][j],step_size/2.))
    f.close()

    # Masked 2D plot: no zeros
    mostlike2d_1 = np.ma.masked_equal(mostlike2d_1,0)
    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X')
    ax.set_ylabel('Pixel Y')
    myplot = plt.pcolor(x, y, mostlike2d_1, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.savefig("plots/KDEPeaks/%s_%s_map_nz_1.pdf" %(assembly,source))

    if source in ["Co","CuInXRF","Am"]:
        mostlike2d_2 = np.ma.masked_equal(mostlike2d_2,0)
        fig, ax = plt.subplots(1,1,figsize=(12, 10))
        ax.set_xlabel('Pixel X')
        ax.set_ylabel('Pixel Y')
        myplot = plt.pcolor(x, y, mostlike2d_2, cmap='jet')
        cbar = fig.colorbar(myplot)
        cbar.ax.tick_params(labelsize=40) 
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(40)
        fig.savefig("plots/KDEPeaks/%s_%s_map_nz_2.pdf" %(assembly,source))

    if source in ["Am"]:
        mostlike2d_3 = np.ma.masked_equal(mostlike2d_3,0)
        fig, ax = plt.subplots(1,1,figsize=(12, 10))
        ax.set_xlabel('Pixel X')
        ax.set_ylabel('Pixel Y')
        myplot = plt.pcolor(x, y, mostlike2d_3, cmap='jet')
        cbar = fig.colorbar(myplot)
        cbar.ax.tick_params(labelsize=40) 
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(40)
        fig.savefig("plots/KDEPeaks/%s_%s_map_nz_3.pdf" %(assembly,source))

    print "finished", assembly, source


def getLimits(assembly,source):

    if assembly == "A06-W0110":
        if source == "Fe": limits = [0,400]
        if source == "Am": limits = [0,1500]
        if source == "Cd": limits = [0,800]
        if source == "CuInXRF": limits = [0,1000]
        if source == "Co": limits = [0,700]
        if source == "CoXRF": limits = [0,400]
        if source == "CrXRF": limits = [0,400]
        if source == "CuXRF": limits = [0,400]
        if source == "FeXRF": limits = [0,400]
        if source == "MnXRF": limits = [0,400]
        if source == "NiXRF": limits = [0,400]
        if source == "TiXRF": limits = [0,200]
        if source == "VXRF": limits = [0,400]

    if assembly == "B06-W0125":
        if source == "Fe": limits = [0,700]
        if source == "Am": limits = [0,2800]
        if source == "Cd": limits = [0,1700]
        if source == "CuInXRF": limits = [0,1700]
        if source == "Co": limits = [0,1100]

    if assembly == "B07-W0125":
        if source == "Fe": limits = [0,500]
        if source == "Am": limits = [0,1500]
        if source == "Cd": limits = [500,900]
        if source == "CuInXRF": limits = [0,1000]
        if source == "Co": limits = [0,600]

    if assembly == "C04-W0110":
        if source == "Fe": limits = [0,500]
        if source == "Am": limits = [0,1400]
        if source == "Cd": limits = [0,800]
        if source == "CuInXRF": limits = [0,700]
        if source == "Co": limits = [0,600]

    if assembly == "D09-W0126":
        if source == "Fe": limits = [0,600]
        if source == "Am": limits = [0,1800]
        if source == "Cd": limits = [0,1200]
        if source == "CuInXRF": limits = [0,1100]
        if source == "Co": limits = [0,800]

    if assembly == "L04-W0125":
        if source == "Fe": limits = [0,500]
        if source == "Am": limits = [0,2000]
        if source == "Cd": limits = [0,1100]
        if source == "CuInXRF": limits = [0,1200]
        if source == "Co": limits = [0,700]

    return limits

def getCuInXRFMidPoint(assembly):

    if assembly == "A06-W0110":
        CuInXRFMidPoint = 450
    if assembly == "B06-W0125":
        CuInXRFMidPoint = 700
    if assembly == "B07-W0125":
        CuInXRFMidPoint = 500
    if assembly == "C04-W0110":
        CuInXRFMidPoint = 400
    if assembly == "D09-W0126":
        CuInXRFMidPoint = 600
    if assembly == "L04-W0125":
        CuInXRFMidPoint = 550

    return CuInXRFMidPoint

limits = getLimits(assembly,source)
CuInXRFMidPoint = getCuInXRFMidPoint(assembly)

findMostLikelyTOT(assembly,source,limits[0],limits[1],CuInXRFMidPoint)
