# will use kde method to finding most likely TOT
# for the global calibration

from optparse import OptionParser
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
import ROOT as R
from os import environ
import getpass
import Constants as C

parser = OptionParser()
parser.add_option("-b", "--assembly",
                  help="Assembly name", dest="ASSEMBLY")

parser.add_option("-p", "--peak",
                  help="Peak name", dest="PEAK")

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

if(options.PEAK):
    peak=options.PEAK
else :
    print "Please specify peak"
    print "choose from", C.known_peaks
    parser.print_help()
    exit()

if peak not in C.known_peaks:
    print "Peak not recognised"
    print "choose from", C.known_peaks
    exit()


def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)


def findMostLikelyTOT(assembly,peak,llim,ulim):

    # get data    
    home = environ['HOME']
    base = "%s/eos/clicdp/data/VertexCalibration" %home
    assembly_start = assembly.split("-")[0]

    if peak == "Fe":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Fe55_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif peak == "Am2" or peak == "Am3":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Am241_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif peak == "Cd":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cd109_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif peak == "Cu" or peak == "In":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cu_In_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/CuIn_%s.root" %(base,assembly,assembly))
    elif peak == "Co1" or peak == "Co2":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Co57_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/Co57_%s.root" %(base,assembly,assembly))

    tree = rootfile.Get("pixels")
    print "got tree"

    username = getpass.getuser()
    junkfile = R.TFile("/tmp/%s/junkfile_%s_%s.root" %(username,assembly,peak),"RECREATE")

    t2 = tree.CopyTree("tot < %i && tot > %i" %(ulim,llim))
    print "copied tree"

    x_grid = np.linspace(llim, ulim, 100)
    ent = t2.GetEntries()
    print "Will load", ent, "entries"
    tot = []
    if ent > 1:
        tot = np.zeros(ent, dtype=float)
        t2.Branch('tot', tot, 'tot/F')
        for i in xrange(ent):
            if i%1000000==0:
                print ".....loading", i
            t2.GetEvent(i)
            tot[i] = t2.tot

        # doesn't work wth my scipy version
        # pdf = kde_scipy(x, x_grid, bandwidth=0.2)
        # workaround
        try:
            print "trying"
            density = gaussian_kde(tot)
            density.covariance_factor = density.silverman_factor
            density._compute_covariance()
            workaround = density(x_grid)

            maxindex = workaround.argmax()
            maxtot = x_grid[maxindex]

            lowerindex = maxindex-1
            while (np.trapz(workaround[lowerindex:maxindex],x=x_grid[lowerindex:maxindex]) < 0.341*2*np.trapz(workaround[0:maxindex],x=x_grid[0:maxindex])) and (lowerindex>0):
                lowerindex=lowerindex-1
            lowersigma = (x_grid[maxindex]-x_grid[lowerindex]) / R.sqrt(ent)

            upperindex = maxindex+1
            while (np.trapz(workaround[maxindex:upperindex],x=x_grid[maxindex:upperindex]) < 0.341*2*np.trapz(workaround[maxindex:99],x=x_grid[maxindex:99])) and (upperindex<99):
                upperindex=upperindex+1
            uppersigma = (x_grid[upperindex]-x_grid[maxindex]) / R.sqrt(ent)

        except np.linalg.linalg.LinAlgError as err:
            maxtot = 0.
            lowersigma = 10000.
            uppersigma = 10000.

    else:
        maxtot = 0.
        lowersigma = 10000.
        uppersigma = 10000.

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    plt.xticks(np.arange(llim,ulim+1,(ulim-llim)/5.))
    ax.set_xlabel('TOT (ADC)')
    ax.text(0.01, 0.99, r'Max: $%i \pm ^{%0.2f} _{%0.2f}$' %(maxtot,uppersigma,lowersigma),
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes,
            fontsize=40)
    ax.hist(tot, bins=100,fc='gray',alpha=0.3,normed=True)
    ax.plot(x_grid, workaround, color='blue', alpha=0.5, lw=3)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    fig.savefig("plots/KDEPeaks/%s_%s_GlobalSpectrum.pdf" %(assembly,peak))

    print "assembly", assembly, "maxtot", maxtot, "pm", uppersigma, lowersigma

    # txt file
    f = open('results/kde/%s_%s_GlobalResults.txt' %(assembly,peak), 'w')
    f.write('%f \t %f \t %f \n' %(maxtot,lowersigma,uppersigma))
    f.close()
    print "finished", assembly, peak


def getLimits(assembly,peak):

    if assembly == "A06-W0110":
        if peak == "Fe": limits = [0,400]
        if peak == "Am3": limits = [900,1500]
        if peak == "Am2": limits = [500,900]
        if peak == "Cd": limits = [400,800]
        if peak == "Cu": limits = [0,500]
        if peak == "In": limits = [500,1000]
        if peak == "Co1": limits = [0,300]
        if peak == "Co2": limits = [350,700]

    if assembly == "B06-W0125":
        if peak == "Fe": limits = [0,700]
        if peak == "Am3": limits = [1800,2800]
        if peak == "Am2": limits = [900,1800]
        if peak == "Cd": limits = [800,1700]
        if peak == "Cu": limits = [0,700]
        if peak == "In": limits = [800,1700]
        if peak == "Co1": limits = [0,600]
        if peak == "Co2": limits = [600,1100]

    if assembly == "B07-W0125":
        if peak == "Fe": limits = [0,500]
        if peak == "Am3": limits = [900,1500]
        if peak == "Am2": limits = [550,1000]
        if peak == "Cd": limits = [500,900]
        if peak == "Cu": limits = [0,400]
        if peak == "In": limits = [400,1000]
        if peak == "Co1": limits = [0,400]
        if peak == "Co2": limits = [350,600]

    if assembly == "C04-W0110":
        if peak == "Fe": limits = [0,500]
        if peak == "Am3": limits = [900,1400]
        if peak == "Am2": limits = [500,900]
        if peak == "Cd": limits = [400,800]
        if peak == "Cu": limits = [0,500]
        if peak == "In": limits = [500,700]
        if peak == "Co1": limits = [0,400]
        if peak == "Co2": limits = [350,600]

    if assembly == "D09-W0126":
        if peak == "Fe": limits = [0,600]
        if peak == "Am3": limits = [1200,1800]
        if peak == "Am2": limits = [600,1200]
        if peak == "Cd": limits = [600,1200]
        if peak == "Cu": limits = [0,600]
        if peak == "In": limits = [600,1100]
        if peak == "Co1": limits = [0,500]
        if peak == "Co2": limits = [500,800]

    if assembly == "L04-W0125":
        if peak == "Fe": limits = [0,500]
        if peak == "Am3": limits = [1000,2000]
        if peak == "Am2": limits = [600,1000]
        if peak == "Cd": limits = [500,1100]
        if peak == "Cu": limits = [0,600]
        if peak == "In": limits = [600,1200]
        if peak == "Co1": limits = [0,400]
        if peak == "Co2": limits = [400,700]

    return limits


limits = getLimits(assembly,peak)
findMostLikelyTOT(assembly,peak,limits[0],limits[1])

