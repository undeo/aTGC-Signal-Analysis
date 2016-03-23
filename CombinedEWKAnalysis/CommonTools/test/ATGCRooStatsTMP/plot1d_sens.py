from ROOT import  *
from array import array
from optparse import OptionParser
import math as math

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf
	


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
parser	= OptionParser()
parser.add_option('--POI',dest='POI',help='parameter of interest')
parser.add_option('--pval',dest='pval',help='value of parameter')
parser.add_option('-c',action='store_true',dest='close',default=False)
(options,args) = parser.parse_args()

POI	= options.POI
pval	= options.pval


def plots():
	path		= '../'
	wsname		= 'higgsCombine1Par_%s%s.MultiDimFit.mH120.root'%(POI,pval)
	print wsname
	fileInATGC	= TFile.Open(path+wsname)
	tree		= fileInATGC.Get('limit')
	NEntries	= tree.GetEntries()


	tree.GetEntry(1)	
	if 'cwww' in wsname:
		par		= 'cwww'
		binlo		= int(tree.cwww%-100)-1
		tree.GetEntry(NEntries-1)
		binhi		= int(tree.cwww%100)+1
	elif 'ccw' in wsname:
		par		= 'ccw'
		binlo		= int(tree.ccw%-100)-1
		tree.GetEntry(NEntries-1)
		binhi		= int(tree.ccw%100)+1
	elif 'cb' in wsname:
		par		= 'cb'
		binlo		= int(tree.cb%-100)-1
		tree.GetEntry(NEntries-1)
		binhi		= int(tree.cb%100)+1


	hist		= TH1F('hist','hist',NEntries,binlo,binhi)
	

	for i in range(NEntries):
		if i%1000==0:
			print i
		tree.GetEntry(i+1)
		if 2*tree.deltaNLL < 4.2:
			hist.SetBinContent(i,2*tree.deltaNLL)
		else:
			hist.SetBinContent(i,-1000)

	line4		= TF1('line4','3.84',binlo-1,binhi+1)
	line4.SetLineStyle(kDashed)
	line4.SetLineWidth(4)
	c1		= TCanvas('c1','c1',1)
	#hist.SetMarkerStyle(2)
	#hist.SetMarkerSize(0.3)
	hist.SetMarkerStyle(8)
	hist.SetMarkerSize(0.5)
	#hist.SetMarkerStyle(1)
	if par == 'cwww':
		hist.GetXaxis().SetTitle('c_{WWW} / \Lambda ^2 (1/ TeV ^2)')
	if par == 'ccw':
		hist.GetXaxis().SetTitle('c_{W} / \Lambda ^2 (1/ TeV ^2)')
	if par == 'cb':
		hist.GetXaxis().SetTitle('c_{B} / \Lambda ^2 (1/ TeV ^2)')



	hist.GetYaxis().SetTitle('2*deltaNLL')
	#hist.GetYaxis().SetRangeUser(0,4.2)
	hist.SetMinimum(-0.1)
	if options.close:
		hist.GetXaxis().SetRangeUser(8,12)
		hist.SetMinimum(-0.05)

	hist.Draw('p')
	line4.Draw('SAME')
	c1.Update()

	raw_input('<>')



plots()

