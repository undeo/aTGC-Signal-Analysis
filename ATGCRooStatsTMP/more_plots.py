from ROOT import  *
from array import array
from optparse import OptionParser
import math as math

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf
	
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)


def plots(wsname,binlo,binhi,nbins):
	path		= '../CombinedEWKAnalysis/CommonTools/test/'
	


	fileInATGC	= TFile.Open(path+wsname)
	treeATGC	= fileInATGC.Get('limit')
	hist		= TH1F('hist','hist',nbins,binlo,binhi)

	for i in range(treeATGC.GetEntries()):
		if i%1000==0:
			print i
		treeATGC.GetEntry(i)
		if 2*limit.deltaNLL < 4.2:
			hist.SetBinContent(i,2*limit.deltaNLL)
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
	if 'cwww' in wsname:
		hist.GetXaxis().SetTitle('c_{WWW}')
	if 'ccw' in wsname:
		hist.GetXaxis().SetTitle('c_{W}')
	if 'cb' in wsname:
		hist.GetXaxis().SetTitle('c_{B}')
	hist.GetYaxis().SetTitle('2*deltaNLL')
	#hist.GetYaxis().SetRangeUser(0,4.2)
	#hist.GetXaxis().SetRangeUser(10,13)
	hist.SetMinimum(0)
	hist.Draw('p')
	line4.Draw('SAME')
	c1.Update()

	raw_input('<>')


def plot2d(wsname):
	path		= '../CombinedEWKAnalysis/CommonTools/test/'
	fileInATGC	= TFile.Open(path+wsname)
	treeATGC	= fileInATGC.Get('limit')
	hist		= TH2F('hist','hist',315,-18,18,315,-30,30)

	j = 0
	for i in range(treeATGC.GetEntries()):
		if i%1000==0:
			print i
		treeATGC.GetEntry(i)
		if 2*limit.deltaNLL<5:
			hist.SetBinContent(i%315,j,2*limit.deltaNLL)
		else:
			hist.SetBinContent(i%315,j,-1000)
		if i%315==0:
			j+=1
	c1 = TCanvas('c1','c1',1)
	hist.SetMinimum(0)
	hist.GetXaxis().SetTitle('c_{WWW}')
	hist.GetYaxis().SetTitle('c_{W}')
	hist.GetZaxis().SetTitle('2*deltaNLL')
	hist.Draw('colz prof')
	c1.Update()
	raw_input('</>')


plots('higgsCombine1Par_cwww1124.MultiDimFit.mH120.root',-20,20,1000)
#plots('higgsCombine1Par_ccw16.MultiDimFit.mH120.root',-30,30,10000)
#plots('higgsCombine1Par_cb42_4.MultiDimFit.mH120.root',-75,75,10000)
#plot2d('higgsCombine2Par_cwww5_ccw8.MultiDimFit.mH120.root')
#plot2d('higgsCombine2Par_cwww1124_ccw0.MultiDimFit.mH120.root')
