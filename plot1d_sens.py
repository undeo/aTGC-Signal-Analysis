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
parser.add_option('-p', action='store_true', default=False, dest='pos')
(options,args) = parser.parse_args()

POI	= options.POI
pval	= 0

pos = options.pos

def plots():
	if pos:
		wsname		= 'higgsCombine1Par_sens%s_final_pos.MultiDimFit.mH120.root'%POI
	else:
		wsname		= 'higgsCombine1Par_sens%s_final_neg.MultiDimFit.mH120.root'%POI
	print wsname
	fileInATGC	= TFile.Open(wsname)
	tree		= fileInATGC.Get('limit')
	NEntries	= tree.GetEntries()


	tree.GetEntry(1)	

	par		= POI
	
	x	= []
	y	= []

	for i in range(NEntries-1):
		if i%1000==0:
			print i
		tree.GetEntry(i+1)
		if 2*tree.deltaNLL < 4.5:
			x.append(tree.GetLeaf(par).GetValue())
			y.append(2*tree.deltaNLL)

	graph	= TGraph(len(x),array('d',x),array('d',y))

	line4		= TF1('line4','3.84',graph.GetXaxis().GetXmin(),graph.GetXaxis().GetXmax())
	line4.SetLineStyle(7)
	line4.SetLineColor(kBlack)
	line4.SetLineWidth(1)
	c1		= TCanvas()
	if par == 'cwww':
		graph.GetXaxis().SetTitle('c_{WWW} / #Lambda^{2} (TeV^{-2})')
	if par == 'ccw':
		graph.GetXaxis().SetTitle('c_{W} / #Lambda^{2} (TeV^{-2})')
	if par == 'cb':
		graph.GetXaxis().SetTitle('c_{B} / #Lambda^{2} (TeV^{-2})')

	graph.GetYaxis().SetTitle('2#DeltaNLL')
	graph.GetYaxis().SetTitleSize(0.05)
	graph.GetYaxis().SetTitleOffset(0.75)
	graph.SetLineWidth(2)
	graph.Draw()
	line4.Draw('SAME')
	c1.Update()
	if pos:
		c1.SaveAs("ATGCRooStatsTMP/docuplots/1dsens_%s_pos.pdf"%POI)
		c1.SaveAs("ATGCRooStatsTMP/docuplots/1dsens_%s_pos.png"%POI)
	else:
		c1.SaveAs("ATGCRooStatsTMP/docuplots/1dsens_%s_neg.pdf"%POI)
		c1.SaveAs("ATGCRooStatsTMP/docuplots/1dsens_%s_neg.png"%POI)

	error = (float(x[1]) - float(x[0]))


	raw_input('<>')



plots()

