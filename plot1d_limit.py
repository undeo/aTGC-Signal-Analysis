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
parser.add_option('--mode',dest='mode', default='linter', help='std, lin or linter')
parser.add_option('--obs', action='store_true', default=False)
(options,args) = parser.parse_args()

POI	= options.POI
pval	= 0
par_latex	= {'cwww' : 'c_{WWW} / #Lambda^{2} (TeV^{-2})', 'ccw' : 'c_{W} / #Lambda^{2} (TeV^{-2})', 'cb' : 'c_{B} / #Lambda^{2} (TeV^{-2})', 'dkz' : '#Delta#kappa_{Z}', 'dg1z' : '#Deltag_{1}^{Z}', 'lZ' : '#lambda_{Z}'}


def plots():
	#wsname		= 'higgsCombine1Par_%s_%s%s.MultiDimFit.mH120.root'%(options.mode,POI,pval)
	wsname_obs	= 'higgsCombine1Par_%s0_obs.MultiDimFit.mH120.root'%POI
	wsname_exp	= 'higgsCombine1Par_%s0_exp.MultiDimFit.mH120.root'%POI
	#wsname_obs	= 'higgsCombine1Par_%s0_obs_WW.MultiDimFit.mH120.root'%POI
	#wsname_exp	= 'higgsCombine1Par_%s0_exp_WW.MultiDimFit.mH120.root'%POI

	print 'reading ' + wsname_obs + ' for observed limits'
	print 'reading ' + wsname_exp + ' for expected limits'

	fileInATGC_obs	= TFile.Open(wsname_obs)
	tree_obs	= fileInATGC_obs.Get('limit')
	NEntries_obs	= tree_obs.GetEntries()
	fileInATGC_exp	= TFile.Open(wsname_exp)
	tree_exp	= fileInATGC_exp.Get('limit')
	NEntries_exp	= tree_exp.GetEntries()


	tree_exp.GetEntry(1)	

	par		= POI
	
	x_exp	= []
	y_exp	= []

	for i in range(NEntries_exp-1):
		if i%1000==0:
			print i
		tree_exp.GetEntry(i+1)
		if 2*tree_exp.deltaNLL < 4.5:
			x_exp.append(tree_exp.GetLeaf(par).GetValue())
			y_exp.append(2*tree_exp.deltaNLL)

	graph_exp	= TGraph(len(x_exp),array('d',x_exp),array('d',y_exp))

	line4		= TF1('line4','3.84',graph_exp.GetXaxis().GetXmin(),graph_exp.GetXaxis().GetXmax())
	line4.SetLineStyle(1)
	line4.SetLineColor(kBlack)
	line4.SetLineWidth(1)
	c1		= TCanvas()
	#c1.SetRightMargin(0.2)

	graph_exp.GetXaxis().SetTitle(par_latex[par])

	for i in range(NEntries_exp/2):
		j=i+1
		tree_exp.GetEntry(j)
		if 2*tree_exp.deltaNLL>3.84 and i<NEntries_exp/2-1:
			continue
		limlo_exp = tree_exp.GetLeaf(par).GetValue()
		break
	for i in range(NEntries_exp/2):
		j=i+ NEntries_exp/2 +1
		tree_exp.GetEntry(j)
		if 2*tree_exp.deltaNLL<3.84 and i<NEntries_exp/2-1:
			continue
		limhi_exp = tree_exp.GetLeaf(par).GetValue()
		break

	linelo_exp = TLine(limlo_exp,0,limlo_exp,4.5)
	linehi_exp = TLine(limhi_exp,0,limhi_exp,4.5)
	linelo_exp.SetLineStyle(7)
	linehi_exp.SetLineStyle(7)
	graph_exp.GetYaxis().SetTitle('2#DeltaNLL')
	graph_exp.GetYaxis().SetTitleSize(0.05)
	graph_exp.GetYaxis().SetTitleOffset(0.75)
	graph_exp.GetYaxis().SetRangeUser(0,4.5)
	graph_exp.GetXaxis().SetTitleSize(0.05)
	graph_exp.GetXaxis().SetTitleOffset(0.75)
	graph_exp.SetLineWidth(2)
	graph_exp.SetLineStyle(7)
	graph_exp.Draw()
	line4.Draw('SAME')
	linelo_exp.Draw('SAME')
	linehi_exp.Draw('SAME')



	tree_obs.GetEntry(1)	

	par		= POI
	
	x_obs	= []
	y_obs	= []

	for i in range(NEntries_obs-1):
		if i%1000==0:
			print i
		tree_obs.GetEntry(i+1)
		if 2*tree_obs.deltaNLL < 4.5:
			x_obs.append(tree_obs.GetLeaf(par).GetValue())
			y_obs.append(2*tree_obs.deltaNLL)

	graph_obs	= TGraph(len(x_obs),array('d',x_obs),array('d',y_obs))


	for i in range(NEntries_obs/2):
		j=i+1
		tree_obs.GetEntry(j)
		if 2*tree_obs.deltaNLL>3.84 and i<NEntries_obs/2-1:
			continue
		limlo_obs = tree_obs.GetLeaf(par).GetValue()
		break
	for i in range(NEntries_obs/2):
		j=i+ NEntries_obs/2 +1
		tree_obs.GetEntry(j)
		if 2*tree_obs.deltaNLL<3.84 and i<NEntries_obs/2-1:
			continue
		limhi_obs = tree_obs.GetLeaf(par).GetValue()
		break

	linelo_obs = TLine(limlo_obs,0,limlo_obs,4.5)
	linehi_obs = TLine(limhi_obs,0,limhi_obs,4.5)
	linelo_obs.SetLineStyle(1)
	linehi_obs.SetLineStyle(1)
	linelo_obs.SetLineColor(kBlue)
	linehi_obs.SetLineColor(kBlue)
	graph_obs.SetLineStyle(1)
	graph_obs.SetLineWidth(2)
	graph_obs.SetLineColor(kBlue)
	graph_obs.Draw('SAME')
	linelo_obs.Draw('SAME')
	linehi_obs.Draw('SAME')

	#leg	= TLegend(0.82,0.5,0.89,0.95)
	#leg.SetFillColor(0)
	#leg.SetBorderSize(0)
	#leg.AddEntry(graph_exp,'expected limit 95% C.L.','l')
	#leg.Draw("SAME")


	c1.Update()
	c1.SaveAs("ATGCRooStatsTMP/docuplots/1dlimit_%s.pdf"%POI)
	c1.SaveAs("ATGCRooStatsTMP/docuplots/1dlimit_%s.png"%POI)
	error = (float(x_obs[1]) - float(x_obs[0]))

	if par in ['cwww','ccw','cb']:
		print '95%% C.L. obs. limit on %s: [%s,%s] +- %s'%(par,round(limlo_obs,2),round(limhi_obs,2),round(error,2))
	else:
		print '95%% C.L. obs. limit on %s: [%s,%s] +- %s'%(par,round(limlo_obs,4),round(limhi_obs,4),round(error,4))

	raw_input('<>')



plots()

