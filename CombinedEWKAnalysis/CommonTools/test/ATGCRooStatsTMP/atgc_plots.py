from ROOT import  *
from array import array
import math as math

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf




def make_input(ch = 'ele'):

	nbins4fit	= 20
	binlo		= 1000
	binhi		= 5000


	#read ATGC and new tree



	cwww_pos_hist	= TH1F('cwww_pos_hist','cwww_pos_hist',nbins4fit,binlo,binhi)
	cwww_neg_hist	= TH1F('cwww_neg_hist','cwww_neg_hist',nbins4fit,binlo,binhi)
	ccw_pos_hist	= TH1F('ccw_pos_hist','ccw_pos_hist',nbins4fit,binlo,binhi)
	ccw_neg_hist	= TH1F('ccw_neg_hist','ccw_neg_hist',nbins4fit,binlo,binhi)
	cb_pos_hist	= TH1F('cb_pos_hist','cb_pos_hist',nbins4fit,binlo,binhi)
	cb_neg_hist	= TH1F('cb_neg_hist','cb_neg_hist',nbins4fit,binlo,binhi)
	SMhist		= TH1F('SMhist','SMhist',nbins4fit,binlo,binhi)
	all3_hist	= TH1F('all3_hist','all3_hist',nbins4fit,binlo,binhi)
	for ch in ['ele','mu']:
		fileInATGC	= TFile.Open('Output/ATGC-Tree_%s.root'%ch)
		treeATGC	= fileInATGC.Get('BasicTree')
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
			if i%10000==0:
				print i
		    	if (treeATGC.c_wwwl == 12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0):
				cwww_pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			if (treeATGC.c_wwwl == 0 and treeATGC.c_wl == 20 and treeATGC.c_bl == 0):
				ccw_pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			if (treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == 60):
		      		cb_pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
		    	if (treeATGC.c_wwwl == -12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0):
				cwww_neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			if (treeATGC.c_wwwl == 0 and treeATGC.c_wl == -20 and treeATGC.c_bl == 0):
				ccw_neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			if (treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == -60):
		      		cb_neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			if treeATGC.c_wwwl == 0 and treeATGC.c_wl ==0 and treeATGC.c_bl == 0:
		    		SMhist.Fill(treeATGC.m_lvj,treeATGC.weight)
			if treeATGC.c_wwwl == -12 and treeATGC.c_wl ==-20 and treeATGC.c_bl == -60:
		    		all3_hist.Fill(treeATGC.m_lvj,treeATGC.weight)


	gStyle.SetOptStat(0)
	gStyle.SetOptTitle(0)
	gStyle.SetLineWidth(2)
	#gStyle.SetMarkerSize(0.75)
	canvas		= TCanvas('test','test',1)
	canvas.cd()
	canvas.SetLogy()

	leg = TLegend(0.15,0.11,0.5,0.5)
	leg.SetFillColor(kWhite)
	leg.SetBorderSize(0)

	all3_hist.SetMinimum(8e-4)	
	all3_hist.SetMaximum(500)
	all3_hist.GetXaxis().SetTitle('M_{WW} (GeV)')
	all3_hist.GetYaxis().SetTitle('N_{events}')
	all3_hist.SetMarkerColor(kRed)
	all3_hist.SetLineColor(kRed)
	all3_hist.SetLineWidth(3)
	all3_hist.Sumw2()
	all3_hist.Draw('HIST')
	cwww_neg_hist.SetLineColor(kOrange)
	cwww_neg_hist.Sumw2()
	cwww_neg_hist.SetLineWidth(3)
	cwww_neg_hist.Draw('HIST SAME')
	ccw_neg_hist.SetLineColor(kGreen+1)
	ccw_neg_hist.SetLineWidth(3)
	ccw_neg_hist.Draw('HIST SAME')
	cb_neg_hist.SetLineColor(kBlack)
	cb_neg_hist.SetLineWidth(3)
	cb_neg_hist.Draw('HIST SAME')


	SMhist.SetLineColor(kBlue)
	SMhist.Sumw2()
	SMhist.SetLineWidth(3)
	SMhist.Draw('HIST SAME')


	leg.SetHeader('parameters in 1/TeV^2')
	leg.AddEntry(all3_hist,'c_{WWW}/\Lambda ^2= -12 , c_{W}/\Lambda ^2= -20 , c_{B}/\Lambda ^2= -60 ','l')
	leg.AddEntry(cb_neg_hist,'c_{WWW}/\Lambda ^2= 0 ,\t c_{W}/\Lambda ^2= 0 ,\t c_{B}/\Lambda ^2= -60 ','l')
	leg.AddEntry(ccw_neg_hist,'c_{WWW}/\Lambda ^2= 0 ,\t c_{W}/\Lambda ^2= -20 , c_{B}/\Lambda ^2= 0 ','l')
	leg.AddEntry(cwww_neg_hist,'c_{WWW}/\Lambda ^2= -12 ,\t c_{W}/\Lambda ^2= 0 ,\t c_{B}/\Lambda ^2= 0 ','l')
	leg.AddEntry(SMhist,'c_{WWW}/\Lambda ^2= 0 ,\t c_{W}/\Lambda ^2= 0 , c_{B}/\Lambda ^2= 0 ','l')
	leg.Draw('SAME')

	canvas.Update()
	raw_input("@@")
	

  
make_input('ele')
