from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import random
import os

gSystem.Load('%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so'%os.environ['CMSSW_BASE'])
from ROOT import RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf



POI		= ['cwww','ccw','cb']
par_max 	= {'cwww' : 12, 'ccw' : 20, 'cb' : 60}#atgc points
par_titles 	= {'cwww' : '#frac{c_{WWW}}{#Lambda^{2}}', 'ccw' : '#frac{c_{W}}{#Lambda^{2}}', 'cb' : '#frac{c_{B}}{#Lambda^{2}}'}#latex titles 

parser	= OptionParser()
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='recreate input trees')
parser.add_option('--n2', action='store_true', dest='newhists', default=False, help='recreate histograms')
parser.add_option('-p', '--plots', action='store_true', dest='make_plots', default=False, help='make plots')
parser.add_option('--savep', action='store_true', dest='savep', default=False, help='save plots')
parser.add_option('--p2', action='store_true', dest='make_plots2', default=False, help='make parabel plot')
parser.add_option('--p3', action='store_true', dest='make_plots3', default=False, help='plot')
parser.add_option('--cat', dest='cat', default='WWWZ', help='category, WW or WZ, defines signal region')
parser.add_option('--std', action='store_true', dest='std', default=False, help='standard mode')
parser.add_option('--linter', action='store_true', dest='linter', default=False, help='include SM- and aTGC-interference')
parser.add_option('--oldfunc', action='store_true', dest='oldfunc', default=False, help='use simple exponential for aTGC')
parser.add_option('-b', action='store_true', dest='batch', default=False, help='batch mode')
parser.add_option('-c', '--ch', dest='chan', default='elmu', help='channel, el, mu or elmu')
parser.add_option('--makeratios', action='store_true', dest='make_ratios', default=False, help='make new WW/WZ ratios')
parser.add_option('--noatgcint', action='store_true', dest='noatgcint', default=False, help='set atgc-interference coefficients to zero')
parser.add_option('--noSMint', action='store_true', dest='noSMint', default=False, help='set SM interference coefficients to zero (after fits!) !')
parser.add_option('--printatgc', action='store_true', default=False, help='print atgc-interference contribution')
parser.add_option('--yieldplots', action='store_true', default=False, help='make plots of relative yields')
parser.add_option('--linfitplot', action='store_true', default=False, help='make plot of SM interference fit')
parser.add_option('--lo', dest='mlvj_lo', default=900)
parser.add_option('--hi', dest='mlvj_hi', default=3500)
parser.add_option('--atgc', action='store_true', dest='atgc', default=False, help='use anomalous coupling parametrization instead of EFT')


(options,args) = parser.parse_args()

signal_category	= options.cat



gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
if options.batch:
	gROOT.SetBatch(True)

if not os.path.isdir('Output'):
	os.system('mkdir Output')
if not os.path.isdir('docuplots'):
	os.system('mkdir docuplots')



def make_ATGCtree(ch='el', signal_cat='WW'):

	if ch=='el':
		METCUT	= 80.
	elif ch=='mu':
		METCUT	= 40.
	else:
		raise RuntimeError('no such channel %s'%ch)


	if signal_cat == 'WW':
		sigreg 	= 'lo'
		mj_lo	= 65.
		mj_hi	= 85.
	elif signal_cat == 'WZ':
		sigreg	= 'hi'
		mj_lo	= 85.
		mj_hi	= 105
	elif signal_cat == 'WV':
		sigreg	= 'all'
		mj_lo	= 65.
		mj_hi	= 105
	else:
		raise RuntimeError('no such category: %s'%signal_cat)

	for categ in ['WW','WZ']:

		print 'reading %s-aTGC_%s.root'%(categ,ch)

		fileInATGC	= TFile.Open('Input/%s-aTGC_%s.root'%(categ,ch))
		treeInATGC	= fileInATGC.Get('treeDumper/BasicTree')
		fileOutATGC	= TFile('Output/ATGC-Tree_%s_%s_%s.root'%(categ,sigreg,ch), 'recreate')
		treeATGC_tmp 	= treeInATGC.CloneTree(0)

		weight		= array('f',[1])
		branch_weight	= treeATGC_tmp.Branch('totEventWeight',weight,'totEventWeight')
		c_wwwl		= array('f',[1])
		branch_c_wwwl	= treeATGC_tmp.Branch('c_wwwl',c_wwwl,'c_wwwl')
		c_wl		= array('f',[1])
		branch_c_wl	= treeATGC_tmp.Branch('c_wl',c_wl,'c_wl')	
		c_bl		= array('f',[1])
		branch_c_bl	= treeATGC_tmp.Branch('c_bl',c_bl,'c_bl')

		lumi_tmp 	= 2300.

		for i in range(treeInATGC.GetEntries()):
			if i%10000==0:
				print i
			treeInATGC.GetEntry(i)
			#apply cuts
			if treeInATGC.jet_pt>200. and treeInATGC.jet_tau2tau1<0.6 \
			and treeInATGC.Mjpruned<mj_hi and treeInATGC.Mjpruned>mj_lo \
			and treeInATGC.W_pt>200. and treeInATGC.deltaR_LeptonWJet>math.pi/2. \
			and abs(treeInATGC.deltaPhi_WJetMet)>2. and abs(treeInATGC.deltaPhi_WJetWlep)>2.\
			and treeInATGC.nbtag==0 and treeInATGC.pfMET>METCUT:
				weight_part	= 1/20. * lumi_tmp * treeInATGC.totWeight
				#SM
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[7] * weight_part
				treeATGC_tmp.Fill()
				#cwww
				c_wwwl[0] = 12; c_wl[0]	= 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[0] * weight_part
				treeATGC_tmp.Fill()
				c_wwwl[0] = -12; c_wl[0] = 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[1] * weight_part
				treeATGC_tmp.Fill()
				#ccw
				c_wwwl[0] = 0; c_wl[0] = 20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[2] * weight_part
				treeATGC_tmp.Fill()
				c_wwwl[0] = 0; c_wl[0] = -20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[3] * weight_part
				treeATGC_tmp.Fill()
				#cb
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = 60;
				weight[0] = treeInATGC.aTGCWeights[4] * weight_part
				treeATGC_tmp.Fill()
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = -60;
				weight[0] = treeInATGC.aTGCWeights[5] * weight_part
				treeATGC_tmp.Fill()
				#all3
				c_wwwl[0] = -12; c_wl[0] = -20; c_bl[0] = -60;
				weight[0] = treeInATGC.aTGCWeights[6] * weight_part
				treeATGC_tmp.Fill()
		treeATGC_tmp.Write()


	fileOut		= TFile('Output/ATGC-Tree_%s_%s.root'%(sigreg,ch), 'recreate')
	mergetree	= TChain('BasicTree')
	#make one tree for each WW and WZ
	mergetree.Add('Output/ATGC-Tree_WW_%s_%s.root'%(sigreg,ch))
	mergetree.Add('Output/ATGC-Tree_WZ_%s_%s.root'%(sigreg,ch))
	newtree		= mergetree.Clone()
	newtree.Write()

	print '--------> Written to file ' + fileOut.GetName()


def make_plots(rrv_x,wtmp,ch,cat,channel,fitres,binlo,binhi,binlo4fit):
        
        can             = []
	can2		= []
        plots	        = []
	plots2		= []
	pads		= []
	
	rrv_x.setVal(2000)

	for i in range(3):
		rrv_x.setRange(binlo4fit,binhi)
                p       = rrv_x.frame(binlo,binhi)
		p2	= rrv_x.frame(binlo,binhi)
                c       = TCanvas(POI[i]+'-',POI[i]+'-',1)
		c.cd()
		pad1	= TPad('pad1_%s'%POI[i],'pad1_%s'%POI[i],0.,0.25,1.,1.)
		pad2	= TPad('pad2_%s'%POI[i],'pad2_%s'%POI[i],0.,0.02,1.,0.25)
		c2	= TCanvas(POI[i]+'+',POI[i]+'+',1)
		c2.cd()
		pad3	= TPad('pad3_%s'%POI[i],'pad3_%s'%POI[i],0.,0.25,1.,1.)
		pad4	= TPad('pad4_%s'%POI[i],'pad4_%s'%POI[i],0.,0.02,1.,0.25)
		p2pads	= [pad1,pad2,pad3,pad4]
                can.append(c)
		can2.append(c2)
                plots.append(p)
		plots2.append(p2)
		pads.append(p2pads)

	for i in range(3):

		can[i].cd()
		pads[i][0].Draw()
		pads[i][1].Draw()
		pads[i][0].SetLeftMargin(0.1)
		pads[i][1].SetLeftMargin(0.1)
		
		norm = wtmp.function('normfactor_3d_4fit_%s'%channel)

		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		wtmp.data('SMdatahist4fit').plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))
		normvalSM	= norm.getVal() * wtmp.data('SMdatahist4fit').sumEntries()

		wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots[i],RooFit.LineColor(kBlack),RooFit.Normalization(normvalSM, RooAbsReal.NumEvent))

		wtmp.data('neg_datahist4fit_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))
		wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		normvalneg = norm.getVal() * wtmp.data('SMdatahist4fit').sumEntries()
		wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots[i],RooFit.LineColor(kBlue),RooFit.Normalization(normvalneg, RooAbsReal.NumEvent))
	
		pullhist = plots[i].pullHist('h_neg_datahist4fit_%s'%POI[i],'aTGC_model_%s_Norm[rrv_mass_lvj]'%channel)
		
		if ch == 'el':
			plotmin = 1e-4
			#plotmin = 0.5
			plotmax = 50
			if signal_category == 'WZ':
				plotmin = 3e-4
		elif ch == 'mu':
			plotmin = 1e-2
			plotmax = 1e2
			if signal_category == 'WZ':
				plotmin = 1e-3
				plotmax = 50
		if signal_category=='WV':
			plotmin	= 1e-2
			plotmax	= 1.5e2
		plots[i].GetYaxis().SetRangeUser(plotmin,plotmax)
		pads[i][0].cd()
		pads[i][0].SetLogy()
		plots[i].SetTitle('')
		plots[i].GetYaxis().SetTitle('arbitrary units')
		plots[i].Draw()
		
		pads[i][1].cd()
		ratio_style = TH1D('ratio_style','ratio_style',(binhi-binlo)/100,binlo,binhi)
		ratio_style.SetMarkerStyle(21)
		ratio_style.SetMaximum(3)
		ratio_style.SetMinimum(-3)
		ratio_style.GetYaxis().SetNdivisions(7)
		ratio_style.GetYaxis().SetTitle('#frac{MC-Fit}{error}')
		ratio_style.GetYaxis().SetLabelSize(0.125)
		ratio_style.GetYaxis().SetTitleSize(0.2)
		ratio_style.GetYaxis().SetTitleOffset(0.2)
		ratio_style.Draw("")
		pullhist.SetLineColor(kBlue)
		pullhist.Draw("SAME E1")


		can[i].Update()
		if options.savep:
			can[i].SaveAs('docuplots/%s_neg_%s.pdf'%(POI[i],channel))
			can[i].SaveAs('docuplots/%s_neg_%s.png'%(POI[i],channel))
		

		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		wtmp.data('SMdatahist4fit').plotOn(plots2[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))
		wtmp.data('pos_datahist4fit_%s'%POI[i]).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))

		wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],RooFit.LineColor(kBlack),RooFit.Normalization(normvalSM, RooAbsReal.NumEvent))
		wtmp.var(POI[i]).setVal(par_max[POI[i]])
		normvalpos = norm.getVal() * wtmp.data('SMdatahist4fit').sumEntries()

        	wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],RooFit.LineColor(kBlue),RooFit.Normalization(normvalpos, RooAbsReal.NumEvent))
		
		wtmp.data('pos_datahist4fit_%s'%POI[i]).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		plots2[i].GetYaxis().SetRangeUser(plotmin,plotmax)
		plots2[i].GetYaxis().SetTitle('arbitrary units')
		can2[i].cd()
		pads[i][2].Draw()
		pads[i][3].Draw()
		pads[i][2].SetLeftMargin(0.1)
		pads[i][3].SetLeftMargin(0.1)
		plots2[i].SetTitle('')
		pads[i][2].SetLogy()
		pads[i][2].cd()
		plots2[i].Draw()
		pullhist2 = plots2[i].pullHist('h_pos_datahist4fit_%s'%POI[i],'aTGC_model_%s_Norm[rrv_mass_lvj]'%channel)
		pads[i][3].cd()
		ratio_style.Draw("")
		pullhist2.SetLineColor(kBlue)
		pullhist2.Draw("E1")

		can2[i].Update()
		if options.savep:
			can2[i].SaveAs('docuplots/%s_pos_%s.pdf'%(POI[i],channel))
			can2[i].SaveAs('docuplots/%s_pos_%s.png'%(POI[i],channel))
		
	if not options.batch:
		raw_input('plots plotted')
	for i in range(3):
		can[i].Delete()
		can2[i].Delete()
		for j in range(4):
			pads[i][j].Delete()


def make_input(ch = 'el', signal_cat = 'WW', binlo = 900, binhi = 3500):
	if signal_cat=='WV' and not options.oldfunc:
		raw_input("error functions not optimized for WV!")
	if signal_cat == 'WW':
		sigreg 	= 'lo'
		mj_lo	= 65.
		mj_hi	= 85.
	elif signal_cat == 'WZ':
		sigreg	= 'hi'
		mj_lo	= 85.
		mj_hi	= 105
	elif signal_cat == 'WV':
		sigreg	= 'all'
		mj_lo	= 65.
		mj_hi	= 105.
	
	if signal_cat=='WV':
		binlo4fit	= 900
	else:
		binlo4fit	= 600

	if not options.linter and not options.std:
		raise RuntimeError('no options chosen!')
	if binlo4fit > binlo:
		raise RuntimeError('cuts dont work for binlo < %s !'%binlo4fit)

	cat		= signal_cat
	channel		= cat+"_"+ch
	nbins		= (binhi-binlo)/100
	nbins4fit	= (binhi-binlo4fit)/100
	WS		= RooWorkspace("WS")

	#read data
	fileInData	= TFile.Open('Input/treeEDBR_data_xww_%s.root'%(ch))
	tree_data	= fileInData.Get('tree')  


	#read or make ATGC and new tree
	if options.newtrees:
		make_ATGCtree(ch, cat)
	if cat == 'WW':
		fileInATGC	= TFile.Open('Output/ATGC-Tree_lo_%s.root'%ch)
	elif cat == 'WZ':
		fileInATGC	= TFile.Open('Output/ATGC-Tree_hi_%s.root'%ch)
	elif cat == 'WV':
		fileInATGC	= TFile.Open('Output/ATGC-Tree_all_%s.root'%ch)
	else:
		raise RuntimeError('no such category: %s'%cat)
	treeATGC	= fileInATGC.Get('BasicTree')

	#prepare variables, parameters and temporary workspace
	wtmp		= RooWorkspace('wtmp')
	fitresults	= []
	##nuisance parameter to change all slope parameters by certain percentage (bigger for cb in WZ-cateogry)
	eps		= RooRealVar('slope_nuis','slope_nuis',1,0,2)
	eps.setConstant(kTRUE)
	eps4cbWZ	= RooFormulaVar('rel_slope_nuis4cbWZ','rel_slope_nuis4cbWZ','1+3*(@0-1)',RooArgList(eps))
	a1_4fit		= RooRealVar('a_SM_4fit_%s'%channel,'a_SM_4fit_%s'%channel,-0.1,-2,0)
	a1		= RooFormulaVar('a_SM_%s'%channel,'a_SM_%s'%channel,'@0*@1',RooArgList(a1_4fit,eps))

	cwww		= RooRealVar('cwww','cwww',0,-120,120);
	ccw		= RooRealVar('ccw','ccw',0,-200,200);
	cb		= RooRealVar('cb','cb',0,-600,600);
	cwww.setConstant(kTRUE);
	ccw.setConstant(kTRUE);
	cb.setConstant(kTRUE);

	##read workspace containing background pdfs
	fileInWs	= TFile.Open('Input/wwlvj_%s_HP%s_workspace.root'%(ch[:2],cat[1]))
	w		= fileInWs.Get('workspace4limit_')
	rrv_mass_lvj	= w.var('rrv_mass_lvj')
	rrv_mass_lvj.SetTitle('M_{WV}')
	rrv_mass_lvj.setRange(binlo,binhi)


	#create histograms for signal yields and save them to a file -> has to be done only once per MWW cut range
	if options.newhists:
		print '##########making histograms for aTGC working points##########'
		hists4scale	= []
		for para in POI:
			pos_hist	= TH1F('c_pos_hist4fit_%s'%para,'c_pos_hist4fit_%s'%para,nbins4fit,binlo4fit,binhi)
			pos_hist.Sumw2(kTRUE)
			pos_histhi	= TH1F('c_pos_hist_%s'%para,'c_pos_hist_%s'%para,nbins,binlo,binhi)
			pos_histhi.Sumw2(kTRUE)
			neg_hist	= TH1F('c_neg_hist4fit_%s'%para,'c_neg_hist4fit_%s'%para,nbins4fit,binlo4fit,binhi)
			neg_hist.Sumw2(kTRUE)
			neg_histhi	= TH1F('c_neg_hist_%s'%para,'c_neg_hist_%s'%para,nbins,binlo,binhi)
			neg_histhi.Sumw2(kTRUE)
			dif_hist	= TH1F('c_dif_hist4fit_%s'%para,'c_dif_hist4fit_%s'%para,nbins4fit,binlo4fit,binhi)
			dif_hist.Sumw2(kTRUE)
			hists4scale.append(pos_hist)
			hists4scale.append(neg_hist)
			hists4scale.append(dif_hist)
			hists4scale.append(pos_histhi)
			hists4scale.append(neg_histhi)
		SM_hist4fit	= TH1F('c_SM_hist4fit','c_SM_hist4fit',nbins4fit,binlo4fit,binhi)
		SM_hist4fit.Sumw2(kTRUE)
		SM_hist		= TH1F('c_SM_hist','c_SM_hist',nbins,binlo,binhi)
		SM_hist.Sumw2(kTRUE)
		histall3	= TH1F('c_histall3','c_histall3',nbins4fit,binlo4fit,binhi)
		histall3.Sumw2(kTRUE)
		dif_hist_ccw	= TH1F('c_dif_hist_ccw','dif_hist_ccw',nbins,binlo,binhi)
		dif_hist_ccw.Sumw2(kTRUE)
		dif_hist_cb	= TH1F('c_dif_hist_cb','dif_hist_cb',nbins,binlo,binhi)
		dif_hist_cb.Sumw2(kTRUE)
		fileOutHist	= TFile.Open('Output/hists4scale_%s_%s_%s.root'%(channel,binlo4fit,binhi),'recreate')
		for i in range(treeATGC.GetEntries()):
			if i%25000==0:
				print str(i) + '/' + str(treeATGC.GetEntries())
			treeATGC.GetEntry(i)
			if treeATGC.MWW > binlo4fit and treeATGC.MWW < binhi:
				if treeATGC.c_wwwl == 12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0 :
					hists4scale[0].Fill(treeATGC.MWW,treeATGC.totEventWeight)
					if treeATGC.MWW > binlo :
						hists4scale[3].Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if treeATGC.c_wwwl == -12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0 :
					hists4scale[1].Fill(treeATGC.MWW,treeATGC.totEventWeight)
					if treeATGC.MWW > binlo :
						hists4scale[4].Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if treeATGC.c_wwwl == 0 and treeATGC.c_wl == 20 and treeATGC.c_bl == 0 :
					hists4scale[5].Fill(treeATGC.MWW,treeATGC.totEventWeight)
					if treeATGC.MWW > binlo :
						hists4scale[8].Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if treeATGC.c_wwwl == 0 and treeATGC.c_wl == -20 and treeATGC.c_bl == 0 :
					hists4scale[6].Fill(treeATGC.MWW,treeATGC.totEventWeight)
					if treeATGC.MWW > binlo :
						hists4scale[9].Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == 60 :
					hists4scale[10].Fill(treeATGC.MWW,treeATGC.totEventWeight)
					if treeATGC.MWW > binlo :
						hists4scale[13].Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == -60 :
					hists4scale[11].Fill(treeATGC.MWW,treeATGC.totEventWeight)
					if treeATGC.MWW > binlo :
						hists4scale[14].Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0 :
					SM_hist4fit.Fill(treeATGC.MWW,treeATGC.totEventWeight)
					if treeATGC.MWW > binlo :
						SM_hist.Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if treeATGC.c_wwwl == -12 and treeATGC.c_wl == -20 and treeATGC.c_bl ==-60:
					histall3.Fill(treeATGC.MWW,treeATGC.totEventWeight)
		for i in range(nbins4fit):
			hists4scale[2].SetBinContent(i,hists4scale[0].GetBinContent(i)-hists4scale[1].GetBinContent(i))
			hists4scale[7].SetBinContent(i,hists4scale[5].GetBinContent(i)-hists4scale[6].GetBinContent(i))
			hists4scale[12].SetBinContent(i,hists4scale[10].GetBinContent(i)-hists4scale[11].GetBinContent(i))
		for i in range(nbins):
			dif_hist_ccw.SetBinContent(i,hists4scale[8].GetBinContent(i)-hists4scale[9].GetBinContent(i))
			dif_hist_cb.SetBinContent(i,hists4scale[13].GetBinContent(i)-hists4scale[14].GetBinContent(i))
		for i in range(len(hists4scale)):
			hists4scale[i].Write()
			hists4scale[i].Print()
		SM_hist4fit.Write()
		SM_hist.Write()
		histall3.Write()
		dif_hist_ccw.Write()
		dif_hist_cb.Write()
		print '-------> Written to file Output/%s'%fileOutHist.GetName()
		fileOutHist.Close()


	#get SM histogram and make RooDataHist
	fileInHist	= TFile.Open('Output/hists4scale_%s_%s_%s.root'%(channel,binlo4fit,binhi))
	rrv_mass_lvj.setRange(binlo4fit,binhi)
	SMdatahist4fit	= RooDataHist('SMdatahist4fit','SMdatahist4fit',RooArgList(rrv_mass_lvj),fileInHist.Get('c_SM_hist4fit'))
	rrv_mass_lvj.setRange(binlo,binhi)
	SMdatahist	= RooDataHist('SMdatahist','SMdatahist',RooArgList(rrv_mass_lvj),fileInHist.Get('c_SM_hist'))
	rrv_mass_lvj.setRange(binlo4fit,binhi)
	fileInHist.Close()

	#make SM pdf
	SMPdf		= RooExponential('SMPdf_%s'%channel,'SMPdf_%s'%channel,rrv_mass_lvj,a1)
	##actual fit to determine SM shape parameter
	fitresSM	= SMPdf.fitTo(SMdatahist4fit, RooFit.SumW2Error(kTRUE), RooFit.Save(kTRUE))
	fitresults.append(fitresSM)
	a1_4fit.setConstant(kTRUE)
	N_SM		= RooRealVar('N_SM_%s'%channel,'N_SM_%s'%channel,SMdatahist4fit.sumEntries())
	N_SM.setConstant(kTRUE)

	getattr(wtmp,'import')(cwww);
	getattr(wtmp,'import')(ccw);
	getattr(wtmp,'import')(cb);
	getattr(wtmp,'import')(eps4cbWZ)
	getattr(wtmp,'import')(SMdatahist4fit)
	getattr(wtmp,'import')(SMdatahist)
	getattr(wtmp,'import')(N_SM)

	#define parameter ranges for each channel

	if cat=='WW':
		Erf_width_cwww	= RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,500.,100.,1000)
		if ch=='el':
			Erf_width_ccw	= RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,500.,100.,2500.)
			Erf_width_cb	= RooRealVar('Erf_width_cb_%s'%channel,'Erf_width_cb_%s'%channel,500.,100.,5000.)
		elif ch=='mu':
			Erf_width_ccw	= RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,300.,10.,5000.)
			Erf_width_cb	= RooRealVar('Erf_width_cb_%s'%channel,'Erf_width_cb_%s'%channel,300.,10.,5000.)
		Erf_offset_cwww	= RooRealVar('Erf_offset_cwww_%s'%channel,'Erf_offset_cwww_%s'%channel,1000.,600.,1500.)
		Erf_offset_ccw	= RooRealVar('Erf_offset_ccw_%s'%channel,'Erf_offset_ccw_%s'%channel,1000.,600.,1500.)
		Erf_offset_cb	= RooRealVar('Erf_offset_cb_%s'%channel,'Erf_offset_cb_%s'%channel,1000.,600.,1500.)
	elif cat=='WZ' or (ch=='mu' and cat=='WV'):
		if ch=='el':
			Erf_width_cwww	= RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,1500.,375.,2000.)
			Erf_width_ccw	= RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,1500.,375.,2000.)
		elif ch=='mu':
			Erf_width_cwww	= RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,500.,375.,600.)
			Erf_width_ccw	= RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,500.,375.,600.)
		Erf_width_cb	= RooRealVar('Erf_width_cb_%s'%channel,'Erf_width_cb_%s'%channel,0.,0.,1.)#not used in WZ
		Erf_offset_cwww	= RooRealVar('Erf_offset_cwww_%s'%channel,'Erf_offset_cwww_%s'%channel,1000.,500.,1500.)
		Erf_offset_ccw	= RooRealVar('Erf_offset_ccw_%s'%channel,'Erf_offset_ccw_%s'%channel,1000.,500.,1500.)
		Erf_offset_cb	= RooRealVar('Erf_offset_cb_%s'%channel,'Erf_offset_cb_%s'%channel,1000.,500.,1500.)
	elif cat=='WV':
		Erf_width_cwww	= RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,900.,-2000.,2000)
		if ch=='el':
			Erf_width_ccw	= RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,900.,-2000.,2500.)
			Erf_width_cb	= RooRealVar('Erf_width_cb_%s'%channel,'Erf_width_cb_%s'%channel,900.,-2000.,2500.)
		elif ch=='mu':
			Erf_width_ccw	= RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,300.,10.,5000.)
			Erf_width_cb	= RooRealVar('Erf_width_cb_%s'%channel,'Erf_width_cb_%s'%channel,300.,10.,5000.)
		Erf_offset_cwww	= RooRealVar('Erf_offset_cwww_%s'%channel,'Erf_offset_cwww_%s'%channel,900.,800.,1500.)
		Erf_offset_ccw	= RooRealVar('Erf_offset_ccw_%s'%channel,'Erf_offset_ccw_%s'%channel,900.,800.,1500.)
		Erf_offset_cb	= RooRealVar('Erf_offset_cb_%s'%channel,'Erf_offset_cb_%s'%channel,900.,800.,1500.)
	Erf_offset_cwww.setConstant(kTRUE)
	Erf_width_cwww.setConstant(kTRUE)
	Erf_offset_ccw.setConstant(kTRUE)
	Erf_width_ccw.setConstant(kTRUE)
	Erf_offset_cb.setConstant(kTRUE)
	Erf_width_cb.setConstant(kTRUE)
	getattr(wtmp,'import')(Erf_width_cwww)
	getattr(wtmp,'import')(Erf_offset_cwww)
	getattr(wtmp,'import')(Erf_width_ccw)
	getattr(wtmp,'import')(Erf_offset_ccw)
	getattr(wtmp,'import')(Erf_width_cb)
	getattr(wtmp,'import')(Erf_offset_cb)




	for i in range(len(POI)):	
		fileInHist	= TFile.Open('Output/hists4scale_%s_%s_%s.root'%(channel,binlo4fit,binhi))
		rrv_mass_lvj.setRange(binlo4fit,binhi)
		pos_datahist4fit	= RooDataHist('pos_datahist4fit_%s'%POI[i],'pos_datahist4fit_%s'%POI[i],RooArgList(rrv_mass_lvj),fileInHist.Get('c_pos_hist4fit_%s'%POI[i]))
		neg_datahist4fit	= RooDataHist('neg_datahist4fit_%s'%POI[i],'neg_datahist4fit_%s'%POI[i],RooArgList(rrv_mass_lvj),fileInHist.Get('c_neg_hist4fit_%s'%POI[i]))
		rrv_mass_lvj.setRange(binlo,binhi)		
		pos_datahist		= RooDataHist('pos_datahist_%s'%POI[i],'pos_datahist_%s'%POI[i],RooArgList(rrv_mass_lvj),fileInHist.Get('c_pos_hist_%s'%POI[i]))
		neg_datahist		= RooDataHist('neg_datahist_%s'%POI[i],'neg_datahist_%s'%POI[i],RooArgList(rrv_mass_lvj),fileInHist.Get('c_neg_hist_%s'%POI[i]))
		if cat=='WV' and POI[i]!='cwww':
			dif_datahist		= RooDataHist('dif_datahist_%s'%POI[i],'dif_datahist_%s'%POI[i],RooArgList(rrv_mass_lvj),fileInHist.Get('c_dif_hist_%s'%POI[i]))
			getattr(wtmp,'import')(dif_datahist)
		rrv_mass_lvj.setRange(binlo4fit,binhi)
		dif_datahist4fit	= RooDataHist('dif_datahist4fit_%s'%POI[i],'dif_datahist4fit_%s'%POI[i],RooArgList(rrv_mass_lvj),fileInHist.Get('c_dif_hist4fit_%s'%POI[i]))
		fileInHist.Close()


		getattr(wtmp,'import')(pos_datahist4fit)
		getattr(wtmp,'import')(neg_datahist4fit)
		getattr(wtmp,'import')(pos_datahist)
		getattr(wtmp,'import')(neg_datahist)
		getattr(wtmp,'import')(dif_datahist4fit)
		getattr(WS,'import')(pos_datahist4fit)
		getattr(WS,'import')(neg_datahist4fit)
		getattr(WS,'import')(pos_datahist)
		getattr(WS,'import')(neg_datahist)
		getattr(WS,'import')(dif_datahist4fit)
#		#get scaling parabel from yields in actual cut range 
		hist4scale = TH1F('hist4scale_%s'%POI[i],'hist4scale_%s'%POI[i],3,-1.5*par_max[POI[i]],1.5*par_max[POI[i]])
		hist4scale.SetBinContent(1,neg_datahist.sumEntries()/SMdatahist.sumEntries())
		hist4scale.SetBinContent(2,1)
		hist4scale.SetBinContent(3,pos_datahist.sumEntries()/SMdatahist.sumEntries())
		#fit parabel
		gROOT.SetBatch(True)
		hist4scale.Fit('pol2')
		if options.yieldplots:
			gROOT.SetBatch(False)
			cc1 = TCanvas()
			cc1.cd()
			hist4scale.GetXaxis().SetTitle(par_titles[POI[i]]+' (TeV^{-2})')
			hist4scale.GetYaxis().SetTitle('N_{events}^{SM+%s} / N_{events}^{SM}'%par_titles[POI[i]])
			hist4scale.GetYaxis().SetTitleSize(0.04)
			hist4scale.Draw()
			cc1.Update()
			#if not options.batch:
			#	raw_input(POI[i])
			cc1.SaveAs('docuplots/yields_%s_%s.pdf'%(POI[i],channel))
			cc1.SaveAs('docuplots/yields_%s_%s.png'%(POI[i],channel))
			cc1.Close()
		if not options.batch:
			gROOT.SetBatch(False)
		fitfunc		= hist4scale.GetFunction('pol2')
		par0		= RooRealVar('par0_%s_%s'%(POI[i],channel),'par0_%s_%s'%(POI[i],channel),fitfunc.GetParameter(0)); 		par0.setConstant(kTRUE);
		par1		= RooRealVar('par1_%s_%s'%(POI[i],channel),'par1_%s_%s'%(POI[i],channel),fitfunc.GetParameter(1)); 		par1.setConstant(kTRUE);
		par2		= RooRealVar('par2_%s_%s'%(POI[i],channel),'par2_%s_%s'%(POI[i],channel),fitfunc.GetParameter(2)); 		par2.setConstant(kTRUE);
		if options.noSMint:
			par1.setVal(0)
		#again for range used for fit, only needed for plots
		hist4fit4plot	=TH1F('hist4fit4plot_%s'%POI[i],'hist4fit4plot_%s'%POI[i],3,-1.5*par_max[POI[i]],1.5*par_max[POI[i]])
		hist4fit4plot.SetBinContent(1,neg_datahist4fit.sumEntries()/SMdatahist4fit.sumEntries())
		hist4fit4plot.SetBinContent(2,1)
		hist4fit4plot.SetBinContent(3,pos_datahist4fit.sumEntries()/SMdatahist4fit.sumEntries())
		gROOT.SetBatch(True)
		hist4fit4plot.Fit('pol2')
		if not options.batch:
			gROOT.SetBatch(False)
		if options.yieldplots:
			cc2		= TCanvas()
			cc2.cd()
			hist4fit4plot.Draw()
			cc2.Draw()
			cc2.Update()
			#if not options.batch:
			#	raw_input(POI[i])
			cc2.Close()
		fitfunc4fit	= hist4fit4plot.GetFunction('pol2')
		par0_4fit	= RooRealVar('par0_4fit_%s_%s'%(POI[i],channel),'par0_4fit_%s_%s'%(POI[i],channel),fitfunc4fit.GetParameter(0)); 		par0_4fit.setConstant(kTRUE);
		par1_4fit	= RooRealVar('par1_4fit_%s_%s'%(POI[i],channel),'par1_4fit_%s_%s'%(POI[i],channel),fitfunc4fit.GetParameter(1)); 		par1_4fit.setConstant(kTRUE);
		par2_4fit	= RooRealVar('par2_4fit_%s_%s'%(POI[i],channel),'par2_4fit_%s_%s'%(POI[i],channel),fitfunc4fit.GetParameter(2)); 		par2_4fit.setConstant(kTRUE);
#i
		N_pos_tmp 	= pos_datahist4fit.sumEntries()
		N_neg_tmp	= neg_datahist4fit.sumEntries()
		if options.std:
			N_quad		= RooRealVar('N_quad_%s_%s'%(POI[i],channel),'N_quad_%s_%s'%(POI[i],channel), N_pos_tmp-N_SM.getVal() )
		else:
			N_quad		= RooRealVar('N_quad_%s_%s'%(POI[i],channel),'N_quad_%s_%s'%(POI[i],channel), ((N_pos_tmp+N_neg_tmp)/2)-N_SM.getVal() )
		N_lin		= RooRealVar('N_lin_%s_%s'%(POI[i],channel),'N_lin_%s_%s'%(POI[i],channel), (N_pos_tmp-N_neg_tmp)/2 )

		#scaleshape is the relative change to SM
		scaleshape	= RooFormulaVar('scaleshape_%s_%s'%(POI[i],channel),'scaleshape_%s_%s'%(POI[i],channel),\
					'(@0+@1*@3+@2*@3**2)-1',\
					RooArgList(par0,par1,par2,wtmp.var(POI[i])))
		scaleshape4fit = RooFormulaVar('scaleshape4fit_%s_%s'%(POI[i],channel),'scaleshape4fit_%s_%s'%(POI[i],channel),\
					'(@0+@1*@3+@2*@3**2)-1',\
					RooArgList(par0_4fit,par1_4fit,par2_4fit,wtmp.var(POI[i])))			
		if cat=='WWWV':
			a2_4fit		= RooRealVar('a_quad_4fit_%s_%s'%(POI[i],channel),'a_quad_4fit_%s_%s'%(POI[i],channel),-0.0001,-0.1,0.)
			a3_4fit		= RooRealVar('a_lin_4fit_%s_%s'%(POI[i],channel),'a_lin_4fit_%s_%s'%(POI[i],channel),-0.0001,-0.1,0.)
		else:
			a2_4fit		= RooRealVar('a_quad_4fit_%s_%s'%(POI[i],channel),'a_quad_4fit_%s_%s'%(POI[i],channel),-0.001,-0.01,0.)
			a3_4fit		= RooRealVar('a_lin_4fit_%s_%s'%(POI[i],channel),'a_lin_4fit_%s_%s'%(POI[i],channel),-0.001,-0.01,0.)
		a2_4fit.setConstant(kTRUE)
		a3_4fit.setConstant(kTRUE)
		##bigger uncertainty for cb in WZ-category
		if (cat=='WZ' and POI[i]=='cb'):
			a2		= RooFormulaVar('a_quad_nuis_%s_%s'%(POI[i],channel),'a_quad_nuis_%s_%s'%(POI[i],channel),'@0*@1',RooArgList(a2_4fit,eps4cbWZ))
			a3		= RooFormulaVar('a_lin_nuis_%s_%s'%(POI[i],channel),'a_lin_nuis_%s_%s'%(POI[i],channel),'@0*@1',RooArgList(a3_4fit,eps4cbWZ))
			cPdf_quad	= RooExponential('Pdf_quad_%s_%s'%(POI[i],channel),'Pdf_quad_%s_%s'%(POI[i],channel),rrv_mass_lvj,a2)
		else:
			a2		= RooFormulaVar('a_quad_nuis_%s_%s'%(POI[i],channel),'a_quad_nuis_%s_%s'%(POI[i],channel),'@0*@1',RooArgList(a2_4fit,eps))
			a3		= RooFormulaVar('a_lin_nuis_%s_%s'%(POI[i],channel),'a_lin_nuis_%s_%s'%(POI[i],channel),'@0*@1',RooArgList(a3_4fit,eps))
			if options.oldfunc:
				cPdf_quad	= RooExponential('Pdf_quad_%s_%s'%(POI[i],channel),'Pdf_quad_%s_%s'%(POI[i],channel),rrv_mass_lvj,a2)			
			else:
				cPdf_quad	= RooErfExpPdf('Pdf_quad_%s_%s'%(POI[i],channel),'Pdf_quad_%s_%s'%(POI[i],channel),rrv_mass_lvj,a2,wtmp.var('Erf_offset_%s_%s'%(POI[i],channel)),wtmp.var('Erf_width_%s_%s'%(POI[i],channel)))
		cPdf_lin	= RooExponential('Pdf_lin_%s_%s'%(POI[i],channel),'Pdf_lin_%s_%s'%(POI[i],channel),rrv_mass_lvj,a3)

		getattr(wtmp,'import')(cPdf_quad,RooFit.RecycleConflictNodes())
		getattr(wtmp,'import')(cPdf_lin,RooFit.RecycleConflictNodes())
		getattr(wtmp,'import')(N_quad)
		getattr(wtmp,'import')(N_lin)
		getattr(wtmp,'import')(scaleshape)
		getattr(wtmp,'import')(scaleshape4fit)

		wtmp.Print()

	fileInHist.Close()

	#make model
	paralist	= RooArgList(N_SM)


	#include SM- and aTGC-interference
	if options.linter:
		#get ratio WW/WZ in WW- and WZ-category for each aTGC-parameter
		if options.make_ratios:
			WWinWW,WZinWW,WWinWZ,WZinWZ = 0,0,0,0
			print 'reading Output/ATGC-Tree_WW_lo_%s.root for WW/WZ ratios'%(ch)
			fileInWWinWW	= TFile.Open('Output/ATGC-Tree_WW_lo_%s.root'%(ch))
			treeInWWinWW	= fileInWWinWW.Get('BasicTree')
			treeInWWinWW.SetBranchStatus('*',0)
			treeInWWinWW.SetBranchStatus('totEventWeight',1)
			treeInWWinWW.SetBranchStatus('MWW',1)
			for i in range(treeInWWinWW.GetEntries()):
				treeInWWinWW.GetEntry(i)
				if treeInWWinWW.MWW > binlo4fit:
					WWinWW		+= treeInWWinWW.totEventWeight
			print 'reading Output/ATGC-Tree_WW_hi_%s.root for WW/WZ ratios'%(ch)
			fileInWWinWZ	= TFile.Open('Output/ATGC-Tree_WW_hi_%s.root'%(ch))
			treeInWWinWZ	= fileInWWinWZ.Get('BasicTree')
			treeInWWinWZ.SetBranchStatus('*',0)
			treeInWWinWZ.SetBranchStatus('totEventWeight',1)
			treeInWWinWZ.SetBranchStatus('MWW',1)
			for i in range(treeInWWinWZ.GetEntries()):
				treeInWWinWZ.GetEntry(i)
				if treeInWWinWZ.MWW > binlo4fit:
					WWinWZ		+= treeInWWinWZ.totEventWeight
			print 'reading Output/ATGC-Tree_WZ_lo_%s.root for WW/WZ ratios'%(ch)
			fileInWZinWW	= TFile.Open('Output/ATGC-Tree_WZ_lo_%s.root'%(ch))
			treeInWZinWW	= fileInWZinWW.Get('BasicTree')
			treeInWZinWW.SetBranchStatus('*',0)
			treeInWZinWW.SetBranchStatus('totEventWeight',1)
			treeInWZinWW.SetBranchStatus('MWW',1)
			for i in range(treeInWZinWW.GetEntries()):
				treeInWZinWW.GetEntry(i)
				if treeInWZinWW.MWW > binlo4fit:
					WZinWW		+= treeInWZinWW.totEventWeight
			print 'reading Output/ATGC-Tree_WZ_hi_%s.root for WW/WZ ratios'%(ch)
			fileInWZinWZ	= TFile.Open('Output/ATGC-Tree_WZ_hi_%s.root'%(ch))
			treeInWZinWZ	= fileInWZinWZ.Get('BasicTree')
			treeInWZinWZ.SetBranchStatus('*',0)
			treeInWZinWZ.SetBranchStatus('totEventWeight',1)
			treeInWZinWZ.SetBranchStatus('MWW',1)
			for i in range(treeInWZinWZ.GetEntries()):
				treeInWZinWZ.GetEntry(i)
				if treeInWZinWZ.MWW > binlo4fit:
					WZinWZ		+= treeInWZinWZ.totEventWeight
			ratio_WW_reg	= (WWinWW)/(WWinWW+WZinWW)
			ratio_WZ_reg	= (WWinWZ)/(WWinWZ+WZinWZ)
			WWWZ_ratios	= {'WW' : ratio_WW_reg, 'WZ' : ratio_WZ_reg}
		else:
			#pre-approval values:
			#if ch=='el':
			#	WWWZ_ratios 	= {'WW' : 0.85870, 'WZ' : 0.47618}
			#elif ch=='mu':
			#	WWWZ_ratios	= {'WW' : 0.86257, 'WZ' : 0.51776}
			if ch=='el':
				WWWZ_ratios 	= {'WW' : 0.8046, 'WZ' : 0.5091, 'WV' : 0.5}
			elif ch=='mu':
				WWWZ_ratios	= {'WW' : 0.8615, 'WZ' : 0.5159, 'WV' : 0.5}

		if options.make_ratios:
			print WWWZ_ratios
			raw_input('..')
		WWratio		= WWWZ_ratios[cat]
		WZratio		= 1-WWWZ_ratios[cat]

		#read results of generator level study
		fileInterWW	= TFile.Open('Input/genlevel_WW_%s.root'%ch)
		w2WW		= fileInterWW.Get('w2')
		fileInterWZ	= TFile.Open('Input/genlevel_WZ_%s.root'%ch)
		w2WZ		= fileInterWZ.Get('w2')
		#make weighted mean with WW/WZ ratio	
		a5WW		= w2WW.var('a5').getVal()
		a5WZ		= w2WZ.var('a5').getVal()
		a5val		= WWratio*a5WW + WZratio*a5WZ
		a5_tmp		= RooRealVar('a_cwww_ccw_%s'%channel,'a_cwww_ccw_%s'%channel,a5val)
		a7WW		= w2WW.var('a7').getVal()
		a7WZ		= w2WZ.var('a7').getVal()
		a7val		= WWratio*a7WW + WZratio*a7WZ
		a7_tmp		= RooRealVar('a_ccw_cb_%s'%channel,'a_ccw_cb_%s'%channel,a7val)
		a5_tmp.setConstant(kTRUE)
		a7_tmp.setConstant(kTRUE)
		#bigger uncertainty for c_B in WZ category
		if cat=='WZ':
			a5		= RooFormulaVar('a_cwww_ccw_nuis_%s'%channel,'a_cwww_ccw_nuis_%s'%channel,'@0*@1',RooArgList(a5_tmp,eps))
			a7		= RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,eps4cbWZ))
		else:
			a5		= RooFormulaVar('a_cwww_ccw_nuis_%s'%channel,'a_cwww_ccw_nuis_%s'%channel,'@0*@1',RooArgList(a5_tmp,eps))
			a7		= RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,eps))

		NSMWW	= w2WW.var('N_SM4fit').getVal()
		NSMWZ	= w2WZ.var('N_SM4fit').getVal()
		NSM		= WWratio*NSMWW + WZratio*NSMWZ

		Pdf_cwww_ccw	= RooExponential('Pdf_cwww_ccw_%s'%channel,'Pdf_cwww_ccw_%s'%channel,rrv_mass_lvj,a5)
		Pdf_ccw_cb	= RooExponential('Pdf_ccw_cb_%s'%channel,'Pdf_ccw_cb_%s'%channel,rrv_mass_lvj,a7)

		#get factor to scale number of events to simulation level (ratio of N_events for all atgc-parameters negative)	
		fileInHist	= TFile.Open('Output/hists4scale_%s_%s_%s.root'%(channel,binlo4fit,binhi))
		hist_all3	= fileInHist.Get('c_histall3')
		hist_all3.Print()
		hist_all3.Sumw2(kTRUE)
		datahist_all3	= RooDataHist('datahist_all3','datahist_all3',RooArgList(rrv_mass_lvj),hist_all3)
		getattr(wtmp,'import')(datahist_all3)
		fileInHist.Close()
		N_4norm		= WWratio*w2WW.var('N_4norm4fit').getVal() + WZratio*w2WZ.var('N_4norm4fit').getVal()
		if options.noatgcint:
			cf	= 0
		else:
			cf		= datahist_all3.sumEntries() / N_4norm


		#define more coefficients, make weighted means
		N1220	= WWratio*w2WW.var('N_cwww_ccw__12__204fit').getVal() + WZratio*w2WZ.var('N_cwww_ccw__12__204fit').getVal()
		N2060	= WWratio*w2WW.var('N_ccw_cb__20__604fit').getVal() + WZratio*w2WZ.var('N_ccw_cb__20__604fit').getVal()
		N12	= WWratio*w2WW.var('N_cwww_124fit').getVal() + WZratio*w2WZ.var('N_cwww_124fit').getVal()
		N12_	= WWratio*w2WW.var('N_cwww__124fit').getVal() + WZratio*w2WZ.var('N_cwww__124fit').getVal()
		N20	= WWratio*w2WW.var('N_ccw_204fit').getVal() + WZratio*w2WZ.var('N_ccw_204fit').getVal()
		N20_	= WWratio*w2WW.var('N_ccw__204fit').getVal() + WZratio*w2WZ.var('N_ccw__204fit').getVal()
		N60	= WWratio*w2WW.var('N_cb_604fit').getVal() + WZratio*w2WZ.var('N_cb_604fit').getVal()
		N60_	= WWratio*w2WW.var('N_cb__604fit').getVal() + WZratio*w2WZ.var('N_cb__604fit').getVal()


		##no cwww-cb interference
		N_cwww_ccw	= RooRealVar('N_cwww_ccw_%s'%channel,'N_cwww_ccw_%s'%channel,\
						cf*((N1220+NSM)-(N12+N20)))
		N_ccw_cb	= RooRealVar('N_ccw_cb_%s'%channel,'N_ccw_cb_%s'%channel,\
						cf*((N2060+NSM)-(N20+N60)))

		paralist.add(RooArgList(wtmp.function('N_quad_%s_%s'%(POI[0],channel)),wtmp.var('cwww'),\
					wtmp.function('N_quad_%s_%s'%(POI[1],channel)),wtmp.function('N_lin_%s_%s'%(POI[1],channel)),wtmp.var('ccw'),\
					wtmp.function('N_quad_%s_%s'%(POI[2],channel)),wtmp.function('N_lin_%s_%s'%(POI[2],channel)),wtmp.var('cb')))
		paralist.add(RooArgList(N_cwww_ccw,N_ccw_cb))
		#neglect SM interference for cwww
		cwww_s		= '+@1*(@2/12)**2'
		ccw_s		= '+@3*(@5/20)**2+@4*(@5/20)'
		cb_s		= '+@6*(@8/60)**2+@7*(@8/60)'
		cwww_ccw_s	= '+@9*(@2/12)*(@5/20)'
		ccw_cb_s	= '+@10*(@5/20)*(@8/60)'
		Pdf_norm	= RooFormulaVar( 'Pdf_norm_%s'%channel,'Pdf_norm_%s'%channel,\
							'@0'+cwww_s+ccw_s+cb_s+cwww_ccw_s+ccw_cb_s,paralist)
		paralistN	= RooArgList()
		for i in range(11):
			paralistN.add(RooArgList(paralist.at(i)))
		paralistN.add(RooArgList(Pdf_norm))

		N1		= RooFormulaVar( 'N1_%s'%channel,'N1_%s'%channel,'@0/@11',paralistN )
		N2		= RooFormulaVar( 'N2_%s'%channel,'N2_%s'%channel,'(@1*(@2/12)**2)/@11',paralistN )
		#no SM-interference for c_WWW
		N4		= RooFormulaVar( 'N4_%s'%channel,'N4_%s'%channel,'(@3*(@5/20)**2)/@11',paralistN )
		N5		= RooFormulaVar( 'N5_%s'%channel,'N5_%s'%channel,'(@4*(@5/20))/@11',paralistN )
		N6		= RooFormulaVar( 'N6_%s'%channel,'N6_%s'%channel,'(@6*(@8/60)**2)/@11',paralistN )
		N7		= RooFormulaVar( 'N7_%s'%channel,'N7_%s'%channel,'(@7*(@8/60))/@11',paralistN )
		#no aTGC-interference for c_WWW/c_B
		N8		= RooFormulaVar( 'N8_%s'%channel,'N8_%s'%channel,'(@9*(@2/12)*(@5/20))/@11',paralistN )
		N10		= RooFormulaVar( 'N10_%s'%channel,'N10_%s'%channel,'(@10*(@5/20)*(@8/60))/@11',paralistN )

		N_list		= RooArgList(N1,N2,N4,N5,N6,N7)
		N_list.add(RooArgList(N8,N10))
		Pdf_list	= RooArgList(SMPdf)
		Pdf_list.add(RooArgList(wtmp.pdf('Pdf_quad_%s_%s'%(POI[0],channel)),\
					wtmp.pdf('Pdf_quad_%s_%s'%(POI[1],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[1],channel)),\
					wtmp.pdf('Pdf_quad_%s_%s'%(POI[2],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[2],channel))))
		Pdf_list.add(RooArgList(Pdf_cwww_ccw,Pdf_ccw_cb))
		model		= RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list)
		N_list.Print()
		Pdf_list.Print()
		model.Print()

		scale_list	= RooArgList(wtmp.function('scaleshape_%s_%s'%(POI[0],channel)),\
						wtmp.function('scaleshape_%s_%s'%(POI[1],channel)),\
						wtmp.function('scaleshape_%s_%s'%(POI[2],channel)))
		normfactor_3d	= RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list)
		scale_list4fit		= RooArgList(wtmp.function('scaleshape4fit_%s_%s'%(POI[0],channel)),\
						wtmp.function('scaleshape4fit_%s_%s'%(POI[1],channel)),\
						wtmp.function('scaleshape4fit_%s_%s'%(POI[2],channel)))
		normfactor_3d_4fit	= RooFormulaVar('normfactor_3d_4fit_%s'%channel,'normfactor_3d_4fit_%s'%channel,'1+@0+@1+@2',scale_list4fit)
		wtmp.Print()

		#fit 3 pdfs
		rrv_mass_lvj.setRange(binlo4fit,binhi)
		for i in range(3):
			for j in range(3):
				wtmp.var(POI[j]).setVal(0)
			wtmp.var(POI[i]).setVal(par_max[POI[i]])

			wtmp.var('a_lin_4fit_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			#fit SM-interference first
			if POI[i] != 'cwww':
				N_SM_tmp = N_SM.getVal()
				N_SM.setVal(0)
				N_quad_tmp = wtmp.var('N_quad_%s_%s'%(POI[i],channel)).getVal()
				wtmp.var('N_quad_%s_%s'%(POI[i],channel)).setVal(0)
				if cat!='WV':
					fitres1		= model.fitTo(wtmp.data('dif_datahist4fit_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
				else:
					fitres1		= model.fitTo(wtmp.data('dif_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))

				fitresults.append(fitres1)

				if options.linfitplot:
					c = TCanvas('c','c',1)
					c.cd();
					p = rrv_mass_lvj.frame()
					wtmp.data('dif_datahist4fit_%s'%POI[i]).plotOn(p,RooFit.DrawOption('E1'))
					wtmp.data('dif_datahist4fit_%s'%POI[i]).Print("V")
					model.plotOn(p,RooFit.LineWidth(1),RooFit.LineColor(kBlack))
					p.GetYaxis().SetRangeUser(0.001,100)
					p.Draw()
					c.SetLogy()
					c.Draw()
					c.Update()
					raw_input("zzz")
					c.Close()

				wtmp.var('a_lin_4fit_%s_%s'%(POI[i],channel)).setConstant(kTRUE)
				N_SM.setVal(N_SM_tmp)
				w_tmp = wtmp.var('N_quad_%s_%s'%(POI[i],channel)).setVal(N_quad_tmp)

			wtmp.var('a_quad_4fit_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			wtmp.var('Erf_offset_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			wtmp.var('Erf_width_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			model.Print()
			for j in range(8):
				print N_list.at(j).GetName() + ' : ' + str(N_list.at(j).getVal())
			for j in range(11):
				print paralist.at(i).GetName() + ' : ' + str(paralist.at(i).getVal())
			if cat=='WV':
				rrv_mass_lvj.setRange(900,3500)
				fitres2		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
			else:
				fitres2		= model.fitTo(wtmp.data('pos_datahist4fit_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
			fitresults.append(fitres2)
			wtmp.var('a_quad_4fit_%s_%s'%(POI[i],channel)).setConstant(kTRUE)
			wtmp.var('Erf_offset_%s_%s'%(POI[i],channel)).setConstant(kTRUE)
			wtmp.var('Erf_width_%s_%s'%(POI[i],channel)).setConstant(kTRUE)
			if POI[i]!='cwww':
				fitres1.Print()
			fitres2.Print()
			#raw_input("CC")
		for i in range(3):
			wtmp.var(POI[i]).setVal(0)
		
		if options.noSMint:
			wtmp.var('N_lin_ccw_%s'%channel).setVal(0)
			wtmp.var('N_lin_cb_%s'%channel).setVal(0)
		if options.atgc:
			#go from EFT parametrization to Lagrangian approach
		
			Z_mass		= 0.0911876
			W_mass		= 0.080385
			G_F		= 11.663787
			g_weak		= math.sqrt((8*G_F*W_mass**2)/(math.sqrt(2)))
			theta_W		= math.acos(W_mass/Z_mass)
			tan_theta_W 	= math.tan(theta_W)
			sin_theta_W	= math.sin(theta_W)
	
			coeff_cb1	= RooRealVar('coeff_cb1','coeff_cb1',2/(tan_theta_W*tan_theta_W*Z_mass*Z_mass))
			coeff_cb2	= RooRealVar('coeff_cb2','coeff_cb2',2/(sin_theta_W*sin_theta_W*W_mass*W_mass))
			coeff_ccw	= RooRealVar('coeff_ccw','coeff_ccw',2/Z_mass*Z_mass)
			coeff_cwww	= RooRealVar('coeff_cwww','coeff_cwww',2/(3*g_weak*g_weak*W_mass*W_mass))
	
			dg1z		= RooRealVar('dg1z','dg1z',0,-1,1)
			lZ		= RooRealVar('lZ','lZ',0,-1,1)
			dkz		= RooRealVar('dkz','dkz',0,-1,1)
			dg1z.setConstant(kTRUE)
			lZ.setConstant(kTRUE)
			dkz.setConstant(kTRUE)

			cwww		= RooFormulaVar('cwww_atgc','cwww_atgc','@0*@1',RooArgList(lZ,coeff_cwww))
			ccw		= RooFormulaVar('ccw_atgc','ccw_atgc','@0*@1',RooArgList(dg1z,coeff_ccw))
			cb		= RooFormulaVar('cb_atgc','cb_atgc','@0*@1-@2*@3',RooArgList(dg1z,coeff_cb1,dkz,coeff_cb2))
			atgc_pars	= RooArgList(cwww,ccw,cb)
			cwww.Print()
			ccw.Print()
			cb.Print()
			N_list_atgc	= RooArgList()
			scale_list_atgc = RooArgList()

			for i in range(8):
				customize_N	= RooCustomizer(N_list.at(i),'customize_N')
				for j in range(3):
					customize_N.replaceArg(wtmp.var(POI[j]),atgc_pars.at(j))
				N_list_atgc.add(customize_N.build())
				N_list_atgc.at(i).SetName(N_list.at(i).GetName())
			model_atgc	= RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list_atgc)
			for i in range(3):
				customize_scale	= RooCustomizer(scale_list.at(i),'customize_scale')
				for j in range(3):
					customize_scale.replaceArg(wtmp.var(POI[j]),atgc_pars.at(j))
				scale_list_atgc.add(customize_scale.build())
				
			normfactor_3d	= RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list_atgc)
		
			getattr(wtmp,'import')(model_atgc)
			getattr(WS,'import')(model_atgc)
			getattr(wtmp,'import')(normfactor_3d,RooFit.RecycleConflictNodes())
			getattr(WS,'import')(normfactor_3d,RooFit.RecycleConflictNodes())
			WS.Print()
			raw_input(channel)
		else:
			model.Print()
			getattr(WS,'import')(normfactor_3d)	
			getattr(wtmp,'import')(normfactor_3d)
			getattr(wtmp,'import')(normfactor_3d_4fit)
			getattr(wtmp,'import')(model)
			getattr(WS,'import')(model)

		
		#print coefficients to see contribution for all atgc-parameter positive
		for i in range(3):
			wtmp.var(POI[i]).setVal(par_max[POI[i]])
		for i in range(11):
			print paralist.at(i).GetName() + ' : ' + str(paralist.at(i).getVal())

		#print fitresults
		for i in range(8):
			print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())

		a5.Print()
		a7.Print()

	#no interference
	if options.std:
		raw_input('not sure if this works, probably not...')
		paralist.add(RooArgList(wtmp.function('N_quad_%s_%s'%(POI[0],channel)),wtmp.var('cwww'),\
					wtmp.function('N_quad_%s_%s'%(POI[1],channel)),wtmp.var('ccw'),\
					wtmp.function('N_quad_%s_%s'%(POI[2],channel)),wtmp.var('cb')))
		N1		= RooFormulaVar( 'N1_%s'%channel,'N1_%s'%channel,'@0/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)
		N2		= RooFormulaVar( 'N2_%s'%channel,'N2_%s'%channel,'@1*(@2/12)**2/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)
		N3		= RooFormulaVar( 'N3_%s'%channel,'N3_%s'%channel,'@3*(@4/20)**2/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)
		N4		= RooFormulaVar( 'N4_%s'%channel,'N4_%s'%channel,'@5*(@6/60)**2/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)

		N_list		= RooArgList(N1,N2,N3,N4)
		Pdf_list	= RooArgList(SMPdf,wtmp.pdf('Pdf_quad_%s_%s'%(POI[0],channel)),wtmp.pdf('Pdf_quad_%s_%s'%(POI[1],channel)),wtmp.pdf('Pdf_quad_%s_%s'%(POI[2],channel)))
		model		= RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list)
		model.Print()

		scale_list	= RooArgList(wtmp.function('scaleshape_%s_%s'%(POI[0],channel)),\
						wtmp.function('scaleshape_%s_%s'%(POI[1],channel)),\
						wtmp.function('scaleshape_%s_%s'%(POI[2],channel)))
		normfactor_3d	= RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list)
		getattr(WS,'import')(normfactor_3d)	
		getattr(wtmp,'import')(normfactor_3d)	
		wtmp.Print()

		#fit 3 pdfs
		for i in range(3):
			for j in range(3):
				wtmp.var(POI[j]).setVal(0)
			wtmp.var(POI[i]).setVal(par_max[POI[i]])
			wtmp.var('a_quad_4fit_%s_%s_%s'%(POI[i],cat,ch)).setConstant(kFALSE)
			fitres		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
			fitresults.append(fitres)
			wtmp.var('a_quad_4fit_%s_%s_%s'%(POI[i],cat,ch)).setConstant(kTRUE)

		for i in range(4):
			fitresults[i].Print()
	
		model.Print()	
		for i in range(3):
			wtmp.var(POI[i]).setVal(0)
		getattr(wtmp,'import')(model)

	#read, rename and write bkg pdfs and bkg rates

	fileInWs	= TFile.Open('Input/wwlvj_%s_HP%s_workspace.root'%(ch[:2],cat[1]))
	w		= fileInWs.Get('workspace4limit_') 
	w.pdf('STop_xww_%s_HP%s'%(ch,cat[1])).SetName('STop')
	w.pdf('TTbar_xww_%s_HP%s'%(ch,cat[1])).SetName('TTbar')
	w.pdf('WJets_xww_%s_HP%s'%(ch,cat[1])).SetName('WJets')
	w.var('rrv_number_VV_xww_%s_mj_prefit_signal_region'%ch).SetName('rate_VV')
	w.var('rrv_number_STop_xww_%s_mj_postfit_signal_region'%ch).SetName('rate_STop') 
	w.var('ttbar_norm_signal_region_postfit').SetName('rate_TTbar')
	w.var('rrv_number_WJets0_xww_%s_mj_postfit_signal_region'%ch).SetName('rate_WJets')
	getattr(WS,'import')(w.pdf('STop'))
	getattr(WS,'import')(w.pdf('TTbar'))
	getattr(WS,'import')(w.pdf('WJets'))	
	getattr(WS,'import')(w.var('rate_VV'))
	getattr(WS,'import')(w.var('rate_STop'))
	getattr(WS,'import')(w.var('rate_TTbar'))
	getattr(WS,'import')(w.var('rate_WJets'))




	path	='%s/src/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'%os.environ["CMSSW_BASE"]
	output 	= TFile('%s/%s.root'%(path,channel),'recreate')
	
	data_obs_hist	= TH1F('data_obs_hist','data_obs_hist',nbins,binlo,binhi)
	data_obs_tree	= TTree('data_obs_tree','data_obs_tree')
	MWW_data	= array('f',[1])
	data_obs_tree.Branch('observable',MWW_data,'observable')
	data_obs_hist.Print("all")

	for i in range(tree_data.GetEntries()):
	  	tree_data.GetEntry(i)
		if tree_data.Mjpruned<mj_hi and tree_data.Mjpruned>mj_lo and tree_data.MWW > binlo and tree_data.MWW < binhi:
			MWW_data[0] = tree_data.MWW 
			data_obs_tree.Fill()
			data_obs_hist.Fill(MWW_data[0])
	data_obs_tree.Write()
	data_obs_hist.Write()

	WS.SetName('w')
	WS.Write();
	output.Close()
	print 'Write to file ' + output.GetName()


	#make plots
	if options.make_plots:
		make_plots(rrv_mass_lvj,wtmp,ch,cat,channel,fitresults,600,binhi,binlo4fit)
		#make_plots(rrv_mass_lvj,wtmp,ch,cat,channel,fitresults,binlo,binhi,binlo4fit)
	if options.make_plots2:
		#cross check
		canvas		= TCanvas(POI[0]+','+POI[1]+','+POI[2] , POI[0]+','+POI[1]+','+POI[2],1)
		pad11	= TPad('pad1_all','pad1_all',0.,0.25,1.,1.)
		pad22	= TPad('pad2_all','pad2_all',0.,0.02,1.,0.25)
		canvas.cd()
		pad11.Draw()
		pad22.Draw()
		pad11.SetLogy()
		pad11.cd()
		
		rrv_mass_lvj.Print()
		rrv_mass_lvj.setRange('all',600,3500)
		print model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),'all').getVal()
		rrv_mass_lvj.setRange(600,3500)
		p4		= rrv_mass_lvj.frame(600,3500)

		#wtmp.data('SMdatahist').plotOn(p4,RooFit.MarkerColor(kBlack),RooFit.MarkerSize(0.75))

		#datahist_all3.plotOn(p4,RooFit.LineColor(kBlue),RooFit.LineWidth(2),RooFit.DrawOption('E'))
		#for i in range(3):
		#		wtmp.var(POI[i]).setVal(0)
		#model.plotOn(p4,RooFit.LineColor(kBlack),RooFit.Normalization(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*wtmp.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))
		for i in range(3):
			wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		normval 	= RooRealVar('normval','normval',wtmp.function('normfactor_3d_4fit_%s'%channel).getVal()*wtmp.data('SMdatahist4fit').sumEntries())
		datahist_all3.plotOn(p4,RooFit.LineColor(kBlue),RooFit.LineWidth(2),RooFit.DrawOption('E'))
		model.plotOn(p4,RooFit.LineColor(kBlue),RooFit.Normalization(normval.getVal(),RooAbsReal.NumEvent))
		wtmp.data('SMdatahist4fit').plotOn(p4,RooFit.LineColor(kBlack),RooFit.LineWidth(2),RooFit.DrawOption('E'))

		#normval = wtmp.function('normfactor_3d_4fit_%s'%channel).getVal()*wtmp.data('SMdatahist4fit').sumEntries()
		print str(normval.getVal()) + ' / ' + str(datahist_all3.sumEntries())
		print str((normval.getVal()/datahist_all3.sumEntries() -1)*100) + ' %'
		print wtmp.function('normfactor_3d_4fit_%s'%channel).getVal(), ' , ' , wtmp.data('SMdatahist4fit').sumEntries(), ' , ' ,datahist_all3.sumEntries() 
		
		plotmin = 0.5
		#if signal_category == 'WZ':
		#	plotmin = 1e-3
		p4.GetYaxis().SetRangeUser(plotmin,p4.GetMaximum()*2)
		#p4.GetXaxis().SetRangeUser(900,3500)
		#canvas.SetLogy()
		pullhist2 = p4.pullHist('h_datahist_all3','aTGC_model_%s_%s_Norm[rrv_mass_lvj]'%(cat,ch))
		pad22.cd()
		data_dif = TH1D('data_dif2','data_dif2',26,900,3500)
		data_dif.Sumw2()
		data_dif.SetMarkerStyle(21)
		data_dif.SetMaximum(3)
		data_dif.SetMinimum(-3)
		data_dif.GetYaxis().SetNdivisions(7)
		data_dif.GetYaxis().SetTitle('#frac{MC-Fit}{error}')
		data_dif.GetYaxis().SetLabelSize(0.125)
		data_dif.GetYaxis().SetTitleSize(0.05)
		data_dif.Draw("")
		pullhist2.SetLineColor(kBlue)
		pullhisttmp.SetLineColor(kRed)
		pullhist2.Draw("SAME E1")

		canvas.Update()
		canvas.SaveAs('docuplots/atgc3_%s.pdf'%channel)
		canvas.SaveAs('docuplots/atgc3_%s.png'%channel)
		
		raw_input('cross check plot plotted')
		canvas.Close()
	
	for i in range(len(fitresults)):
		fitresults[i].Print()
	if options.noatgcint:
		fileOut4plot = TFile.Open('docuplots/make_plots/%s_noatgcint.root'%channel,'RECREATE')
	else:
		fileOut4plot = TFile.Open('docuplots/make_plots/%s.root'%channel,'RECREATE')
	wtmp.Write()
	fileOut4plot.Close()

	if options.printatgc:
		wtmp.var('cwww').setVal(12)
		wtmp.var('ccw').setVal(20)
		wtmp.var('cb').setVal(0)
		print 'cwww and ccw positive:'
		for i in range(8):
			print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())
		wtmp.var('cwww').setVal(12)
		wtmp.var('ccw').setVal(0)
		wtmp.var('cb').setVal(60)
		print 'cwww and cb positive:'
		for i in range(8):
			print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())
		wtmp.var('cwww').setVal(0)
		wtmp.var('ccw').setVal(20)
		wtmp.var('cb').setVal(60)
		print 'ccw and cb positive:'
		for i in range(8):
			print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())
		raw_input(channel)


if options.cat=='WWWZ':
	if options.chan=='elmu':
		make_input('el','WW',int(options.mlvj_lo),int(options.mlvj_hi))
		make_input('mu','WW',int(options.mlvj_lo),int(options.mlvj_hi))
		make_input('el','WZ',int(options.mlvj_lo),int(options.mlvj_hi))
		make_input('mu','WZ',int(options.mlvj_lo),int(options.mlvj_hi))
	else:
		make_input(options.chan,'WW',int(options.mlvj_lo),int(options.mlvj_hi))
		make_input(options.chan,'WZ',int(options.mlvj_lo),int(options.mlvj_hi))
else:
	if options.chan=='elmu':
		make_input('el',options.cat,int(options.mlvj_lo),int(options.mlvj_hi))
		make_input('mu',options.cat,int(options.mlvj_lo),int(options.mlvj_hi))
	else:
		make_input(options.chan,options.cat,int(options.mlvj_lo),int(options.mlvj_hi))

