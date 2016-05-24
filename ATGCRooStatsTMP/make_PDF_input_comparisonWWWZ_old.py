from ROOT import  *
from array import array
from optparse import OptionParser
import math as math

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

POI	=	[]
par_max = {'cwww' : 12, 'ccw' : 20, 'cb' : 60}

parser	= OptionParser()
parser.add_option('--POIs', dest='parameters', help='define parameters of interest')
parser.add_option('-p', '--plots', action='store_true', dest='do_plots', default=False, help='make plots')
parser.add_option('--cat', dest='cat', default='WW', help='category, WW or WZ')
(options,args) = parser.parse_args()
for parname in options.parameters.split(','):
	POI.append(parname)

if options.cat == 'WW':
	sigreg = 'lo'
elif options.cat == 'WZ':
	sigreg = 'hi'

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)





def make_input(ch = 'el'):

	cat		= options.cat
	nbins4fit	= 20
	binlo		= 1000
	binhi		= 3500
	WS		= RooWorkspace("WS")

	#read bkg
	fileInWs	= TFile.Open('Input/wwlvj_%s_HP%s_workspace.root'%(ch[:2],cat[1]))
	w		= fileInWs.Get('workspace4limit_') 
	w.pdf('STop_xww_%s_HP%s'%(ch[:2],cat[1])).SetName('STop')
	w.pdf('TTbar_xww_%s_HP%s'%(ch[:2],cat[1])).SetName('TTbar')
	w.pdf('WJets_xww_%s_HP%s'%(ch[:2],cat[1])).SetName('WJets')
	rrv_mass_lvj	= w.var('rrv_mass_lvj')
	rrv_mass_lvj.setRange(binlo,binhi) 

	#read data
	#fileInData	= TFile.Open('Input/treeEDBR_data_xww_%s.root'%(ch))
	#tree_tmp	= fileInData.Get('tree')  
	data_obs = TH1F('data_obs','data_obs',nbins4fit,binlo,binhi)
	#for i in range(tree_tmp.GetEntries()):
	#  	tree_tmp.GetEntry(i)
	    	#data_obs.Fill(tree_tmp.m_lvj)
		#data still blinded!
	#	data_obs.Fill(i)

	wtmp		= RooWorkspace('wtmp')
	if 'cwww' in POI:	
		cwww		= RooRealVar('cwww','cwww',0,-120,120); 	cwww.setConstant(kTRUE);	getattr(wtmp,'import')(cwww);
		cwww12		= RooFormulaVar('cwww12','cwww12','(@0/12)**2',RooArgList(cwww));		getattr(wtmp,'import')(cwww12);
	if 'ccw' in POI:
		ccw		= RooRealVar('ccw','ccw',0,-200,200);		ccw.setConstant(kTRUE);		getattr(wtmp,'import')(ccw);
		ccw20		= RooFormulaVar('ccw20','ccw20','(@0/20)**2',RooArgList(ccw));			getattr(wtmp,'import')(ccw20);
	if 'cb' in POI:
		cb		= RooRealVar('cb','cb',0,-600,600);		cb.setConstant(kTRUE);		getattr(wtmp,'import')(cb);	
		cb60		= RooFormulaVar('cb60','cb60','(@0/60)**2',RooArgList(cb));			getattr(wtmp,'import')(cb60);


	#read or make ATGC and new tree


	SMhist		= TH1F('SMhist','SMhist',nbins4fit,binlo,binhi)
	SMhist.Sumw2(kTRUE)
	SMhistWW	= TH1F('SMhist_WW','SMhist_WW',nbins4fit,binlo,binhi)
	SMhistWW.Sumw2(kTRUE)
	SMhistWZ	= TH1F('SMhist_WZ','SMhist_WZ',nbins4fit,binlo,binhi)
	SMhistWZ.Sumw2(kTRUE)
	for categ in ['WW','WZ']:
		fileInATGC	= TFile.Open('Output/ATGC-Tree_%s_%s_%s.root'%(categ,sigreg,ch))
		treeATGC	= fileInATGC.Get('BasicTree')
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
			if treeATGC.c_wwwl == 0 and treeATGC.c_wl ==0 and treeATGC.c_bl == 0:
		    		SMhist.Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if categ == 'WW':
					SMhistWW.Fill(treeATGC.MWW,treeATGC.totEventWeight)
				if categ == 'WZ':
					SMhistWZ.Fill(treeATGC.MWW,treeATGC.totEventWeight)
	SMdatahist	= RooDataHist('SMdatahist%s'%sigreg,'SMdatahist%s'%sigreg,RooArgList(rrv_mass_lvj),SMhist)
	SMdatahistWW	= RooDataHist('SMdatahist%s_WW'%sigreg,'SMdatahist%s_WW'%sigreg,RooArgList(rrv_mass_lvj),SMhistWW)
	SMdatahistWZ	= RooDataHist('SMdatahist%s_WZ'%sigreg,'SMdatahist%s_WZ'%sigreg,RooArgList(rrv_mass_lvj),SMhistWZ)
	N_SM		= RooRealVar('N_SM','N_SM',SMdatahist.sumEntries())
	N_SM.setConstant(kTRUE)
	N_SM_WW		= RooRealVar('N_SM_WW','N_SM_WW',SMdatahistWW.sumEntries())
	N_SM_WW.setConstant(kTRUE)
	N_SM_WZ		= RooRealVar('N_SM_WZ','N_SM_WZ',SMdatahistWZ.sumEntries())
	N_SM_WZ.setConstant(kTRUE)
	getattr(wtmp,'import')(N_SM_WW)
	getattr(wtmp,'import')(N_SM_WZ)
	getattr(wtmp,'import')(N_SM)
	getattr(wtmp,'import')(SMdatahist)


	a1		= RooRealVar('a_%s'%(ch),'a_%s'%(ch),-0.1,-2,0)
	SMPdf		= RooExponential('SMPdf','SMPdf',rrv_mass_lvj,a1)

	getattr(wtmp,'import')(SMPdf)



	for para in POI:
		gROOT.SetBatch(kTRUE)
		neg_hist_all	= TH1F('c_neg_hist_all','c_neg_hist_all',nbins4fit,binlo,binhi)
		neg_hist_all.Sumw2(kTRUE)
		pos_hist_all	= TH1F('c_pos_hist_all','c_pos_hist_all',nbins4fit,binlo,binhi)
		pos_hist_all.Sumw2(kTRUE)
		for categ in ['WW','WZ']:
			fileInATGC	= TFile.Open('Output/ATGC-Tree_%s_%s_%s.root'%(categ,sigreg,ch))
			treeATGC	= fileInATGC.Get('BasicTree')
			for i in range(treeATGC.GetEntries()):
				treeATGC.GetEntry(i)
			    	if (para == 'cwww' and treeATGC.c_wwwl == 12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
				or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 20 and treeATGC.c_bl == 0)\
				or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == 60):
			      		pos_hist_all.Fill(treeATGC.MWW,treeATGC.totEventWeight)
			    	if (para == 'cwww' and treeATGC.c_wwwl == -12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
				or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == -20 and treeATGC.c_bl == 0)\
				or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == -60):
			      		neg_hist_all.Fill(treeATGC.MWW,treeATGC.totEventWeight)
		pos_datahist_all= RooDataHist('pos_datahist_all_%s'%(para),'pos_datahist_all_%s'%(para),RooArgList(rrv_mass_lvj),pos_hist_all)
		neg_datahist_all= RooDataHist('neg_datahist_all_%s'%(para),'neg_datahist_all_%s'%(para),RooArgList(rrv_mass_lvj),neg_hist_all)
		getattr(wtmp,'import')(pos_datahist_all)
		getattr(wtmp,'import')(neg_datahist_all)

		hist4fit	= TH1F('hist4fit','hist4fit',3,-1.5*par_max[para],1.5*par_max[para])
		hist4fit.SetBinContent(1,wtmp.data('neg_datahist_all_%s'%para).sumEntries()/wtmp.var('N_SM').getVal())
		hist4fit.SetBinContent(2,1)
		hist4fit.SetBinContent(3,wtmp.data('pos_datahist_all_%s'%para).sumEntries()/wtmp.var('N_SM').getVal())
		hist4fit.Fit('pol2','OQ')
		fitfunc		= hist4fit.GetFunction('pol2')
		p0		= RooRealVar('p0_%s_%s'%(para,ch),'p0_%s_%s'%(para,ch),fitfunc.GetParameter(0))
 		p0.setConstant(kTRUE)
		p1		= RooRealVar('p1_%s_%s'%(para,ch),'p1_%s_%s'%(para,ch),fitfunc.GetParameter(1))
 		p1.setConstant(kTRUE)
		p2		= RooRealVar('p2_%s_%s'%(para,ch),'p2_%s_%s'%(para,ch),fitfunc.GetParameter(2))
		p2.setConstant(kTRUE)
		scaleshape	= RooFormulaVar('scaleshape_%s_%s'%(para,ch),'scaleshape_%s_%s'%(para,ch),\
						'@0+@1*@3+@2*@3**2',\
						RooArgList(p0,p1,p2,wtmp.var(para)))
		getattr(wtmp,'import')(scaleshape)
		gROOT.SetBatch(kFALSE)


	for categ in ['WW','WZ']:
		fileInATGC	= TFile.Open('Output/ATGC-Tree_%s_%s_%s.root'%(categ,sigreg,ch))
		treeATGC	= fileInATGC.Get('BasicTree')
		for para in POI:
			neg_hist	= TH1F('c_neg_hist','c_neg_hist',nbins4fit,binlo,binhi)
			neg_hist.Sumw2(kTRUE)
			for i in range(treeATGC.GetEntries()):
				treeATGC.GetEntry(i)
			    	if (para == 'cwww' and treeATGC.c_wwwl == -12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
				or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == -20 and treeATGC.c_bl == 0)\
				or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == -60):
			      		neg_hist.Fill(treeATGC.MWW,treeATGC.totEventWeight)
			neg_datahist	= RooDataHist('neg_datahist_%s_%s'%(categ,para),'neg_datahist_%s_%s'%(categ,para),RooArgList(rrv_mass_lvj),neg_hist)
			N_diff		= RooRealVar('N_diff_%s'%(categ),'N_diff_%s'%(categ),neg_datahist.sumEntries() - wtmp.var('N_SM_%s'%categ).getVal())		
			b2		= RooRealVar('b_%s_%s_%s'%(categ,para,ch),'b_%s_%s_%s'%(categ,para,ch),-0.001,-2,0)
			b2.setConstant(kTRUE)
			cPdf		= RooExponential('Pdf_%s_%s'%(categ,para),'Pdf_%s_%s'%(categ,para),rrv_mass_lvj,b2)
			getattr(wtmp,'import')(cPdf)
			getattr(wtmp,'import')(N_diff)
			getattr(wtmp,'import')(neg_datahist)
			neg_hist.Reset()
		

	wtmp.Print()

	#make and fit model
	if len(POI) == 1:
		 
		normfactor_1d 	= RooFormulaVar('normfactor_1d_%s'%ch,'normfactor_1d_%s'%ch,\
						'1+(@0-1)',\
						RooArgList(wtmp.function('scaleshape_%s_%s'%(POI[0],ch))))


		neg_datahist_all= wtmp.data('neg_datahist_all_%s'%POI[0])
		paralist	= RooArgList(wtmp.var('N_SM'),\
						wtmp.function('N_diff_WW'),\
						wtmp.function('N_diff_WZ'),\
						wtmp.function(POI[0]+str(par_max[POI[0]])))
		N1_1d		= RooFormulaVar('N1_1d','N1_1d','@0/(@0+@1*@3+@2*@3)',paralist)	
		N2_1d		= RooFormulaVar('N2_1d','N2_1d','@1*@3/(@0+@1*@3+@2*@3)',paralist)
		N3_1d		= RooFormulaVar('N3_1d','N3_1d','@2*@3/(@0+@1*@3+@2*@3)',paralist)
		model 		= RooAddPdf('aTGC_model','aTGC_model',\
					    RooArgList(wtmp.pdf('SMPdf'),wtmp.pdf('Pdf_WW_%s'%POI[0]),wtmp.pdf('Pdf_WZ_%s'%POI[0])),\
					    RooArgList(N1_1d,N2_1d,N3_1d) )
		#SM fit
		wtmp.var(POI[0]).setVal(0)
		fitresSM	= model.fitTo(wtmp.data('SMdatahist%s'%sigreg),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		wtmp.var('a_%s'%ch).setConstant(kTRUE)
		plot = rrv_mass_lvj.frame()
		wtmp.data('SMdatahist%s'%sigreg).plotOn(plot,RooFit.MarkerColor(kBlue))
		model.plotOn(plot)

		#atgc fit
		wtmp.var(POI[0]).setVal(-par_max[POI[0]])
		wtmp.var('b_WW_%s_%s'%(POI[0],ch)).setConstant(kFALSE)
		wtmp.var('b_WZ_%s_%s'%(POI[0],ch)).setConstant(kFALSE)
		fitres2		= model.fitTo(neg_datahist_all, RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		wtmp.var('b_WW_%s_%s'%(POI[0],ch)).setConstant(kTRUE)
		wtmp.var('b_WZ_%s_%s'%(POI[0],ch)).setConstant(kTRUE)
		neg_datahist_all.plotOn(plot, RooFit.MarkerColor(kRed))

		normval 	= normfactor_1d.getVal() * N_SM.getVal()
		model.plotOn(plot, RooFit.LineColor(kRed), RooFit.Normalization(normval,RooAbsReal.NumEvent))

		wtmp.data('neg_datahist_WW_%s'%POI[0]).plotOn(plot,RooFit.MarkerColor(kGreen+1),RooFit.MarkerStyle(25))
		SMdatahistWW.plotOn(plot,RooFit.MarkerColor(kGreen+1),RooFit.MarkerStyle(26))
		wtmp.data('neg_datahist_WZ_%s'%POI[0]).plotOn(plot,RooFit.MarkerColor(kOrange+1),RooFit.MarkerStyle(25))
		SMdatahistWZ.plotOn(plot,RooFit.MarkerColor(kOrange+1),RooFit.MarkerStyle(26))

		leg = TLegend(0.75,0.75,0.9,0.9)
		leg.AddEntry(plot.getObject(0),'SM','p')
		leg.AddEntry(plot.getObject(1),'SM-Model','l')
		leg.AddEntry(plot.getObject(2),'atgc','p')
		leg.AddEntry(plot.getObject(3),'atgc-model','l')
		leg.AddEntry(plot.getObject(4),'WW-atgc','p')
		leg.AddEntry(plot.getObject(5),'WW-SM','p')
		leg.AddEntry(plot.getObject(6),'WZ-atgc','p')
		leg.AddEntry(plot.getObject(7),'WZ-SM','p')

		plot.remove('aTGC_model')
		plot.remove('aTGC_model')

		plot.addObject(leg)

		fitresSM.Print()
		fitres2.Print()

		model.Print()
		wtmp.pdf('Pdf_WW_%s'%POI[0]).Print()
		wtmp.pdf('Pdf_WZ_%s'%POI[0]).Print()
		print normval
		print neg_datahist_all.sumEntries()
		print str((normval / neg_datahist_all.sumEntries() -1) * 100) + '%'
	

		c1 = TCanvas('c1','c1',1)
		c1.cd()
		c1.SetLogy()
		plot.GetYaxis().SetRangeUser(1e-4,100)
		plot.Draw()
		c1.Update()

		raw_input("S")
		exit(0)

		wtmp.var('b_WW_%s_%s'%(POI[0],ch)).setConstant(kFALSE)
		wtmp.var(POI[0]).setVal(-par_max[POI[0]])			
		fitresWW	= model.fitTo(wtmp.data('neg_datahist_WW_%s'%POI[0]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		wtmp.var('b_WW_%s_%s'%(POI[0],ch)).setConstant(kTRUE)



		fitresWW.Print()

                getattr(WS,'import')(normfactor_1d)
		getattr(wtmp,'import')(normfactor_1d)
		getattr(wtmp,'import')(model)
		if options.do_plots:
			make_plots(rrv_mass_lvj,wtmp,ch)

	if len(POI) == 2:
		raise RuntimeError('not supported atm')

	if len(POI) == 3:
		raise RuntimeError('not supported atm')



  
make_input('el')
make_input('mu')
