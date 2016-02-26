from ROOT import  *
from array import array
from optparse import OptionParser

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

POI	=	[]

parser	= OptionParser()
parser.add_option('--POIs', dest='parameters', help='define parameters of interest')
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='make new input trees')
parser.add_option('-p', '--plots', action='store_true', dest='do_plots', default=False, help='make plots')
(options,args) = parser.parse_args()
for parname in options.parameters.split(','):
	POI.append(parname)


def make_ATGCtree(fileInATGC='', ch='ele'):

	fileInATGC	= TFile.Open(fileInATGC)
	treeInATGC	= fileInATGC.Get('treeDumper/BasicTree')
	fileOutATGC	= TFile('Output/ATGC-Tree_%s.root'%ch, 'recreate')
	treeInATGC.SetBranchStatus('*',0)
	treeATGC_tmp 	= treeInATGC.CloneTree(0)
	treeInATGC.SetBranchStatus('aTGCWeights',1); treeInATGC.SetBranchStatus('m_lvj',1); treeInATGC.SetBranchStatus('PUweight',1);
	weight		= array('f',[1]);	branch_weight	= treeATGC_tmp.Branch('weight',weight,'weight');
	m_lvj		= array('f',[1]);	branch_m_lvj	= treeATGC_tmp.Branch('m_lvj',m_lvj,'m_lvj');
	c_wwwl		= array('f',[1]);	branch_c_wwwl	= treeATGC_tmp.Branch('c_wwwl',c_wwwl,'c_wwwl');
	c_wl		= array('f',[1]);	branch_c_wl	= treeATGC_tmp.Branch('c_wl',c_wl,'c_wl');	
	c_bl		= array('f',[1]);	branch_c_bl	= treeATGC_tmp.Branch('c_bl',c_bl,'c_bl');

	lumi_tmp 	= 2093.917403402
	for i in range(treeInATGC.GetEntries()):
		if i%10000==0:
			print i
		treeInATGC.GetEntry(i)
		m_lvj[0]	= treeInATGC.m_lvj
		weight_part	= 1/20. * lumi_tmp * treeInATGC.PUweight
		c_wwwl[0]	= 0
		c_wl[0]		= 0
		c_bl[0]		= 0
		weight[0]	= treeInATGC.aTGCWeights[7] * weight_part
		treeATGC_tmp.Fill()
		if 'cwww' in POI:
			c_wwwl[0] = 12; c_wl[0]	= 0; c_bl[0] = 0;
			weight[0] = treeInATGC.aTGCWeights[0] * weight_part
			treeATGC_tmp.Fill()
			c_wwwl[0] = -12; c_wl[0] = 0; c_bl[0] = 0;
			weight[0] = treeInATGC.aTGCWeights[1] * weight_part
			treeATGC_tmp.Fill()
		if 'ccw' in POI:
			c_wwwl[0] = 0; c_wl[0] = 20; c_bl[0] = 0;
			weight[0] = treeInATGC.aTGCWeights[2] * weight_part
			treeATGC_tmp.Fill()
			c_wwwl[0] = 0; c_wl[0] = -20; c_bl[0] = 0;
			weight[0] = treeInATGC.aTGCWeights[3] * weight_part
			treeATGC_tmp.Fill()
		if 'cb' in POI:
			c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = 60;
			weight[0] = treeInATGC.aTGCWeights[4] * weight_part
			treeATGC_tmp.Fill()
			c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = -60;
			weight[0] = treeInATGC.aTGCWeights[5] * weight_part
			treeATGC_tmp.Fill()
		if 'cwww' in POI and 'ccw' in POI and 'cb' in POI:
			c_wwwl[0] = 12; c_wl[0] = 20; c_bl[0] = 60;
			weight[0] = treeInATGC.aTGCWeights[6] * weight_part
			treeATGC_tmp.Fill()
		
	treeATGC_tmp.Write()
	treeATGC_tmp.Print()
	print '--------> Write to file ' + fileOutATGC.GetName()
	return treeATGC_tmp	


def make_SMBKG_input(ch = 'ele'):

	nbins4fit	= 20
	binlo		= 1000
	binhi		= 3500
	WS		= RooWorkspace("WS")

	#read bkg
	fileInWs	= TFile.Open('Input/wwlvj_BulkG_WW_lvjj_M800_%s_HPW_workspace.root'%(ch))
	w		= fileInWs.Get('workspace4limit_') 
	w.pdf('STop_xww_%s_HPW'%(ch[:2])).SetName('STop')
	w.pdf('TTbar_xww_%s_HPW'%(ch[:2])).SetName('TTbar')
	w.pdf('WJets_xww_%s_HPW'%(ch[:2])).SetName('WJets')
	rrv_mass_lvj	= w.var('rrv_mass_lvj')
	rrv_mass_lvj.setRange(binlo,binhi) 

	#read data
	fileInData	= TFile.Open('Input/treeEDBR_data_xww_%s.root'%(ch))
	tree_tmp	= fileInData.Get('tree')  
	data_obs = TH1F('data_obs','data_obs',nbins4fit,binlo,binhi)
	for i in range(tree_tmp.GetEntries()):
	  	tree_tmp.GetEntry(i)
	    	data_obs.Fill(tree_tmp.m_lvj)

	#read or make ATGC and new tree
	if options.newtrees:
		treeATGC	= make_ATGCtree('Input/WW-aTGC-%s.root'%ch, ch)
		treeATGC.Print()
	else:
		fileInATGC	= TFile.Open('Output/ATGC-Tree_%s.root'%ch)
		treeATGC	= fileInATGC.Get('BasicTree')
		treeATGC.Print()

	#prepare fit
	a1		= RooRealVar('a_SM_%s'%ch,'a_SM_%s'%ch,-0.1,-2,0)	
	cwww		= RooRealVar('cwww','cwww',0,-12,12); 		cwww.setConstant(kTRUE);
	ccw		= RooRealVar('ccw','ccw',0,-20,20);		ccw.setConstant(kTRUE);
	cb		= RooRealVar('cb','cb',0,-60,60);		cb.setConstant(kTRUE);

	#make and fill SM histogram, SM fit
	SMhist		= TH1F('SMhist','SMhist',nbins4fit,binlo,binhi)
	SMPdf		= RooExponential('SMPdf','SMPdf',rrv_mass_lvj,a1)
	for i in range(treeATGC.GetEntries()):
		treeATGC.GetEntry(i)
		if treeATGC.c_wwwl == 0 and treeATGC.c_wl ==0 and treeATGC.c_bl == 0:
	    		SMhist.Fill(treeATGC.m_lvj,treeATGC.weight)
	SMdatahist	= RooDataHist('SMdatahist','SMdatahist',RooArgList(rrv_mass_lvj),SMhist)
	SMPdf.fitTo(SMdatahist, RooFit.SumW2Error(kTRUE));	a1.setConstant(kTRUE)

	#make and fill ATGC histograms
	par_max		= {'cwww' : 12, 'ccw' : 20, 'cb' : 60}
	wtmp		= RooWorkspace('wtmp')
	getattr(wtmp,'import')(SMdatahist)

	for para in POI:
		hist4fit	= TH1F('hist4fit','hist4fit',3,-1.5*par_max[para],1.5*par_max[para])
		pos_hist	= TH1F('c_pos_hist','c_pos_hist',nbins4fit,binlo,binhi)
		neg_hist	= TH1F('c_neg_hist','c_neg_hist',nbins4fit,binlo,binhi)
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (para == 'cwww' and treeATGC.c_wwwl == par_max[para] and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
			or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == par_max[para] and treeATGC.c_bl == 0)\
			or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == par_max[para]):
		      		pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
		    	if (para == 'cwww' and treeATGC.c_wwwl == -par_max[para] and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
			or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == -par_max[para] and treeATGC.c_bl == 0)\
			or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == -par_max[para]):
		      		neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
		pos_datahist	= RooDataHist('pos_datahist_%s'%para,'pos_datahist_%s'%para,RooArgList(rrv_mass_lvj),pos_hist)
		neg_datahist	= RooDataHist('neg_datahist_%s'%para,'neg_datahist_%s'%para,RooArgList(rrv_mass_lvj),neg_hist)
		hist4fit.SetBinContent(1,neg_datahist.sumEntries()/SMdatahist.sumEntries())
		hist4fit.SetBinContent(2,1)
		hist4fit.SetBinContent(3,pos_datahist.sumEntries()/SMdatahist.sumEntries())
		hist4fit.Fit('pol2')
		fitfunc		= hist4fit.GetFunction('pol2')
		par0		= RooRealVar('par0_%s_%s'%(para,ch),'par0_%s_%s'%(para,ch),fitfunc.GetParameter(0)); 		par0.setConstant(kTRUE);
		par1		= RooRealVar('par1_%s_%s'%(para,ch),'par1_%s_%s'%(para,ch),fitfunc.GetParameter(1)); 		par1.setConstant(kTRUE);
		par2		= RooRealVar('par2_%s_%s'%(para,ch),'par2_%s_%s'%(para,ch),fitfunc.GetParameter(2)); 		par2.setConstant(kTRUE);
		N1_ratio	= RooRealVar('N1_ratio_%s'%para,'N1_ratio_%s'%para,pos_datahist.sumEntries()/SMdatahist.sumEntries()-1)
		if para == 'cwww':
			N1		= RooFormulaVar('N1_%s'%para,'N1_%s'%para,'1/(1+@0*(@1/12)**2)',RooArgList(N1_ratio,cwww))
			N2		= RooFormulaVar('N2_%s'%para,'N2_%s'%para,'(@0*(@1/12)**2)/(1+@0*(@1/12)**2)',RooArgList(N1_ratio,cwww))
			scaleshape	= RooFormulaVar('scaleshape_%s'%para,'scaleshape%s'%para,'@0+@1*@3+@2*@3*@3',RooArgList(par0,par1,par2,cwww))
		elif para == 'ccw':
			N1		= RooFormulaVar('N1_%s'%para,'N1_%s'%para,'1/(1+@0*(@1/20)**2)',RooArgList(N1_ratio,ccw))
			N2		= RooFormulaVar('N2_%s'%para,'N2_%s'%para,'(@0*(@1/20)**2)/(1+@0*(@1/20)**2)',RooArgList(N1_ratio,ccw))
			scaleshape	= RooFormulaVar('scaleshape_%s'%para,'scaleshape%s'%para,'@0+@1*@3+@2*@3*@3',RooArgList(par0,par1,par2,ccw))
		elif para == 'cb':
			N1		= RooFormulaVar('N1_%s'%para,'N1_%s'%para,'1/(1+@0*(@1/60)**2)',RooArgList(N1_ratio,cb))
			N2		= RooFormulaVar('N2_%s'%para,'N2_%s'%para,'(@0*(@1/60)**2)/(1+@0*(@1/60)**2)',RooArgList(N1_ratio,cb))
			scaleshape	= RooFormulaVar('scaleshape_%s'%para,'scaleshape%s'%para,'@0+@1*@3+@2*@3*@3',RooArgList(par0,par1,par2,cb))			
		a2		= RooRealVar('a_%s_%s'%(para,ch),'a_%s_%s'%(para,ch),-0.001,-2,0)
		a2.setConstant(kTRUE)
		cwwwPdf		= RooExponential('Pdf_%s'%para,'Pdf_%s'%para,rrv_mass_lvj,a2)
		getattr(wtmp,'import')(cwwwPdf)
		getattr(wtmp,'import')(scaleshape)
		getattr(WS,'import')(scaleshape)
		getattr(wtmp,'import')(N1);getattr(wtmp,'import')(N2);
		getattr(wtmp,'import')(pos_datahist);getattr(wtmp,'import')(neg_datahist)
		wtmp.Print()

	#make and fit model
	if len(POI) == 1:
		wtmp.var(POI[0]).setVal(par_max[POI[0]])
		print wtmp.function('N1_%s'%POI[0]).getVal()
		model 		= RooAddPdf('aTGC_model','aTGC_model',RooArgList(SMPdf,wtmp.pdf('Pdf_%s'%POI[0])),RooArgList(wtmp.function('N1_%s'%POI[0]),wtmp.function('N2_%s'%POI[0])))
		wtmp.var('a_%s_%s'%(POI[0],ch)).setConstant(kFALSE)			
		fitres		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[0]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		wtmp.var('a_%s_%s'%(POI[0],ch)).setConstant(kTRUE)
		fitres.Print()
		if options.do_plots:
			can1		= TCanvas('canvas1',POI[0],1)
			p1 		= rrv_mass_lvj.frame()
			wtmp.data('pos_datahist_%s'%POI[0]).plotOn(p1,RooFit.MarkerColor(kRed))
			model.plotOn(p1,RooFit.LineColor(kRed))
			wtmp.data('SMdatahist').plotOn(p1)
			p1.GetYaxis().SetRangeUser(0.03,75)
			p1.Draw()
			can1.SetLogy()
			raw_input('@@')

	if len(POI) == 2:
		if options.do_plots:
			can1		= TCanvas('canvas1',POI[0],1)
			can2		= TCanvas('canvas2',POI[1],1)
			p1 		= rrv_mass_lvj.frame()
			p2 		= rrv_mass_lvj.frame()
		if 'cwww' in POI:
			cwww12	= RooFormulaVar('cwww12','cwww12','(@0/12)**2',RooArgList(wtmp.var('cwww')))
		if 'ccw' in POI:
			ccw20	= RooFormulaVar('ccw20','ccw20','(@0/20)**2',RooArgList(wtmp.var('ccw')))
		paralist	= RooArgList(wtmp.function('N1_ratio_%s'%POI[0]),cwww12,wtmp.function('N1_ratio_%s'%POI[1]),ccw20)
		N1_2d		= RooFormulaVar('N1_2d','N1_2d','1/(1+@0*@1+@2*@3)',paralist)
		N2_2d		= RooFormulaVar('N2_2d','N2_2d','(@0*@1)/(1+@0*@1+@2*@3)',paralist)
		N3_2d		= RooFormulaVar('N3_2d','N3_2d','(@2*@3)/(1+@0*@1+@2*@3)',paralist)		
		model		= RooAddPdf('aTGC_model','aTGC_model',RooArgList(SMPdf,wtmp.pdf('Pdf_%s'%POI[0]),wtmp.pdf('Pdf_%s'%POI[1])),RooArgList(N1_2d,N2_2d,N3_2d))
		model.Print()
		#fit first pdf
		wtmp.var(POI[0]).setVal(par_max[POI[0]])
		wtmp.var(POI[1]).setVal(0)
		wtmp.var('a_%s_%s'%(POI[0],ch)).setConstant(kFALSE)
		fitres1		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[0]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		wtmp.var('a_%s_%s'%(POI[0],ch)).setConstant(kTRUE)
		if options.do_plots:
			wtmp.data('pos_datahist_%s'%POI[0]).plotOn(p1,RooFit.MarkerColor(kRed))
			model.plotOn(p1,RooFit.LineColor(kRed))
		#fit second pdf
		wtmp.var(POI[1]).setVal(par_max[POI[1]])
		wtmp.var(POI[0]).setVal(0)
		wtmp.var('a_%s_%s'%(POI[1],ch)).setConstant(kFALSE)
		fitres2		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[1]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		wtmp.var('a_%s_%s'%(POI[1],ch)).setConstant(kTRUE);
		fitres1.Print()
		fitres2.Print()
		if options.do_plots:
			wtmp.data('pos_datahist_%s'%POI[1]).plotOn(p2,RooFit.MarkerColor(kRed))
			model.plotOn(p2,RooFit.LineColor(kRed))
			#add SM to plot
			wtmp.data('SMdatahist').plotOn(p1,RooFit.MarkerColor(kBlack))
			wtmp.data('SMdatahist').plotOn(p2,RooFit.MarkerColor(kBlack))
			wtmp.var(POI[0]).setVal(0)
			wtmp.var(POI[1]).setVal(0)
			model.plotOn(p1,RooFit.LineColor(kBlack))
			model.plotOn(p2,RooFit.LineColor(kBlack))
			for i in range(par_max[POI[0]]):
				wtmp.var(POI[0]).setVal(i);		
				wtmp.var(POI[1]).setVal(0)
				normval0	= wtmp.data('SMdatahist').sumEntries()*wtmp.function('scaleshape_%s'%POI[0]).getVal()
				model.plotOn(p1,RooFit.LineWidth(2),RooFit.LineColor(i+1),RooFit.Normalization(normval0,RooAbsReal.NumEvent))
				wtmp.data('SMdatahist').plotOn(p1,RooFit.MarkerColor(kBlack))
				wtmp.data('pos_datahist_%s'%POI[0]).plotOn(p1,RooFit.MarkerColor(kRed))
			for i in range(par_max[POI[1]]):
				wtmp.var(POI[0]).setVal(0);		
				wtmp.var(POI[1]).setVal(i)
				normval1	= wtmp.data('SMdatahist').sumEntries()*wtmp.function('scaleshape_%s'%POI[1]).getVal()
				model.plotOn(p2,RooFit.LineWidth(2),RooFit.LineColor(i+1),RooFit.Normalization(normval1,RooAbsReal.NumEvent))
				wtmp.data('SMdatahist').plotOn(p2,RooFit.MarkerColor(kBlack))
				wtmp.data('pos_datahist_%s'%POI[1]).plotOn(p2,RooFit.MarkerColor(kRed))
			#draw plot
			can1.cd()
			p1.GetYaxis().SetRangeUser(0.03,75)
			p1.Draw()
			can1.SetLogy()
			can1.Update()
			can2.cd()
			p2.GetYaxis().SetRangeUser(0.03,75)
			p2.Draw()
			can2.SetLogy()
			can2.Update()
			raw_input('...')

	if len(POI) == 3:
		if options.do_plots:
			can1		= TCanvas('canvas1',POI[0],1)
			can2		= TCanvas('canvas2',POI[1],1)
			can3		= TCanvas('canvas3',POI[2],1)
			p1 		= rrv_mass_lvj.frame()
			p2 		= rrv_mass_lvj.frame()
			p3 		= rrv_mass_lvj.frame()
			plots		= [p1,p2,p3]
		if 'cwww' in POI:
			cwww12	= RooFormulaVar('cwww12','cwww12','(@0/12)**2',RooArgList(wtmp.var('cwww')))
		if 'ccw' in POI:
			ccw20	= RooFormulaVar('ccw20','ccw20','(@0/20)**2',RooArgList(wtmp.var('ccw')))
		if 'cb' in POI:
			cb60	= RooFormulaVar('cb60','cb60','(@0/60)**2',RooArgList(wtmp.var('cb')))
		paralist	= RooArgList(wtmp.function('N1_ratio_%s'%POI[0]),cwww12,wtmp.function('N1_ratio_%s'%POI[1]),ccw20,wtmp.function('N1_ratio_%s'%POI[2]),cb60)
		N1_3d		= RooFormulaVar('N1_3d','N1_3d','1/(1+@0*@1+@2*@3+@4*@5)',paralist)
		N2_3d		= RooFormulaVar('N2_3d','N2_3d','(@0*@1)/(1+@0*@1+@2*@3+@4*@5)',paralist)
		N3_3d		= RooFormulaVar('N3_3d','N3_3d','(@2*@3)/(1+@0*@1+@2*@3+@4*@5)',paralist)	
		N4_3d		= RooFormulaVar('N4_3d','N4_3d','(@4*@5)/(1+@0*@1+@2*@3+@4*@5)',paralist)	
		model		= RooAddPdf('aTGC_model','aTGC_model',RooArgList(SMPdf,wtmp.pdf('Pdf_%s'%POI[0]),wtmp.pdf('Pdf_%s'%POI[1]),wtmp.pdf('Pdf_%s'%POI[2])),RooArgList(N1_3d,N2_3d,N3_3d,N4_3d))
		model.Print()

		#fit 3 pdfs
		fitresults	= []
		for i in range(3):
			wtmp.var(POI[0]).setVal(0); wtmp.var(POI[1]).setVal(0); wtmp.var(POI[2]).setVal(0);
			wtmp.var(POI[i]).setVal(par_max[i])
			wtmp.var('a_%s_%s'%(POI[i],ch)).setConstant(kFALSE)
			fitres		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
			fitresults.append(fitres)
			wtmp.var('a_%s_%s'%(POI[i],ch)).setConstant(kTRUE)
			if options.do_plots:
				wtmp.data('pos_datahist_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kRed))
				model.plotOn(plots[i],RooFit.LineColor(kRed))
		for j in range(3):
			fitresults[j].Print()


	#import to workspace
	#getattr(w,'import')(model)
	getattr(WS,'import')(model)
	getattr(WS,'import')(w.pdf('STop'))
	getattr(WS,'import')(w.pdf('TTbar'))
	getattr(WS,'import')(w.pdf('WJets'))	
	getattr(WS,'import')(w.var('rate_VV_xww_for_unbin'))
	getattr(WS,'import')(w.var('rate_STop_xww_for_unbin'))
	getattr(WS,'import')(w.var('rate_TTbar_xww_for_unbin'))
	getattr(WS,'import')(w.var('rate_WJets_xww_for_unbin'))

	path	='../../CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'
	scwww 	= '_cwww' if 'cwww' in POI else ''
	sccw	= '_ccw' if 'ccw' in POI else '' 
	scb 	= '_cb' if 'cb' in POI else ''
	output 	= TFile('%s/ch_%s%s%s%s.root'%(path,ch,scwww,sccw,scb),'recreate')

	data_obs.Write()
	WS.SetName('workspace'); WS.Print(); WS.Write();
	output.Close()
	print 'Write to file ' + output.GetName()


  
make_SMBKG_input('ele')
make_SMBKG_input('mu')
