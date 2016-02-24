from ROOT import  *
from array import array

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

#POI	=	['cwww']
POI	=	['cwww','ccw']

def make_SMBKG_input(ch = 'ele', POIs = ['cwww']):
	do_plot		= 1 ; oldtrees	= 0; old = 0;
	nbins4fit	= 20
	binlo		= 1000
	binhi		= 3500
	#read bkg
	fileInWS	= TFile.Open('Input/wwlvj_BulkG_WW_lvjj_M800_%s_HPW_workspace.root'%(ch))
	w		= fileInWS.Get('workspace4limit_') 
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
	if oldtrees:
		fileInATGC	= TFile.Open('Output/ATGC-Tree_%s'%ch)
		treeATGC	= fileInATGC.Get('BasicTree')
	else:
		fileInATGC	= TFile.Open('Input/WW-aTGC-%s.root'%ch)
		treeInATGC	= fileInATGC.Get('treeDumper/BasicTree')
		fileOutATGC	= TFile('Output/ATGC-Tree_%s.root'%ch, 'recreate')
		treeInATGC.SetBranchStatus('*',0)
		treeATGC 	= treeInATGC.CloneTree(0)
		treeInATGC.SetBranchStatus('aTGCWeights',1); treeInATGC.SetBranchStatus('m_lvj',1); treeInATGC.SetBranchStatus('PUweight',1);
		weight		= array('f',[1]);	branch_weight	= treeATGC.Branch('weight',weight,'weight');
		m_lvj		= array('f',[1]);	branch_m_lvj	= treeATGC.Branch('m_lvj',m_lvj,'m_lvj');
		c_wwwl		= array('f',[1]);	branch_c_wwwl	= treeATGC.Branch('c_wwwl',c_wwwl,'c_wwwl')
		c_wl		= array('f',[1]);	branch_c_wl	= treeATGC.Branch('c_wl',c_wl,'c_wl')	
		c_bl		= array('f',[1]);	branch_c_bl	= treeATGC.Branch('c_bl',c_bl,'c_bl')

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
			treeATGC.Fill()
			if 'cwww' in POIs:
				c_wwwl[0] = 12; c_wl[0]	= 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[0] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = -12; c_wl[0] = 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[1] * weight_part
				treeATGC.Fill()
			if 'ccw' in POIs:
				c_wwwl[0] = 0; c_wl[0] = 20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[2] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = 0; c_wl[0] = -20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[3] * weight_part
				treeATGC.Fill()
			if 'cb' in POIs:
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = 60;
				weight[0] = treeInATGC.aTGCWeights[4] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = -60;
				weight[0] = treeInATGC.aTGCWeights[5] * weight_part
				treeATGC.Fill()
		treeATGC.Write()
		print '--------> Write to file ' + fileOutATGC.GetName()


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
	ws		= RooWorkspace('ws')
	getattr(ws,'import')(SMdatahist)
	for para in POIs:
		hist4fit	= TH1F('hist4fit','hist4fit',3,-1.5*par_max[para],1.5*par_max[para])
		pos_hist	= TH1F('c_pos_hist','c_pos_hist',nbins4fit,binlo,binhi)
		neg_hist	= TH1F('c_neg_hist','c_neg_hist',nbins4fit,binlo,binhi)
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (para == 'cwww' and treeATGC.c_wwwl == par_max[para])\
			or (para == 'ccw' and treeATGC.c_wl == par_max[para])\
			or (para == 'cb' and treeATGC.c_bl == par_max[para]):
		      		pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
		    	if (para == 'cwww' and treeATGC.c_wwwl == -par_max[para])\
			or (para == 'ccw' and treeATGC.c_wl == -par_max[para])\
			or (para == 'cb' and treeATGC.c_bl == -par_max[para]):
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
		#elif para == 'cb':
		#	poi	= cb60

		a2		= RooRealVar('a_%s_%s'%(para,ch),'a_%s_%s'%(para,ch),-0.001,-2,0)
		cwwwPdf		= RooExponential('Pdf_%s'%para,'Pdf_%s'%para,rrv_mass_lvj,a2)
		getattr(ws,'import')(cwwwPdf)
		getattr(ws,'import')(scaleshape)
		getattr(ws,'import')(N1);getattr(ws,'import')(N2);
		getattr(ws,'import')(pos_datahist);getattr(ws,'import')(neg_datahist)
		ws.Print()



	#make and fit model
	if len(POIs) == 1:
		ws.var(POIs[0]).setVal(par_max[POIs[0]])
		print ws.function('N1_%s'%POIs[0]).getVal()
		model 		= RooAddPdf('aTGC_model','aTGC_model',RooArgList(SMPdf,ws.pdf('Pdf_%s'%POIs[0])),RooArgList(ws.function('N1_%s'%POIs[0]),ws.function('N2_%s'%POIs[0])))			
		fitres		= model.fitTo(ws.data('pos_datahist_%s'%POIs[0]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		ws.var('a_%s_%s'%(POIs[0],ch)).setConstant(kTRUE)
		fitres.Print()
		if do_plot:
			can1		= TCanvas('canvas1',POIs[0],1)
			p1 		= rrv_mass_lvj.frame()
			ws.data('pos_datahist_%s'%POIs[0]).plotOn(p1,RooFit.MarkerColor(par_max[POIs[0]]+1))
			model.plotOn(p1,RooFit.LineColor(par_max[POIs[0]]+1))
			ws.data('SMdatahist').plotOn(p1)
			p1.GetYaxis().SetRangeUser(0.03,75)
			p1.Draw()
			can1.SetLogy()
			raw_input('@@')


	if len(POIs) == 2:
		if do_plot:
			can1		= TCanvas('canvas1',POIs[0],1)
			can2		= TCanvas('canvas2',POIs[1],1)
			p1 		= rrv_mass_lvj.frame()
			p2 		= rrv_mass_lvj.frame()
		if 'cwww' in POIs:
			cwww12	= RooFormulaVar('cwww12','cwww12','(@0/12)**2',RooArgList(ws.var('cwww')))
		if 'ccw' in POIs:
			ccw20	= RooFormulaVar('ccw20','ccw20','(@0/20)**2',RooArgList(ws.var('ccw')))
		paralist	= RooArgList(ws.function('N1_ratio_%s'%POIs[0]),cwww12,ws.function('N1_ratio_%s'%POIs[1]),ccw20)
		N1_2d		= RooFormulaVar('N1_2d','N1_2d','1/(1+@0*@1+@2)',paralist)
		N2_2d		= RooFormulaVar('N2_2d','N2_2d','(@0*@1)/(1+@0*@1+@2*@3)',paralist)
		N3_2d		= RooFormulaVar('N3_2d','N3_2d','(@2*@3)/(1+@0*@1+@2*@3)',paralist)		
		model		= RooAddPdf('aTGC_model','aTGC_model',RooArgList(SMPdf,ws.pdf('Pdf_%s'%POIs[0]),ws.pdf('Pdf_%s'%POIs[1])),RooArgList(N1_2d,N2_2d,N3_2d))
		model.Print()
		#fit first pdf
		ws.var(POIs[0]).setVal(par_max[POIs[0]])
		ws.var(POIs[1]).setVal(0)
		ws.var('a_%s_%s'%(POIs[1],ch)).setConstant(kTRUE)
		fitres1		= model.fitTo(ws.data('pos_datahist_%s'%POIs[0]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		ws.var('a_%s_%s'%(POIs[0],ch)).setConstant(kTRUE)
		if do_plot:
			ws.data('pos_datahist_%s'%POIs[0]).plotOn(p1,RooFit.MarkerColor(par_max[POIs[0]]+1))
			model.plotOn(p1,RooFit.LineColor(par_max[POIs[0]]+1))
		ws.var('a_%s_%s'%(POIs[1],ch)).setConstant(kFALSE)
		#fit second pdf
		ws.var(POIs[1]).setVal(par_max[POIs[1]])
		ws.var(POIs[0]).setVal(0)
		fitres2		= model.fitTo(ws.data('pos_datahist_%s'%POIs[1]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		ws.var('a_%s_%s'%(POIs[1],ch)).setConstant(kTRUE);
		fitres1.Print()
		fitres2.Print()
		if do_plot:
			ws.data('pos_datahist_%s'%POIs[1]).plotOn(p2,RooFit.MarkerColor(par_max[POIs[1]]+1))
			model.plotOn(p2,RooFit.LineColor(par_max[POIs[1]]+1))
			#add SM to plot
			ws.data('SMdatahist').plotOn(p1,RooFit.MarkerColor(kBlack))
			ws.data('SMdatahist').plotOn(p2,RooFit.MarkerColor(kBlack))
			ws.var(POIs[0]).setVal(0)
			ws.var(POIs[1]).setVal(0)
			model.plotOn(p1,RooFit.LineColor(kBlack))
			model.plotOn(p2,RooFit.LineColor(kBlack))
			for i in range(par_max[POIs[0]]):
				ws.var(POIs[0]).setVal(i);		
				ws.var(POIs[1]).setVal(0)
				normval0	= ws.data('SMdatahist').sumEntries()*ws.function('scaleshape_%s'%POIs[0]).getVal()
				model.plotOn(p1,RooFit.LineWidth(2),RooFit.LineColor(i+1),RooFit.Normalization(normval0,RooAbsReal.NumEvent))
				ws.data('SMdatahist').plotOn(p1,RooFit.MarkerColor(kBlack))
				ws.data('pos_datahist_%s'%POIs[0]).plotOn(p1,RooFit.MarkerColor(par_max[POIs[0]]+1))
			for i in range(par_max[POIs[1]]):
				ws.var(POIs[0]).setVal(0);		
				ws.var(POIs[1]).setVal(i)
				normval1	= ws.data('SMdatahist').sumEntries()*ws.function('scaleshape_%s'%POIs[1]).getVal()
				model.plotOn(p2,RooFit.LineWidth(2),RooFit.LineColor(i+1),RooFit.Normalization(normval1,RooAbsReal.NumEvent))
				ws.data('SMdatahist').plotOn(p2,RooFit.MarkerColor(kBlack))
				ws.data('pos_datahist_%s'%POIs[1]).plotOn(p2,RooFit.MarkerColor(par_max[POIs[1]]+1))

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

	if old:

		#import to ws
		getattr(w,'import')(model)
	

		path	='../../CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'
		cwww 	= '_cwww' if 'cwww' in POIs else ''
		ccw = '_ccw' if 'ccw' in POIs else '' 
		cb = '_cb' if 'cb' in POIs else ''
		output 	= TFile('%s/ch_%s%s%s%s.root'%(path,ch,cwww,ccw,cb),'recreate')

		data_obs.Write()
		w.SetName('workspace'); w.Print(); w.Write();
		output.Close()
		print 'Write to file ' + output.GetName()

  
make_SMBKG_input('ele',POI)
#make_SMBKG_input('mu',POI)
#make_SMBKG_input('ele',POI)
#make_SMBKG_input('mu',POI)
