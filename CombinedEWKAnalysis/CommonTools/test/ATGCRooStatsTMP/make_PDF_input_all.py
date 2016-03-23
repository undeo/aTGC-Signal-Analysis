from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import random

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

POI	=	['cwww','ccw','cb']
par_max = {'cwww' : 12, 'ccw' : 20, 'cb' : 60}

parser	= OptionParser()
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='make new input trees')
parser.add_option('-p', '--plots', action='store_true', dest='do_plots', default=False, help='make plots')
parser.add_option('--plots2', action='store_true', dest='do_plots2', default=False, help='make parabel plot')
parser.add_option('--cat', dest='cat', default='WW', help='category, WW or WZ')
(options,args) = parser.parse_args()


if options.cat == 'WW':
	sigreg 	= 'lo'
	mj_lo	= 65.
	mj_hi	= 85.
elif options.cat == 'WZ':
	sigreg	= 'hi'
	mj_lo	= 85.
	mj_hi	= 105
else:
	raise RuntimeError('cateogry not supported!')

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def make_ATGCtree(ch='ele'):

	for categ in ['WW','WZ']:
		fileInATGC	= TFile.Open('Input/%s-aTGC-%s.root'%(categ,ch))
		treeInATGC	= fileInATGC.Get('treeDumper/BasicTree')
		fileOutATGC	= TFile('Output/ATGC-Tree_%s_%s_%s.root'%(categ,sigreg,ch), 'recreate')
		treeATGC_tmp 	= treeInATGC.CloneTree(0)

		weight		= array('f',[1])
		branch_weight	= treeATGC_tmp.Branch('weight',weight,'weight')
		c_wwwl		= array('f',[1])
		branch_c_wwwl	= treeATGC_tmp.Branch('c_wwwl',c_wwwl,'c_wwwl')
		c_wl		= array('f',[1])
		branch_c_wl	= treeATGC_tmp.Branch('c_wl',c_wl,'c_wl')	
		c_bl		= array('f',[1])
		branch_c_bl	= treeATGC_tmp.Branch('c_bl',c_bl,'c_bl')

		lumi_tmp 	= 2093.917403402
		for i in range(treeInATGC.GetEntries()):
			if i%10000==0:
				print i
			treeInATGC.GetEntry(i)
			if treeInATGC.jet_pt>200. and treeInATGC.jet_tau2tau1<0.6 \
			and treeInATGC.Mjpruned<mj_hi and treeInATGC.Mjpruned>mj_lo \
			and treeInATGC.W_pt>200. and treeInATGC.deltaR_LeptonWJet>math.pi/2. \
			and abs(treeInATGC.deltaPhi_WJetMet)>2. and abs(treeInATGC.deltaPhi_WJetWlep)>2.\
			and treeInATGC.MWW>1000 and treeInATGC.MWW<3500:
				weight_part	= 1/20. * lumi_tmp * treeInATGC.puweight * (treeInATGC.genweight/abs(treeInATGC.genweight))
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
	mergetree.Add('Output/ATGC-Tree_WW_%s_%s.root'%(sigreg,ch))
	mergetree.Add('Output/ATGC-Tree_WZ_%s_%s.root'%(sigreg,ch))
	newtree		= mergetree.Clone()
	newtree.Write()
	newtree.Print()

	print '--------> Write to file ' + fileOut.GetName()
	
	

def make_plots(rrv_x,wtmp,ch):
        
        can             = []
        plots	        = []
	
	for i in range(len(POI)):
                c       = TCanvas(POI[i],POI[i],1)
                p       = rrv_x.frame()
                can.append(c)
                plots.append(p)
	for i in range(len(POI)):
		for j in range(len(POI)):
			wtmp.var(POI[j]).setVal(0)
		wtmp.data('SMdatahist').plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		wtmp.data('neg_datahist_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		normvalSM	= wtmp.function('normfactor_%sd_%s'%(len(POI),ch)).getVal() * wtmp.data('SMdatahist').sumEntries()
		wtmp.pdf('aTGC_model').plotOn(plots[i],\
					      RooFit.LineColor(kBlack),\
					      RooFit.Normalization(normvalSM, RooAbsReal.NumEvent))
		wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		normval		= wtmp.function('normfactor_%sd_%s'%(len(POI),ch)).getVal() * wtmp.data('SMdatahist').sumEntries()
        	wtmp.pdf('aTGC_model').plotOn(plots[i],\
					      RooFit.LineColor(kBlue),\
					      RooFit.Normalization(normval, RooAbsReal.NumEvent))
		for j in range(len(POI)):
			wtmp.var(POI[j]).setVal(0)
		#for j in range(15):
		#	wtmp.var(POI[i]).setVal(-par_max[POI[i]]/10 * j+2)
		#	normval	= wtmp.function('normfactor_%sd_%s'%(len(POI),ch)).getVal() * wtmp.data('SMdatahist').sumEntries()
		#	wtmp.pdf('aTGC_model').plotOn(plots[i],\
		#				      RooFit.LineColor(kGray+1),\
		#				      RooFit.LineWidth(1),\
		#				      RooFit.Normalization(normval,RooAbsReal.NumEvent))
		wtmp.data('SMdatahist').plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		wtmp.data('neg_datahist_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		if ch == 'ele':
			plotmin = 1e-2
			plotmax = 50
			if options.cat == 'WZ':
				plotmin = 3e-4
		elif ch == 'mu':
			plotmin = 1e-2
			plotmax = 1e2
			if options.cat == 'WZ':
				plotmin = 1e-3
				plotmax = 50
		plots[i].GetYaxis().SetRangeUser(plotmin,plotmax)
		can[i].cd()
		can[i].SetLogy()	
		plots[i].SetTitle('')
		plots[i].Draw()
		can[i].Update()
		can[i].SaveAs('docuplots/%s_%s_%s.pdf'%(POI[i],options.cat,ch))
	raw_input('plots plotted')



def make_input(ch = 'ele'):

	cat		= options.cat	
	nbins4fit	= 20
	binlo		= 1000
	binhi		= 3500
	WS		= RooWorkspace("WS")

	#read data
	fileInData	= TFile.Open('Input/treeEDBR_data_xww_%s.root'%(ch))
	tree_tmp	= fileInData.Get('tree')  
	data_obs = TH1F('data_obs','data_obs',nbins4fit,binlo,binhi)
	for i in range(tree_tmp.GetEntries()):
	  	tree_tmp.GetEntry(i)
	    	#add cuts to data!
		#data still blinded!
		data_obs.Fill(random.random())

	#read or make ATGC and new tree
	if options.newtrees:
		make_ATGCtree(ch)
	fileInATGC	= TFile.Open('Output/ATGC-Tree_%s_%s.root'%(sigreg,ch))
	treeATGC	= fileInATGC.Get('BasicTree')

	#prepare variables, parameters and temporary workspace
	wtmp		= RooWorkspace('wtmp')
	a1		= RooRealVar('a_SM_%s'%ch,'a_SM_%s'%ch,-0.1,-2,0)
	if 'cwww' in POI:	
		cwww		= RooRealVar('cwww','cwww',0,-120,120); 	cwww.setConstant(kTRUE);	getattr(wtmp,'import')(cwww);
		cwww12		= RooFormulaVar('cwww12','cwww12','(@0/12)**2',RooArgList(cwww));		getattr(wtmp,'import')(cwww12);
	if 'ccw' in POI:
		ccw		= RooRealVar('ccw','ccw',0,-200,200);		ccw.setConstant(kTRUE);		getattr(wtmp,'import')(ccw);
		ccw20		= RooFormulaVar('ccw20','ccw20','(@0/20)**2',RooArgList(ccw));			getattr(wtmp,'import')(ccw20);
	if 'cb' in POI:
		cb		= RooRealVar('cb','cb',0,-600,600);		cb.setConstant(kTRUE);		getattr(wtmp,'import')(cb);	
		cb60		= RooFormulaVar('cb60','cb60','(@0/60)**2',RooArgList(cb));			getattr(wtmp,'import')(cb60);

	#make and fill SM histogram, SM fit
	SMhist		= TH1F('SMhist','SMhist',nbins4fit,binlo,binhi)
	SMhist.Sumw2(kTRUE)

	fileInWs	= TFile.Open('Input/wwlvj_BulkG_WW_lvjj_M800_%s_HPW_workspace.root'%ch[:2])
	w		= fileInWs.Get('workspace4limit_')
	rrv_mass_lvj	= w.var('rrv_mass_lvj')
	rrv_mass_lvj.setRange(binlo,binhi)

	SMPdf		= RooExponential('SMPdf','SMPdf',rrv_mass_lvj,a1)
	for i in range(treeATGC.GetEntries()):
		treeATGC.GetEntry(i)
		if treeATGC.c_wwwl == 0 and treeATGC.c_wl ==0 and treeATGC.c_bl == 0:
	    		SMhist.Fill(treeATGC.MWW,treeATGC.weight)
	SMdatahist	= RooDataHist('SMdatahist','SMdatahist',RooArgList(rrv_mass_lvj),SMhist)
	fitresSM	= SMPdf.fitTo(SMdatahist, RooFit.SumW2Error(kTRUE))
	a1.setConstant(kTRUE)
	N_SM		= RooRealVar('N_SM','N_SM',SMdatahist.sumEntries())
	N_SM.setConstant(kTRUE)
	getattr(wtmp,'import')(SMdatahist)
	getattr(wtmp,'import')(N_SM)

	#make and fill ATGC histograms, fit quadratic scales

	for para in POI:
		hist4fit	= TH1F('hist4fit','hist4fit',3,-1.5*par_max[para],1.5*par_max[para])
		pos_hist	= TH1F('c_pos_hist','c_pos_hist',nbins4fit,binlo,binhi)
		pos_hist.Sumw2(kTRUE)
		neg_hist	= TH1F('c_neg_hist','c_neg_hist',nbins4fit,binlo,binhi)
		neg_hist.Sumw2(kTRUE)

		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (para == 'cwww' and treeATGC.c_wwwl == par_max[para] and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
			or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == par_max[para] and treeATGC.c_bl == 0)\
			or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == par_max[para]):
		      		pos_hist.Fill(treeATGC.MWW,treeATGC.weight)
		    	if (para == 'cwww' and treeATGC.c_wwwl == -par_max[para] and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
			or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == -par_max[para] and treeATGC.c_bl == 0)\
			or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == -par_max[para]):
		      		neg_hist.Fill(treeATGC.MWW,treeATGC.weight)

		pos_datahist	= RooDataHist('pos_datahist_%s'%para,'pos_datahist_%s'%para,RooArgList(rrv_mass_lvj),pos_hist)
		neg_datahist	= RooDataHist('neg_datahist_%s'%para,'neg_datahist_%s'%para,RooArgList(rrv_mass_lvj),neg_hist)

		hist4fit.SetBinContent(1,neg_datahist.sumEntries()/N_SM.getVal())
		hist4fit.SetBinContent(2,1)
		hist4fit.SetBinContent(3,pos_datahist.sumEntries()/N_SM.getVal())
		dummy		= TCanvas('dummy','dummy',1)
		dummy.cd()
		hist4fit.Fit('pol2','OQ')
		dummy.Close()
		fitfunc		= hist4fit.GetFunction('pol2')
		par0		= RooRealVar('par0_%s_%s'%(para,ch),'par0_%s_%s'%(para,ch),fitfunc.GetParameter(0)); 		par0.setConstant(kTRUE);
		par1		= RooRealVar('par1_%s_%s'%(para,ch),'par1_%s_%s'%(para,ch),fitfunc.GetParameter(1)); 		par1.setConstant(kTRUE);
		par2		= RooRealVar('par2_%s_%s'%(para,ch),'par2_%s_%s'%(para,ch),fitfunc.GetParameter(2)); 		par2.setConstant(kTRUE);
		N1_ratio	= RooRealVar('N1_ratio_%s_%s'%(para,ch),'N1_ratio_%s_%s'%(para,ch),neg_datahist.sumEntries()/SMdatahist.sumEntries()-1)		
		scaleshape	= RooFormulaVar('scaleshape_%s_%s'%(para,ch),'scaleshape_%s_%s'%(para,ch),\
						'@0+@1*@3+@2*@3**2',\
						RooArgList(par0,par1,par2,wtmp.var(para)))			
		a2		= RooRealVar('a_%s_%s'%(para,ch),'a_%s_%s'%(para,ch),-0.001,-2,0)
		a2.setConstant(kTRUE)
		cPdf		= RooExponential('Pdf_%s'%para,'Pdf_%s'%para,rrv_mass_lvj,a2)
		getattr(wtmp,'import')(cPdf)
		getattr(wtmp,'import')(scaleshape)
		getattr(wtmp,'import')(N1_ratio)
		getattr(wtmp,'import')(pos_datahist)
		getattr(wtmp,'import')(neg_datahist)
		wtmp.Print()
		

	#make model

	paralist	= RooArgList(wtmp.function('N1_ratio_%s_%s'%(POI[0],ch)),wtmp.function(POI[0]+str(par_max[POI[0]])),\
				     wtmp.function('N1_ratio_%s_%s'%(POI[1],ch)),wtmp.function(POI[1]+str(par_max[POI[1]])),\
				     wtmp.function('N1_ratio_%s_%s'%(POI[2],ch)),wtmp.function(POI[2]+str(par_max[POI[2]])))
	N1_3d		= RooFormulaVar('N1_3d','N1_3d','1/(1+@0*@1+@2*@3+@4*@5)',paralist)
	N2_3d		= RooFormulaVar('N2_3d','N2_3d','(@0*@1)/(1+@0*@1+@2*@3+@4*@5)',paralist)
	N3_3d		= RooFormulaVar('N3_3d','N3_3d','(@2*@3)/(1+@0*@1+@2*@3+@4*@5)',paralist)	
	N4_3d		= RooFormulaVar('N4_3d','N4_3d','(@4*@5)/(1+@0*@1+@2*@3+@4*@5)',paralist)
	N_list		= RooArgList(N1_3d,N2_3d,N3_3d,N4_3d)
	Pdf_list	= RooArgList(SMPdf,wtmp.pdf('Pdf_%s'%POI[0]),wtmp.pdf('Pdf_%s'%POI[1]),wtmp.pdf('Pdf_%s'%POI[2]))
	model		= RooAddPdf('aTGC_model','aTGC_model', Pdf_list, N_list)

	scale_list	= RooArgList(wtmp.function('scaleshape_%s_%s'%(POI[0],ch)),\
					wtmp.function('scaleshape_%s_%s'%(POI[1],ch)),\
					wtmp.function('scaleshape_%s_%s'%(POI[2],ch)))
	normfactor_3d	= RooFormulaVar('normfactor_3d_%s'%ch,'normfactor_3d_%s'%ch,'1+(@0-1)+(@1-1)+(@2-1)',scale_list)
					
	getattr(WS,'import')(normfactor_3d)	
	getattr(wtmp,'import')(normfactor_3d)	
	model.Print()
	getattr(wtmp,'import')(model)
	wtmp.Print()
	
	#fit 3 pdfs
	fitresults	= []
	for i in range(3):
		wtmp.var(POI[0]).setVal(0); wtmp.var(POI[1]).setVal(0); wtmp.var(POI[2]).setVal(0);
		wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		wtmp.var('a_%s_%s'%(POI[i],ch)).setConstant(kFALSE)
		fitres		= model.fitTo(wtmp.data('neg_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		fitresults.append(fitres)
		wtmp.var('a_%s_%s'%(POI[i],ch)).setConstant(kTRUE)
	for i in range(3):
		fitresults[i].Print()
	
	if options.do_plots:
		make_plots(rrv_mass_lvj,wtmp,ch)
	if options.do_plots2:
		#cross check
		canvas		= TCanvas(POI[0]+','+POI[1]+','+POI[2],POI[0]+','+POI[1]+','+POI[2],1)
		p4		= rrv_mass_lvj.frame()
		hist_all3	= TH1F('hist_all3','hist_all3',nbins4fit,binlo,binhi)
		hist_all3.Sumw2(kTRUE)
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (treeATGC.c_wwwl == -12 and treeATGC.c_wl == -20 and treeATGC.c_bl == -60):
		      		hist_all3.Fill(treeATGC.MWW,treeATGC.weight)
		datahist_all3	= RooDataHist('datahist_all3','datahist_all3',RooArgList(rrv_mass_lvj),hist_all3)
		wtmp.data('SMdatahist').plotOn(p4,RooFit.MarkerColor(kBlack),RooFit.MarkerSize(0.75))

		datahist_all3.plotOn(p4,RooFit.MarkerColor(kBlue),RooFit.MarkerSize(0.75))
		for i in range(3):
				wtmp.var(POI[i]).setVal(0)
		model.plotOn(p4,RooFit.LineColor(kBlack),RooFit.Normalization(wtmp.function('normfactor_3d_%s'%ch).getVal()*wtmp.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))
		for i in range(3):
				wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		model.plotOn(p4,RooFit.LineColor(kBlue),RooFit.Normalization(wtmp.function('normfactor_3d_%s'%ch).getVal()*wtmp.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))

		wtmp.data('SMdatahist').plotOn(p4,RooFit.MarkerColor(kBlack),RooFit.MarkerSize(0.75))
		datahist_all3.plotOn(p4,RooFit.MarkerColor(kBlue),RooFit.MarkerSize(0.75))

		for i in range(3):
			wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		normval = wtmp.function('normfactor_3d_%s'%ch).getVal()*wtmp.data('SMdatahist').sumEntries()
		print str(normval) + ' / ' + str(datahist_all3.sumEntries())
		print str((normval/datahist_all3.sumEntries() -1)*100) + ' %'

		canvas.cd()
		plotmin = 1e-2
		if options.cat == 'WZ':
			plotmin = 1e-3
		p4.GetYaxis().SetRangeUser(plotmin,p4.GetMaximum()*3)
		canvas.SetLogy()
		p4.Draw()
		canvas.SetLogy
		canvas.Update()
		canvas.SaveAs('docuplots/atgc3_%s_%s.pdf'%(options.cat,ch))
		raw_input('cross check plot plotted')


	#read, rename and write bkg pdfs and bkg rates

	fileInWs	= TFile.Open('Input/wwlvj_BulkG_WW_lvjj_M800_%s_HP%s_workspace.root'%(ch[:2],cat[1]))
	w		= fileInWs.Get('workspace4limit_') 
	w.pdf('STop_xww_%s_HP%s'%(ch[:2],cat[1])).SetName('STop')
	w.pdf('TTbar_xww_%s_HP%s'%(ch[:2],cat[1])).SetName('TTbar')
	w.pdf('WJets_xww_%s_HP%s'%(ch[:2],cat[1])).SetName('WJets')
	w.var('rate_VV_xww_for_unbin').SetName('rate_VV')
	w.var('rate_STop_xww_for_unbin').SetName('rate_STop') 
	w.var('rate_TTbar_xww_for_unbin').SetName('rate_TTbar')
	w.var('rate_WJets_xww_for_unbin').SetName('rate_WJets')
	getattr(WS,'import')(w.pdf('STop'))
	getattr(WS,'import')(w.pdf('TTbar'))
	getattr(WS,'import')(w.pdf('WJets'))	
	getattr(WS,'import')(w.var('rate_VV'))
	getattr(WS,'import')(w.var('rate_STop'))
	getattr(WS,'import')(w.var('rate_TTbar'))
	getattr(WS,'import')(w.var('rate_WJets'))

	#import to workspace

	getattr(WS,'import')(model)

	path	='/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/src/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'
	output 	= TFile('%s/ch_%s_%s.root'%(path,cat,ch),'recreate')

	data_obs.Write()
	WS.SetName('workspace'); WS.Print(); WS.Write();
	output.Close()
	print 'Write to file ' + output.GetName()

  
make_input('ele')
make_input('mu')
