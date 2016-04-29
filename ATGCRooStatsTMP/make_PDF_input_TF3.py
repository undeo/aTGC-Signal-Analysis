from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import random

gSystem.Load('/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so')

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

POI	=	['cwww','ccw','cb']
par_max = {'cwww' : 12, 'ccw' : 20, 'cb' : 60}#atgc points

parser	= OptionParser()
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='recreate input trees')
parser.add_option('-p', '--plots', action='store_true', dest='do_plots', default=False, help='make plots')
parser.add_option('--plots2', action='store_true', dest='do_plots2', default=False, help='make parabel plot')
parser.add_option('--cat', dest='cat', default='WW', help='category, WW or WZ, defines signal region')
parser.add_option('--lin', action='store_true', dest='lin', default=False, help='add linear dependency')
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

binlo 	= 900.
binhi	= 3500.

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def make_ATGCtree(ch='el'):

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
			and treeInATGC.MWW>binlo and treeInATGC.MWW<binhi\
			and treeInATGC.nbtag==0:
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
	newtree.Print()

	print '--------> Write to file ' + fileOut.GetName()
	
	

def make_plots(rrv_x,wtmp,ch,cat):
        
        can             = []
        plots	        = []
	
	for i in range(len(POI)):
                c       = TCanvas(POI[i],POI[i],1)
                p       = rrv_x.frame()
                can.append(c)
                plots.append(p)
	for i in range(3):
		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		wtmp.data('SMdatahist').plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		wtmp.data('neg_datahist_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		normvalSM	= wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries()
		wtmp.pdf('aTGC_model').plotOn(plots[i],\
					      RooFit.LineColor(kBlack),\
					      RooFit.Normalization(normvalSM, RooAbsReal.NumEvent))
		wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		normval		= wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries()

        	wtmp.pdf('aTGC_model').plotOn(plots[i],\
					      RooFit.LineColor(kBlue),\
					      RooFit.Normalization(normval, RooAbsReal.NumEvent))
		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		for j in range(15):
			wtmp.var(POI[i]).setVal(-par_max[POI[i]]/10 * j+2)
			normval	= wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries()
			wtmp.pdf('aTGC_model').plotOn(plots[i],\
						      RooFit.LineColor(kGray+1),\
						      RooFit.LineWidth(1),\
						      RooFit.Normalization(normval,RooAbsReal.NumEvent))
		wtmp.data('SMdatahist').plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		wtmp.data('neg_datahist_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.MarkerSize(0.75),RooFit.DataError(RooAbsData.SumW2))
		if ch == 'el':
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



def make_input(ch = 'el'):

	cat		= options.cat	
	nbins4fit	= 20
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
	a1		= RooRealVar('a_SM_%s_%s'%(cat,ch),'a_SM_%s_%s'%(ch,cat),-0.1,-2,0)
	cwww		= RooRealVar('cwww','cwww',0,-120,120);
	ccw		= RooRealVar('ccw','ccw',0,-200,200);
	cb		= RooRealVar('cb','cb',0,-600,600);
	cwww.setConstant(kTRUE);
	ccw.setConstant(kTRUE);
	cb.setConstant(kTRUE);
	getattr(wtmp,'import')(cwww);
	getattr(wtmp,'import')(ccw);
	getattr(wtmp,'import')(cb);


	#make and fill SM histogram, SM fit
	SMhist		= TH1F('SMhist','SMhist',nbins4fit,binlo,binhi)
	SMhist.Sumw2(kTRUE)
	##read workspace containing background pdfs
	fileInWs	= TFile.Open('Input/wwlvj_%s_HPW_workspace.root'%ch[:2])
	w		= fileInWs.Get('workspace4limit_')
	rrv_mass_lvj	= w.var('rrv_mass_lvj')
	rrv_mass_lvj.setRange(binlo,binhi)

	SMPdf		= RooExponential('SMPdf_%s_%s'%(cat,ch),'SMPdf_%s_%s'%(cat,ch),rrv_mass_lvj,a1)
	for i in range(treeATGC.GetEntries()):
		treeATGC.GetEntry(i)
		if treeATGC.c_wwwl == 0 and treeATGC.c_wl ==0 and treeATGC.c_bl == 0:
	    		SMhist.Fill(treeATGC.MWW,treeATGC.totEventWeight)
	SMdatahist	= RooDataHist('SMdatahist','SMdatahist',RooArgList(rrv_mass_lvj),SMhist)
	##actual fit to determine SM shape parameter
	fitresSM	= SMPdf.fitTo(SMdatahist, RooFit.SumW2Error(kTRUE))
	a1.setConstant(kTRUE)
	N_SM		= RooRealVar('N_SM_%s_%s'%(cat,ch),'N_SM_%s_%s'%(cat,ch),SMdatahist.sumEntries())
	N_SM.setConstant(kTRUE)

	getattr(wtmp,'import')(SMdatahist)
	getattr(wtmp,'import')(N_SM)


	#make and fill ATGC histograms, fit quadratic scales
	##do this for all 3 parameters
	for para in POI:
		hist4fit = TH1F('hist4fit_%s'%para,'hist4fit_%s'%para,3,-1.5*par_max[para],1.5*par_max[para])
		pos_hist	= TH1F('c_pos_hist','c_pos_hist',nbins4fit,binlo,binhi)
		pos_hist.Sumw2(kTRUE)
		neg_hist	= TH1F('c_neg_hist','c_neg_hist',nbins4fit,binlo,binhi)
		neg_hist.Sumw2(kTRUE)

		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (para == 'cwww' and treeATGC.c_wwwl == 12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
			or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 20 and treeATGC.c_bl == 0)\
			or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == 60):
		      		pos_hist.Fill(treeATGC.MWW,treeATGC.totEventWeight)
		    	if (para == 'cwww' and treeATGC.c_wwwl == -12 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0)\
			or (para == 'ccw' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == -20 and treeATGC.c_bl == 0)\
			or (para == 'cb' and treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == -60):
		      		neg_hist.Fill(treeATGC.MWW,treeATGC.totEventWeight)

		pos_datahist	= RooDataHist('pos_datahist_%s'%para,'pos_datahist_%s'%para,RooArgList(rrv_mass_lvj),pos_hist)
		neg_datahist	= RooDataHist('neg_datahist_%s'%para,'neg_datahist_%s'%para,RooArgList(rrv_mass_lvj),neg_hist)
#
		hist4fit.SetBinContent(1,neg_datahist.sumEntries()/N_SM.getVal())
		hist4fit.SetBinContent(2,1)
		hist4fit.SetBinContent(3,pos_datahist.sumEntries()/N_SM.getVal())
		#fit parabel
		#gROOT.SetBatch(kTRUE)
		cc1 = TCanvas()
		cc1.cd()
		hist4fit.Fit('pol2')
		print str(pos_datahist.sumEntries()) +"/"+str(SMdatahist.sumEntries())+"/"+ str(neg_datahist.sumEntries())
		#raw_input("<")
		#gROOT.SetBatch(kFALSE)
		fitfunc		= hist4fit.GetFunction('pol2')
		par0		= RooRealVar('par0_%s_%s_%s'%(para,cat,ch),'par0_%s_%s_%s'%(para,cat,ch),fitfunc.GetParameter(0)); 		par0.setConstant(kTRUE);
		par1		= RooRealVar('par1_%s_%s_%s'%(para,cat,ch),'par1_%s_%s_%s'%(para,cat,ch),fitfunc.GetParameter(1)); 		par1.setConstant(kTRUE);
		par2		= RooRealVar('par2_%s_%s_%s'%(para,cat,ch),'par2_%s_%s_%s'%(para,cat,ch),fitfunc.GetParameter(2)); 		par2.setConstant(kTRUE);
#
		N_quad		= RooRealVar('N_quad_%s_%s_%s'%(para,cat,ch),'N_quad_%s_%s_%s'%(para,cat,ch), ((pos_datahist.sumEntries()+neg_datahist.sumEntries())/2)-N_SM.getVal() )
		N_lin		= RooRealVar('N_lin_%s_%s_%s'%(para,cat,ch),'N_lin_%s_%s_%s'%(para,cat,ch), (pos_datahist.sumEntries()-neg_datahist.sumEntries())/2 )

		#scaleshape is the relative change to SM		
		scaleshape	= RooFormulaVar('scaleshape_%s_%s_%s'%(para,cat,ch),'scaleshape_%s_%s_%s'%(para,cat,ch),\
						'(@0+@1*@3+@2*@3**2)-1',\
						RooArgList(par0,par1,par2,wtmp.var(para)))			

		a2		= RooRealVar('a_quad_%s_%s_%s'%(para,cat,ch),'a_quad_%s_%s_%s'%(para,cat,ch),-0.001,-1,0.1)
		a2.setConstant(kTRUE)
		cPdf_quad	= RooExponential('Pdf_quad_%s_%s_%s'%(para,cat,ch),'Pdf_quad_%s_%s_%s'%(para,cat,ch),rrv_mass_lvj,a2)
		a3		= RooRealVar('a_lin_%s_%s_%s'%(para,cat,ch),'a_lin_%s_%s_%s'%(para,cat,ch),-0.001,-0.1,0.1)
		a3.setConstant(kTRUE)
		cPdf_lin	= RooExponential('Pdf_lin_%s_%s_%s'%(para,cat,ch),'Pdf_lin_%s_%s_%s'%(para,cat,ch),rrv_mass_lvj,a3)

		getattr(wtmp,'import')(cPdf_quad)
		getattr(wtmp,'import')(cPdf_lin)
		getattr(wtmp,'import')(pos_datahist)
		getattr(wtmp,'import')(neg_datahist)
		getattr(wtmp,'import')(N_quad)
		getattr(wtmp,'import')(N_lin)
		getattr(wtmp,'import')(scaleshape)

		wtmp.Print()
	tmpp='''
	c1 = TCanvas()
	c1.cd()
	wtmp.var("cwww").setRange(-15,15)
	p1 = wtmp.var("cwww").frame()
	wtmp.function("scaleshape_cwww_%s_%s"%(cat,ch)).plotOn(p1,RooFit.LineColor(kBlue))
	p1.Draw()
	c2 = TCanvas()
	c2.cd()
	wtmp.var("ccw").setRange(-30,30)
	p2 = wtmp.var("ccw").frame()
	wtmp.function("scaleshape_ccw_%s_%s"%(cat,ch)).plotOn(p2,RooFit.LineColor(kRed))
	p2.Draw()
	c3 = TCanvas()
	c3.cd()
	wtmp.var("cb").setRange(-90,90)
	p3 = wtmp.var("cb").frame()
	wtmp.function("scaleshape_cb_%s_%s"%(cat,ch)).plotOn(p3,RooFit.LineColor(kGreen))
	p3.Draw()
	raw_input("SS")'''


	#make model
	if options.lin:
		paralist	= RooArgList(N_SM)
		paralist.add(RooArgList(wtmp.function('N_quad_%s_%s_%s'%(POI[0],cat,ch)),wtmp.function('N_lin_%s_%s_%s'%(POI[0],cat,ch)),wtmp.var('cwww'),\
					wtmp.function('N_quad_%s_%s_%s'%(POI[1],cat,ch)),wtmp.function('N_lin_%s_%s_%s'%(POI[1],cat,ch)),wtmp.var('ccw'),\
					wtmp.function('N_quad_%s_%s_%s'%(POI[2],cat,ch)),wtmp.function('N_lin_%s_%s_%s'%(POI[2],cat,ch)),wtmp.var('cb')))
		Pdf_norm	= RooFormulaVar( 'Pdf_norm_%s_%s'%(cat,ch),'Pdf_norm_%s_%s'%(cat,ch),'@0+@1*(@3/12)**2+@2*(@3/12)+@4*(@6/20)**2+@5*(@6/20)+@7*(@9/60)**2+@8*(@9/60)',paralist )
		paralistN1	= RooArgList(Pdf_norm,N_SM) 
		paralistN2	= RooArgList(Pdf_norm,paralist.at(1),paralist.at(3))
		paralistN3	= RooArgList(Pdf_norm,paralist.at(2),paralist.at(3)) 
		paralistN4	= RooArgList(Pdf_norm,paralist.at(4),paralist.at(6))
		paralistN5	= RooArgList(Pdf_norm,paralist.at(5),paralist.at(6))
		paralistN6	= RooArgList(Pdf_norm,paralist.at(7),paralist.at(9))
		paralistN7	= RooArgList(Pdf_norm,paralist.at(8),paralist.at(9))

		N1		= RooFormulaVar( 'N1_%s_%s'%(cat,ch),'N1_%s_%s'%(cat,ch),'@1/@0',paralistN1 )
		N2		= RooFormulaVar( 'N2_%s_%s'%(cat,ch),'N2_%s_%s'%(cat,ch),'(@1*(@2/12)**2)/@0',paralistN2 )
		N3		= RooFormulaVar( 'N3_%s_%s'%(cat,ch),'N3_%s_%s'%(cat,ch),'(@1*(@2/12))/@0',paralistN3 )
		N4		= RooFormulaVar( 'N4_%s_%s'%(cat,ch),'N4_%s_%s'%(cat,ch),'(@1*(@2/20)**2)/@0',paralistN4 )
		N5		= RooFormulaVar( 'N5_%s_%s'%(cat,ch),'N5_%s_%s'%(cat,ch),'(@1*(@2/20))/@0',paralistN5 )
		N6		= RooFormulaVar( 'N6_%s_%s'%(cat,ch),'N6_%s_%s'%(cat,ch),'(@1*(@2/60)**2)/@0',paralistN6 )
		N7		= RooFormulaVar( 'N7_%s_%s'%(cat,ch),'N7_%s_%s'%(cat,ch),'(@1*(@2/60))/@0',paralistN7 )
	

		N_list		= RooArgList(N1,N2,N3,N4,N5,N6,N7)
		Pdf_list	= RooArgList(SMPdf,
						wtmp.pdf('Pdf_quad_%s_%s_%s'%(POI[0],cat,ch)),wtmp.pdf('Pdf_lin_%s_%s_%s'%(POI[0],cat,ch)),\
						wtmp.pdf('Pdf_quad_%s_%s_%s'%(POI[1],cat,ch)),wtmp.pdf('Pdf_lin_%s_%s_%s'%(POI[1],cat,ch)),\
						wtmp.pdf('Pdf_quad_%s_%s_%s'%(POI[2],cat,ch)),wtmp.pdf('Pdf_lin_%s_%s_%s'%(POI[2],cat,ch)))
	else:

		paralist	= RooArgList(N_SM,wtmp.function('N_quad_%s_%s_%s'%(POI[0],cat,ch)),wtmp.var('cwww'),\
						wtmp.function('N_quad_%s_%s_%s'%(POI[1],cat,ch)),wtmp.var('ccw'),\
						wtmp.function('N_quad_%s_%s_%s'%(POI[2],cat,ch)),wtmp.var('cb'))
		N1		= RooFormulaVar( 'N1_%s_%s'%(cat,ch),'N1_%s_%s'%(cat,ch),'@0/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)
		N2		= RooFormulaVar( 'N2_%s_%s'%(cat,ch),'N2_%s_%s'%(cat,ch),'@1*(@2/12)**2/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)
		N3		= RooFormulaVar( 'N3_%s_%s'%(cat,ch),'N3_%s_%s'%(cat,ch),'@3*(@4/20)**2/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)
		N4		= RooFormulaVar( 'N4_%s_%s'%(cat,ch),'N4_%s_%s'%(cat,ch),'@5*(@6/60)**2/(@0+@1*(@2/12)**2+@3*(@4/20)**2+@5*(@6/60)**2)',paralist)

		N_list		= RooArgList(N1,N2,N3,N4)
		Pdf_list	= RooArgList(SMPdf,wtmp.pdf('Pdf_quad_%s_%s_%s'%(POI[0],cat,ch)),wtmp.pdf('Pdf_quad_%s_%s_%s'%(POI[1],cat,ch)),wtmp.pdf('Pdf_quad_%s_%s_%s'%(POI[2],cat,ch)))


	model		= RooAddPdf('aTGC_model','aTGC_model', Pdf_list, N_list)
	model.Print()

	scale_list	= RooArgList(wtmp.function('scaleshape_%s_%s_%s'%(POI[0],cat,ch)),\
					wtmp.function('scaleshape_%s_%s_%s'%(POI[1],cat,ch)),\
					wtmp.function('scaleshape_%s_%s_%s'%(POI[2],cat,ch)))
	normfactor_3d	= RooFormulaVar('normfactor_3d_%s_%s'%(cat,ch),'normfactor_3d_%s_%s'%(cat,ch),'1+@0+@1+@2',scale_list)



	getattr(WS,'import')(normfactor_3d)	
	getattr(wtmp,'import')(normfactor_3d)	
	wtmp.Print()
	
	#fit 3 pdfs
	fitresults	= []
	for i in range(3):
		wtmp.var(POI[0]).setVal(0); wtmp.var(POI[1]).setVal(0); wtmp.var(POI[2]).setVal(0);
		wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		wtmp.var('a_quad_%s_%s_%s'%(POI[i],cat,ch)).setConstant(kFALSE)
		wtmp.var('a_lin_%s_%s_%s'%(POI[i],cat,ch)).setConstant(kFALSE)
		fitres		= model.fitTo(wtmp.data('neg_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
		fitresults.append(fitres)
		wtmp.var('a_quad_%s_%s_%s'%(POI[i],cat,ch)).setConstant(kTRUE)
		wtmp.var('a_lin_%s_%s_%s'%(POI[i],cat,ch)).setConstant(kTRUE)


	for i in range(3):
		fitresults[i].Print()

	model.Print()	
	getattr(wtmp,'import')(model)

	if options.do_plots:
		make_plots(rrv_mass_lvj,wtmp,ch,cat)
	if options.do_plots2:
		#cross check
		canvas		= TCanvas(POI[0]+','+POI[1]+','+POI[2] , POI[0]+','+POI[1]+','+POI[2],1)
		p4		= rrv_mass_lvj.frame()
		hist_all3	= TH1F('hist_all3','hist_all3',nbins4fit,binlo,binhi)
		hist_all3.Sumw2(kTRUE)
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (treeATGC.c_wwwl == -12 and treeATGC.c_wl == -20 and treeATGC.c_bl == -60):
		      		hist_all3.Fill(treeATGC.MWW,treeATGC.totEventWeight)
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

	fileInWs	= TFile.Open('Input/wwlvj_%s_HP%s_workspace.root'%(ch[:2],cat[1]))
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
	output 	= TFile('%s/%s_%s.root'%(path,cat,ch),'recreate')

	data_obs.Write()
	WS.SetName('w'); WS.Print(); WS.Write();
	output.Close()
	print 'Write to file ' + output.GetName()

  
make_input('el')
make_input('mu')
