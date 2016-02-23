from ROOT import  *
from array import array

gSystem.Load("/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so")

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

def make_SMBKG_input(ch = "ele", POIs = ["c_www"], oldtrees = 0):

	nbins4fit	= 20
	binlo		= 1000
	binhi		= 3500
	#read bkg
	fileInWS	= TFile.Open("Input/wwlvj_BulkG_WW_lvjj_M800_%s_HPW_workspace.root"%(ch))
	w		= fileInWS.Get("workspace4limit_") 
	singletop	= w.pdf("STop_xww_%s_HPW"%(ch[:2]))
	ttbar		= w.pdf("TTbar_xww_%s_HPW"%(ch[:2]))
	wjets		= w.pdf("WJets_xww_%s_HPW"%(ch[:2]))
	rrv_mass_lvj	= w.var("rrv_mass_lvj")
	rrv_mass_lvj.setRange(binlo,binhi) 
	#read data
	fileInData	= TFile.Open("Input/treeEDBR_data_xww_%s.root"%(ch))
	tree_tmp	= fileInData.Get("tree")  
	data_obs = TH1F("data_obs","data_obs",nbins4fit,binlo,binhi)
	for i in range(tree_tmp.GetEntries()):
	  	tree_tmp.GetEntry(i)
	    	data_obs.Fill(tree_tmp.m_lvj)
	#read or make ATGC and new tree
	if oldtrees:
		fileInATGC	= TFile.Open("Output/ATGC-Tree_%s"%ch)
		treeATGC	= fileInATGC.Get("BasicTree")
	else:
		fileInATGC	= TFile.Open("Input/WW-aTGC-%s.root"%ch)
		treeInATGC	= fileInATGC.Get("treeDumper/BasicTree")
		fileOutATGC	= TFile("Output/ATGC-Tree_%s.root"%ch, "recreate")
		treeInATGC.SetBranchStatus("*",0)
		treeATGC 	= treeInATGC.CloneTree(0)
		treeInATGC.SetBranchStatus("aTGCWeights",1); treeInATGC.SetBranchStatus("m_lvj",1); treeInATGC.SetBranchStatus("PUweight",1);
		weight		= array('f',[1]);	branch_weight	= treeATGC.Branch("weight",weight,"weight");
		m_lvj		= array('f',[1]);	branch_m_lvj	= treeATGC.Branch("m_lvj",m_lvj,"m_lvj");
		c_wwwl		= array('f',[1]);	branch_c_wwwl	= treeATGC.Branch("c_wwwl",c_wwwl,"c_wwwl")
		c_wl		= array('f',[1]);	branch_c_wl	= treeATGC.Branch("c_wl",c_wl,"c_wl")	
		c_bl		= array('f',[1]);	branch_c_bl	= treeATGC.Branch("c_bl",c_bl,"c_bl")

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
			if "c_www" in POIs:
				c_wwwl[0] = 12; c_wl[0]	= 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[0] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = -12; c_wl[0] = 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[1] * weight_part
				treeATGC.Fill()
			if "c_w" in POIs:
				c_wwwl[0] = 0; c_wl[0] = 20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[2] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = 0; c_wl[0] = -20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[3] * weight_part
				treeATGC.Fill()
			if "c_b" in POIs:
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = 60;
				weight[0] = treeInATGC.aTGCWeights[4] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = -60;
				weight[0] = treeInATGC.aTGCWeights[5] * weight_part
				treeATGC.Fill()
		treeATGC.Write()
		print "--------> Write to file " + fileOutATGC.GetName()


	#prepare fit
	a1		= RooRealVar("a1_%s"%ch,"a1_%s"%ch,-0.1,-2,0)	
	c_www		= RooRealVar("c_www","c_www",0,-12,12); 	c_www.setConstant(kTRUE);
	c_w		= RooRealVar("c_w","c_w",0,-20,20);		c_w.setConstant(kTRUE);
	c_b		= RooRealVar("c_b","c_b",0,-60,60);		c_b.setConstant(kTRUE);
	#make and fill SM histogram, SM fit
	SMhist		= TH1F("SMhist","SMhist",nbins4fit,binlo,binhi)
	SMPdf		= RooExponential("SMPdf","SMPdf",rrv_mass_lvj,a1)
	for i in range(treeATGC.GetEntries()):
		treeATGC.GetEntry(i)
		if treeATGC.c_wwwl == 0 and treeATGC.c_wl ==0 and treeATGC.c_bl == 0:
	    		SMhist.Fill(treeATGC.m_lvj,treeATGC.weight)
	SMdatahist	= RooDataHist("SMdatahist","SMdatahist",RooArgList(rrv_mass_lvj),SMhist)
	SMPdf.fitTo(SMdatahist, RooFit.SumW2Error(kTRUE));	a1.setConstant(kTRUE)
	#make and fill ATGC histograms
	N1_factors	= []
	PDFs		= []
	for para in POIs:
		if para == "cwww":
			par_max = 12
		if para == "cw":
			par_max = 20
		if para == "cb":
			par_max = 60
		hist4fit	= TH1F("hist4fit","hist4fit",3,-1.5*par_max,1.5*par_max)
		c_pos_hist	= TH1F("c_pos_hist","c_pos_hist",nbins4fit,binlo,binhi)
		c_neg_hist	= TH1F("c_neg_hist","c_neg_hist",nbins4fit,binlo,binhi)
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (para == "cwww" and treeATGC.c_wwwl == par_max) or (para == "cw" and treeATGC.c_wl == par_max) or (para == "cb" and treeATGC.c_bl == par_max):
		      		c_pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
		    	if (para == "cwww" and treeATGC.c_wwwl == -par_max) or (para == "cw" and treeATGC.c_wl == -par_max) or (para == "cb" and treeATGC.c_bl == -par_max):
		      		c_neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)

		c_pos_datahist	= RooDataHist("c_pos_datahist","c_pos_datahist",RooArgList(rrv_mass_lvj),c_pos_hist)
		c_neg_datahist	= RooDataHist("c_neg_datahist","c_neg_datahist",RooArgList(rrv_mass_lvj),c_neg_hist)
		hist4fit.SetBinContent(1,c_neg_datahist.sumEntries()/SMdatahist.sumEntries())
		hist4fit.SetBinContent(2,1)
		hist4fit.SetBinContent(3,c_pos_datahist.sumEntries()/SMdatahist.sumEntries())
		hist4fit.Fit("pol2")
		fitfunc		= hist4fit.GetFunction("pol2")
		if para == "cwww":
			par0_cwww		= RooRealVar("par0_cwww","par0_cwww",fitfunc.GetParameter(0)); 		par0_cwww.setConstant(kTRUE);
			par1_cwww		= RooRealVar("par1_cwww","par1_cwww",fitfunc.GetParameter(1)); 		par1_cwww.setConstant(kTRUE);
			par2_cwww		= RooRealVar("par2_cwww","par2_cwww",fitfunc.GetParameter(2)); 		par2_cwww.setConstant(kTRUE);
			scaleshape_cwww		= RooFormulaVar("scaleshape_cwww","scaleshape_cwww","@0+@1*@3+@2*@3*@3",RooArgList(par0_cwww,par1_cwww,par2_cwww,c_www))
			scaleshapemax_cwww	= RooRealVar("scaleshapemax_cwww","scaleshapemax_cwww",hist4fit.GetMaximum());		scaleshapemax_cwww.setConstant(kTRUE);
			N1_cwww			= RooFormulaVar("N1_cwww","N1_cwww","abs((@1-@0)/(@1-1))",RooArgList(scaleshape_cwww,scaleshapemax_cwww))
			N1_factors.append(N1_cwww)
			a2		= RooRealVar("a_%s_%s"%(para,ch),"a2_%s_%s"%(para,ch),-0.1,-2,0)
			cwwwPdf		= RooExponential("cwwwPdf","cwwwPdf",rrv_mass_lvj,a2)
			c_www.setVal(par_max)
			cwwwPdf.fitTo(c_pos_datahist);	a2.setConstant(kTRUE)
			PDFs.append(cwwwPdf)
		if para == "cw":
			par0_cw			= RooRealVar("par0_cw","par0_cw",fitfunc.GetParameter(0)); 		par0_cw.setConstant(kTRUE);
			par1_cw			= RooRealVar("par1_cw","par1_cw",fitfunc.GetParameter(1)); 		par1_cw.setConstant(kTRUE);
			par2_cw			= RooRealVar("par2_cw","par2_cw",fitfunc.GetParameter(2)); 		par2_cw.setConstant(kTRUE);
			scaleshape_cw		= RooFormulaVar("scaleshape_cw","scaleshape_cw","@0+@1*@3+@2*@3*@3",RooArgList(par0_cw,par1_cw,par2_cw,c_w))
			scaleshapemax_cw	= RooRealVar("scaleshapemax_cw","scaleshapemax_cw",hist4fit.GetMaximum());		scaleshapemax_cw.setConstant(kTRUE);
			N1_cw			= RooFormulaVar("N1_cw","N1_cw","abs((@1-@0)/(@1-1))",RooArgList(scaleshape_cw,scaleshapemax_cw))
			N1_factors.append(N1_cw)
			a3		= RooRealVar("a_%s_%s"%(para,ch),"a3_%s_%s"%(para,ch),-0.1,-2,0)
			cwPdf		= RooExponential("cwPdf","cwPdf",rrv_mass_lvj,a3)
			c_w.setVal(par_max)
			cwPdf.fitTo(c_pos_datahist);	a3.setConstant(kTRUE)
			PDFs.append(cwPdf)
		if para == "cb":
			par0_cb			= RooRealVar("par0_cb","par0_cb",fitfunc.GetParameter(0)); 		par0_cb.setConstant(kTRUE);
			par1_cb			= RooRealVar("par1_cb","par1_cb",fitfunc.GetParameter(1)); 		par1_cb.setConstant(kTRUE);
			par2_cb			= RooRealVar("par2_cb","par2_cb",fitfunc.GetParameter(2)); 		par2_cb.setConstant(kTRUE);
			scaleshape_cb		= RooFormulaVar("scaleshape_cb","scaleshape_cb","@0+@1*@3+@2*@3*@3",RooArgList(par0_cb,par1_cb,par2_cb,c_b))
			scaleshapemax_cb	= RooRealVar("scaleshapemax_cb","scaleshapemax_cb",hist4fit.GetMaximum());		scaleshapemax_cb.setConstant(kTRUE);
			N1_cb			= RooFormulaVar("N1_cb","N1_cb","abs((@1-@0)/(@1-1))",RooArgList(scaleshape_cb,scaleshapemax_cb))
			N1_factors.append(N1_cb)
			a4		= RooRealVar("a_%s_%s"%(para,ch),"a4_%s_%s"%(para,ch),-0.1,-2,0)
			cbPdf		= RooExponential("cbPdf","cbPdf",rrv_mass_lvj,a4)
			c_b.setVal(par_max)
			cbPdf.fitTo(c_pos_datahist);	a4.setConstant(kTRUE)
			PDFs.append(cbPdf)

	if len(PDFs) == 1:
		aTGC_pdf = PDFs[0]
	if len(PDFs) == 2:
		aTGC_pdf = RooAddPdf("aTGC-Pdf","aTGC-Pdf",RooArgList(PDFs[0],PDFs[1]))
	if len(PDFs) == 3:
		#aTGC_pdf = RooProdPdf("aTGC-Pdf","aTGC-Pdf",RooArgList(PDFs[0],PDFs[1],PDFs[2]))
		raise RuntimeError("Not implemented yet!")
	print aTGC_pdf
	aTGC_pdf.Print()
	exit(0)







	


	
	
	model 		= RooAddPdf("aTGC-model","aTGC-model",SMPdf,c12Pdf,N1)
	model.Print()

	#plot + canvas
	plot		= rrv_mass_lvj.frame()
	plot2		= c_www.frame()
	can1 		= TCanvas("can1","can1",1)
	can2 		= TCanvas("can2","can2",1)

	#model fit
	c_www.setVal(12)
	res = model.fitTo(cwww12datahist,RooFit.Save(true), RooFit.SumW2Error(kTRUE))
	a2.setConstant(kTRUE)
	cwww12datahist.plotOn(plot, RooFit.MarkerColor(13))
	model.plotOn(plot, RooFit.LineColor(13))
	SMdatahist.plotOn(plot, RooFit.MarkerColor(1))
	res.Print()

	#plot
	can1.cd()
	cwww12datahist.plotOn(plot, RooFit.MarkerColor(13))
	model.plotOn(plot, RooFit.LineColor(13))
	for i in range(13):	
		c_www.setVal(i)
		normval = SMdatahist.sumEntries()*scaleshape.getVal()
		model.plotOn(plot,RooFit.LineColor(i+1), RooFit.LineWidth(1),RooFit.Normalization(normval,RooAbsReal.NumEvent))
	plot.GetYaxis().SetRangeUser(0.03,50)
	can1.SetLogy()
	plot.Draw()
	can1.Update()

	can2.cd()
	scaleshape.plotOn(plot2)
	plot2.Draw()
	can2.Update()

	raw_input("~~")
	 
	#import to ws

	c_www.setVal(0)	
	getattr(w,"import")(model)
	getattr(w,"import")(scaleshape)
	singletop.SetName("STop"); getattr(w,"import")(singletop)
	ttbar.SetName("TTbar"); getattr(w,"import")(ttbar)
	wjets.SetName("WJets"); getattr(w,"import")(wjets) 

	path	="../../CombinedEWKAnalysis/CommonTools/data/anomalousCoupling"
	cwww 	= "_cwww" if "c_www" in POIs else ""; cw = "_cw" if "c_w" in POIs else ""; cb = "_cb" if "c_b" in POIs else "";
	output 	= TFile("%s/ch_%s%s%s%s.root"%(path,ch,cwww,cw,cb),"recreate")

	data_obs.Write()
	w.SetName("workspace"); w.Print(); w.Write();
	output.Close()
	print "Write to file " + output.GetName()

  
make_SMBKG_input("ele",["cwww","cw"],1)
#make_SMBKG_input("mu",["cwww","cw"],0) 
