from ROOT import  *
from array import array

gSystem.Load("/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so")

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

def make_SMBKG_input(ch = "ele", POIs = ["cwww"], oldtrees = 0):
	do_plot		= 0
	nbins4fit	= 20
	binlo		= 1000
	binhi		= 3500
	#read bkg
	fileInWS	= TFile.Open("Input/wwlvj_BulkG_WW_lvjj_M800_%s_HPW_workspace.root"%(ch))
	w		= fileInWS.Get("workspace4limit_") 
	w.pdf("STop_xww_%s_HPW"%(ch[:2])).SetName("STop")
	w.pdf("TTbar_xww_%s_HPW"%(ch[:2])).SetName("TTbar")
	w.pdf("WJets_xww_%s_HPW"%(ch[:2])).SetName("WJets")
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
			if "cwww" in POIs:
				c_wwwl[0] = 12; c_wl[0]	= 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[0] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = -12; c_wl[0] = 0; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[1] * weight_part
				treeATGC.Fill()
			if "ccw" in POIs:
				c_wwwl[0] = 0; c_wl[0] = 20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[2] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = 0; c_wl[0] = -20; c_bl[0] = 0;
				weight[0] = treeInATGC.aTGCWeights[3] * weight_part
				treeATGC.Fill()
			if "cb" in POIs:
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = 60;
				weight[0] = treeInATGC.aTGCWeights[4] * weight_part
				treeATGC.Fill()
				c_wwwl[0] = 0; c_wl[0] = 0; c_bl[0] = -60;
				weight[0] = treeInATGC.aTGCWeights[5] * weight_part
				treeATGC.Fill()
		treeATGC.Write()
		print "--------> Write to file " + fileOutATGC.GetName()


	#prepare fit
	a1		= RooRealVar("a_SM_%s"%ch,"a_SM_%s"%ch,-0.1,-2,0)	
	cwww		= RooRealVar("cwww","cwww",0,-12,12); 		cwww.setConstant(kTRUE);
	ccw		= RooRealVar("ccw","ccw",0,-20,20);		ccw.setConstant(kTRUE);
	cb		= RooRealVar("cb","cb",0,-60,60);		cb.setConstant(kTRUE);
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
	exp_factors	= []
	par_max		= [12,20,60]
	PDFs		= []
	hists		= [SMdatahist]
	for para in POIs:
		if para == "cwww":
			hist4fit_cwww	= TH1F("hist4fit_cwww","hist4fit_cwww",3,-1.5*par_max[0],1.5*par_max[0])
			cwww_pos_hist	= TH1F("c_pos_hist_cwww","c_pos_hist_cwww",nbins4fit,binlo,binhi)
			cwww_neg_hist	= TH1F("c_neg_hist_cwww","c_neg_hist_cwww",nbins4fit,binlo,binhi)
			for i in range(treeATGC.GetEntries()):
				treeATGC.GetEntry(i)
			    	if (para == "cwww" and treeATGC.c_wwwl == par_max[0]):
			      		cwww_pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			    	if (para == "cwww" and treeATGC.c_wwwl == -par_max[0]):
			      		cwww_neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			cwww_pos_datahist	= RooDataHist("cwww_pos_datahist","cwww_pos_datahist",RooArgList(rrv_mass_lvj),cwww_pos_hist)
			cwww_neg_datahist	= RooDataHist("cwww_neg_datahist","cwww_neg_datahist",RooArgList(rrv_mass_lvj),cwww_neg_hist)
			hist4fit_cwww.SetBinContent(1,cwww_neg_datahist.sumEntries()/SMdatahist.sumEntries())
			hist4fit_cwww.SetBinContent(2,1)
			hist4fit_cwww.SetBinContent(3,cwww_pos_datahist.sumEntries()/SMdatahist.sumEntries())
			hist4fit_cwww.Fit("pol2")
			fitfunc_cwww		= hist4fit_cwww.GetFunction("pol2")
			par0_cwww		= RooRealVar("par0_cwww_%s"%ch,"par0_cwww_%s"%ch,fitfunc_cwww.GetParameter(0)); 		par0_cwww.setConstant(kTRUE);
			par1_cwww		= RooRealVar("par1_cwww_%s"%ch,"par1_cwww_%s"%ch,fitfunc_cwww.GetParameter(1)); 		par1_cwww.setConstant(kTRUE);
			par2_cwww		= RooRealVar("par2_cwww_%s"%ch,"par2_cwww_%s"%ch,fitfunc_cwww.GetParameter(2)); 		par2_cwww.setConstant(kTRUE);
			scaleshape_cwww		= RooFormulaVar("scaleshape_cwww","scaleshape_cwww","@0+@1*@3+@2*@3*@3",RooArgList(par0_cwww,par1_cwww,par2_cwww,cwww))
			scaleshapemax_cwww	= RooRealVar("scaleshapemax_cwww","scaleshapemax_cwww",hist4fit_cwww.GetMaximum());		scaleshapemax_cwww.setConstant(kTRUE);
			N1_cwww			= RooFormulaVar("N1_cwww","N1_cwww","abs((@0-1)/(@1-1))",RooArgList(scaleshape_cwww,scaleshapemax_cwww))
			N1_factors.append(N1_cwww)
			a2		= RooRealVar("a_%s_%s"%(para,ch),"a_%s_%s"%(para,ch),-0.001,-2,0)
			cwwwPdf		= RooExponential("cwwwPdf","cwwwPdf",rrv_mass_lvj,a2)
			PDFs.append(cwwwPdf)
			hists.extend([cwww_neg_datahist,cwww_pos_datahist])
			exp_factors.append(a2)
		if para == "ccw":
			par_max.append(20)
			hist4fit_ccw	= TH1F("hist4fit_ccw","hist4fit_ccw",3,-1.5*par_max[1],1.5*par_max[1])
			ccw_pos_hist	= TH1F("c_pos_hist_ccw","c_pos_hist_ccw",nbins4fit,binlo,binhi)
			ccw_neg_hist	= TH1F("c_neg_hist_ccw","c_neg_hist_ccw",nbins4fit,binlo,binhi)
			for i in range(treeATGC.GetEntries()):
				treeATGC.GetEntry(i)
			    	if (para == "ccw" and treeATGC.c_wl == par_max[1]):
			      		ccw_pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			    	if (para == "ccw" and treeATGC.c_wl == -par_max[1]):
			      		ccw_neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			ccw_pos_datahist	= RooDataHist("ccw_pos_datahist","ccw_pos_datahist",RooArgList(rrv_mass_lvj),ccw_pos_hist)
			ccw_neg_datahist	= RooDataHist("ccw_neg_datahist","ccw_neg_datahist",RooArgList(rrv_mass_lvj),ccw_neg_hist)
			hist4fit_ccw.SetBinContent(1,ccw_neg_datahist.sumEntries()/SMdatahist.sumEntries())
			hist4fit_ccw.SetBinContent(2,1)
			hist4fit_ccw.SetBinContent(3,ccw_pos_datahist.sumEntries()/SMdatahist.sumEntries())
			hist4fit_ccw.Fit("pol2")
			fitfunc_ccw		= hist4fit_ccw.GetFunction("pol2")
			par0_ccw		= RooRealVar("par0_ccw_%s"%ch,"par0_ccw_%s"%ch,fitfunc_ccw.GetParameter(0)); 		par0_ccw.setConstant(kTRUE);
			par1_ccw		= RooRealVar("par1_ccw_%s"%ch,"par1_ccw_%s"%ch,fitfunc_ccw.GetParameter(1)); 		par1_ccw.setConstant(kTRUE);
			par2_ccw		= RooRealVar("par2_ccw_%s"%ch,"par2_ccw_%s"%ch,fitfunc_ccw.GetParameter(2)); 		par2_ccw.setConstant(kTRUE);
			scaleshape_ccw		= RooFormulaVar("scaleshape_ccw","scaleshape_ccw","@0+@1*@3+@2*@3*@3",RooArgList(par0_ccw,par1_ccw,par2_ccw,ccw))
			scaleshapemax_ccw	= RooRealVar("scaleshapemax_ccw","scaleshapemax_ccw",hist4fit_ccw.GetMaximum());		scaleshapemax_ccw.setConstant(kTRUE);
			N1_ccw			= RooFormulaVar("N1_ccw","N1_ccw","abs((@0-1)/(@1-1))",RooArgList(scaleshape_ccw,scaleshapemax_ccw))
			N1_factors.append(N1_ccw)
			a3		= RooRealVar("a_%s_%s"%(para,ch),"a_%s_%s"%(para,ch),-0.001,-2,0)
			ccwPdf		= RooExponential("ccwPdf","ccwPdf",rrv_mass_lvj,a3)
			PDFs.append(ccwPdf)
			hists.extend([ccw_neg_datahist,ccw_pos_datahist])
			exp_factors.append(a3)
		if para == "cb":
			par_max.append(60)
			hist4fit_cb	= TH1F("hist4fit_cb","hist4fit_cb",3,-1.5*par_max[2],1.5*par_max[2])
			cb_pos_hist	= TH1F("c_pos_hist_cb","c_pos_hist_cb",nbins4fit,binlo,binhi)
			cb_neg_hist	= TH1F("c_neg_hist_cb","c_neg_hist_cb",nbins4fit,binlo,binhi)
			for i in range(treeATGC.GetEntries()):
				treeATGC.GetEntry(i)
			    	if (para == "cb" and treeATGC.c_bl == par_max[2]):
			      		cb_pos_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			    	if (para == "cb" and treeATGC.c_bl == -par_max[2]):
			      		cb_neg_hist.Fill(treeATGC.m_lvj,treeATGC.weight)
			cb_pos_datahist	= RooDataHist("cb_pos_datahist","cb_pos_datahist",RooArgList(rrv_mass_lvj),cb_pos_hist)
			cb_neg_datahist	= RooDataHist("cb_neg_datahist","cb_neg_datahist",RooArgList(rrv_mass_lvj),cb_neg_hist)
			hist4fit_cb.SetBinContent(1,cb_neg_datahist.sumEntries()/SMdatahist.sumEntries())
			hist4fit_cb.SetBinContent(2,1)
			hist4fit_cb.SetBinContent(3,cb_pos_datahist.sumEntries()/SMdatahist.sumEntries())
			hist4fit_cb.Fit("pol2")
			fitfunc_cb		= hist4fit_cb.GetFunction("pol2")
			par0_cb			= RooRealVar("par0_cb_%s"%ch,"par0_cb_%s"%ch,fitfunc_cb.GetParameter(0)); 		par0_cb.setConstant(kTRUE);
			par1_cb			= RooRealVar("par1_cb_%s"%ch,"par1_cb_%s"%ch,fitfunc_cb.GetParameter(1)); 		par1_cb.setConstant(kTRUE);
			par2_cb			= RooRealVar("par2_cb_%s"%ch,"par2_cb_%s"%ch,fitfunc_cb.GetParameter(2)); 		par2_cb.setConstant(kTRUE);
			scaleshape_cb		= RooFormulaVar("scaleshape_cb","scaleshape_cb","@0+@1*@3+@2*@3*@3",RooArgList(par0_cb,par1_cb,par2_cb,cb))
			scaleshapemax_cb	= RooRealVar("scaleshapemax_cb","scaleshapemax_cb",hist4fit_cb.GetMaximum());		scaleshapemax_cb.setConstant(kTRUE);
			N1_cb			= RooFormulaVar("N1_cb","N1_cb","abs((@0-1)/(@1-1))",RooArgList(scaleshape_cb,scaleshapemax_cb))
			N1_factors.append(N1_cb)
			a4		= RooRealVar("a_%s_%s"%(para,ch),"a_%s_%s"%(para,ch),-0.001,-2,0)
			cbPdf		= RooExponential("cbPdf","cbPdf",rrv_mass_lvj,a4)
			PDFs.append(cbPdf)
			hists.extend([cb_neg_datahist,cb_pos_datahist])
			exp_factors.append(a4)
	#make and fit model
	if len(PDFs) == 1:
		aTGC_pdf 	= PDFs[0]
		model 		= RooAddPdf("aTGC-model","aTGC-model",SMPdf,aTGC_pdf,N1_factors[0])
		fitres		= model.fitTo(hists[2],RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE));	exp_factors[0].setConstant(kTRUE);
		fitres.Print()
	if len(PDFs) == 2:
		if do_plot:
			can1		= TCanvas("canvas1",POIs[0],1)
			can2		= TCanvas("canvas2",POIs[1],1)
			p1 		= rrv_mass_lvj.frame()
			p2 		= rrv_mass_lvj.frame()
		model		= RooAddPdf("aTGC_model_%s"%ch,"aTGC_model_%s"%ch,RooArgList(PDFs[0],PDFs[1],SMPdf),RooArgList(N1_factors[0],N1_factors[1]))
		POI_list	= [N1_factors[0].getParameter("scaleshape_%s"%POIs[0]).getParameter("%s"%POIs[0]),N1_factors[1].getParameter("scaleshape_%s"%POIs[1]).getParameter("%s"%POIs[1])]
		#fit first pdf
		POI_list[0].setVal(par_max[0])
		POI_list[1].setVal(0)
		exp_factors[1].setConstant(kTRUE)
		fitres1		= model.fitTo(hists[2],RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE));	exp_factors[0].setConstant(kTRUE);
		if do_plot:
			hists[2].plotOn(p1,RooFit.MarkerColor(par_max[0]+1))
			model.plotOn(p1,RooFit.LineColor(par_max[0]+1))
		exp_factors[1].setConstant(kFALSE)
		#fit second pdf
		POI_list[1].setVal(par_max[1])
		POI_list[0].setVal(0)
		fitres2		= model.fitTo(hists[4],RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE));	exp_factors[1].setConstant(kTRUE);
		fitres1.Print();	
		fitres2.Print();
		if do_plot:
			hists[4].plotOn(p2,RooFit.MarkerColor(par_max[1]+1))
			model.plotOn(p2,RooFit.LineColor(par_max[1]+1))
			#add SM to plot
			hists[0].plotOn(p1,RooFit.MarkerColor(kBlack))
			hists[0].plotOn(p2,RooFit.MarkerColor(kBlack))
			POI_list[0].setVal(0);		POI_list[1].setVal(0);
			model.plotOn(p1,RooFit.LineColor(kBlack))
			model.plotOn(p2,RooFit.LineColor(kBlack))
			for i in range(par_max[0]):
				POI_list[0].setVal(i);		
				POI_list[1].setVal(0)
				normval0	= SMdatahist.sumEntries()*N1_factors[0].getParameter("scaleshape_%s"%POIs[0]).getVal()
				model.plotOn(p1,RooFit.LineWidth(2),RooFit.LineColor(i+1),RooFit.Normalization(normval0,RooAbsReal.NumEvent))
				hists[0].plotOn(p1,RooFit.MarkerColor(kBlack))
				hists[2].plotOn(p1,RooFit.MarkerColor(par_max[0]+1))
			for i in range(par_max[1]):
				POI_list[1].setVal(i);		
				POI_list[0].setVal(0)
				normval1	= SMdatahist.sumEntries()*N1_factors[1].getParameter("scaleshape_%s"%POIs[1]).getVal()
				model.plotOn(p2,RooFit.LineWidth(2),RooFit.LineColor(i+1),RooFit.Normalization(normval1,RooAbsReal.NumEvent))
				hists[0].plotOn(p2,RooFit.MarkerColor(kBlack))
				hists[4].plotOn(p2,RooFit.MarkerColor(par_max[1]+1))

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
			raw_input("...")

	if len(PDFs) == 3:
		raise RuntimeError("3 parameters not done yet!")


	#import to ws
	getattr(w,"import")(model)

	path	="../../CombinedEWKAnalysis/CommonTools/data/anomalousCoupling"
	cwww 	= "_cwww" if "cwww" in POIs else ""
	ccw = "_ccw" if "ccw" in POIs else "" 
	cb = "_cb" if "cb" in POIs else ""
	output 	= TFile("%s/ch_%s%s%s%s.root"%(path,ch,cwww,ccw,cb),"recreate")

	data_obs.Write()
	w.SetName("workspace"); w.Print(); w.Write();
	output.Close()
	print "Write to file " + output.GetName()

  
make_SMBKG_input("ele",["ccw"],0)
make_SMBKG_input("mu",["ccw"],0)
