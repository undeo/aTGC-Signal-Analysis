from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import random
import os

gSystem.Load('%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so'%os.environ['CMSSW_BASE'])
from ROOT import RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf

#gSystem.Load('Util_cxx.so')
#from ROOT import draw_error_band



POI		= ['cwww','ccw','cb']
par_max 	= {'cwww' : 12, 'ccw' : 20, 'cb' : 60}#atgc points
par_titles 	= {'cwww' : '#frac{c_{WWW}}{#Lambda^{2}}', 'ccw' : '#frac{c_{W}}{#Lambda^{2}}', 'cb' : '#frac{c_{B}}{#Lambda^{2}}'}#latex titles 

parser	= OptionParser()
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='recreate input trees')
parser.add_option('-p', '--plots', action='store_true', dest='make_plots', default=False, help='make plots')
parser.add_option('--p2', action='store_true', dest='make_plots2', default=False, help='make parabel plot')
parser.add_option('--p3', action='store_true', dest='make_plots3', default=False, help='plot')
parser.add_option('--cat', dest='cat', default='WW', help='category, WW or WZ, defines signal region')
parser.add_option('--std', action='store_true', dest='std', default=False, help='standard mode')
parser.add_option('--lin2', action='store_true', dest='lin2', default=False, help='include SM-interference, two slopes')
parser.add_option('--lin1', action='store_true', dest='lin1', default=False, help='include SM-interference, one slope')
parser.add_option('--inter', action='store_true', dest='inter', default=False, help='include aTGC-interference')
parser.add_option('--linter', action='store_true', dest='linter', default=False, help='include SM- and aTGC-interference')
parser.add_option('--eband', action='store_true', dest='eband', default=False, help='draw error band')



(options,args) = parser.parse_args()

signal_category	= options.cat

if signal_category == 'WW':
	sigreg 	= 'lo'
	mj_lo	= 65.
	mj_hi	= 85.
elif signal_category == 'WZ':
	sigreg	= 'hi'
	mj_lo	= 85.
	mj_hi	= 105
else:
	raise RuntimeError('cateogry not supported!')



gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def make_ATGCtree(ch='el',binlo=900,binhi=3500):

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

#tmp='''
def draw_error_band(ws, pdf, rrv_x,rrv_number_events,fitres,fitparname,mplot,channel,par,number_point=100,number_errorband=2000):
	
	if len(fitres)>2:
		raise RuntimeError('max 2 fitresults supported!')

	rand		= TRandom3(1234)

	rrv_x.Print()
	x_min 		= rrv_x.getMin()
	x_max 		= rrv_x.getMax()
	delta_x	 	= (x_max-x_min)/number_point
	width_x		= mplot.getFitRangeBinW()
	n_events_mean	= rrv_number_events.getVal()
	n_events_error	= rrv_number_events.getError()
	fitpartmp	= []
	print fitparname


	rrv_x_set	= RooArgSet(rrv_x)

	syst	= []
	for i in range(number_errorband):
		tmp_graph 	= TGraph(number_point+1)
		for j in range(len(fitres)):
			par_tmp = fitres[j].randomizePars()
			for k in range(len(par_tmp)):
				if i==0:
					print fitparname[j][k] + ' -> ' + par_tmp[k].GetName()
				ws.var(fitparname[j][k]).setVal(par_tmp[k].getVal())
		#n_events_tmp	= rand.Gaus(n_events_mean,n_events_error)
		n_events_tmp	= n_events_mean
		for j in range(number_point+1):
			rrv_x.setVal(x_min+delta_x*j)
			tmp_graph.SetPoint(j,x_min+delta_x*j,n_events_mean*pdf.getVal(RooArgSet(rrv_x))*width_x)
		syst.append(tmp_graph)

	val		= [0]*number_errorband
	e_up		= TGraph(number_point+1)
	e_dn		= TGraph(number_point+1)
	e_up.SetName('error_up_%s_%s'%(channel,par))
	e_dn.SetName('error_dn_%s_%s'%(channel,par))

	for i in range(number_point+1):
		for j in range(number_errorband):
			val[j] = syst[j].GetY()[i]
		val_sort	= sorted(val)
		e_up.SetPoint(i,x_min+delta_x*i,val_sort[int(0.16*number_errorband)])
		e_dn.SetPoint(i,x_min+delta_x*i,val_sort[int(0.84*number_errorband)])
	e_up.SetLineWidth(2)
	e_up.SetLineColor(kMagenta)

	e_dn.SetLineWidth(2)
	e_dn.SetLineColor(kMagenta)

	mplot.addObject(e_up)
	mplot.addObject(e_dn)




def make_plots(rrv_x,wtmp,ch,cat,channel,fitres):
        
        can             = []
	can2		= []
        plots	        = []
	plots2		= []
	
	rrv_x.setVal(2000)

	for i in range(3):
                c       = TCanvas(POI[i]+'-',POI[i]+'-',1)
		c2	= TCanvas(POI[i]+'+',POI[i]+'+',1)
                p       = rrv_x.frame()
		p2	= rrv_x.frame()
                can.append(c)
		can2.append(c2)
                plots.append(p)
		plots2.append(p2)
	for i in range(3):
		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		wtmp.data('SMdatahist').plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		wtmp.data('neg_datahist_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		normval	= RooRealVar('normvalSM','normvalSM',wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())
		normval.setError(0)

		normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())
		wtmp.pdf('aTGC_model_%s_%s'%(cat,ch)).plotOn(plots[i],\
					      RooFit.LineColor(kBlack),\
					      RooFit.Normalization(normval.getVal(), RooAbsReal.NumEvent))

		wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())

        	normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())
		wtmp.pdf('aTGC_model_%s_%s'%(cat,ch)).plotOn(plots[i],\
					      RooFit.LineColor(kBlue),\
					      RooFit.Normalization(normval.getVal(), RooAbsReal.NumEvent))

		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		if options.eband:
        		normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())
			#draw_error_band(wtmp,normval,fitres[0],'a_SM_%s'%channel,plots[i],channel,POI[i])
			draw_error_band(wtmp,wtmp.pdf('aTGC_model_%s'%channel),wtmp.var('rrv_mass_lvj'),normval,[fitres[0]],[['a_SM_4fit_%s'%channel]],plots[i],channel,POI[i])
		wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		if options.eband:
	        	normval.setVal(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries())
			#draw_error_band(wtmp,wtmp.pdf('aTGC_model_%s_%s'%(cat,ch)),wtmp.var('rrv_mass_lvj'),normval,fitres[i+1],['a_quad_%s_%s'%(POI[i],channel)],plots[i],channel,POI[i])
			if options.lin1 or options.inter or options.std:
				float_pars = [['a_SM_4fit_%s'%channel],['a_quad_4fit_%s_%s'%(POI[i],channel)]]
			if options.lin2 or options.linter:
				float_pars = [['a_SM_4fit_%s'%channel],['a_lin_4fit_%s_%s'%(POI[i],channel),'a_quad_4fit_%s_%s'%(POI[i],channel)]]
			draw_error_band(wtmp,wtmp.pdf('aTGC_model_%s'%channel),wtmp.var('rrv_mass_lvj'),normval,[fitres[0],fitres[i+1]],float_pars,plots[i],channel,POI[i])


		tmp='''
		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		for j in range(15):
			wtmp.var(POI[i]).setVal(-par_max[POI[i]]/10 * j+2)
			normval.setVal(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries())
			wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots[i],\
						      RooFit.LineColor(kGray+1),\
						      RooFit.LineWidth(1),\
						      RooFit.Normalization(normval.getVal(),RooAbsReal.NumEvent))
		'''
		wtmp.data('SMdatahist').plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		wtmp.data('neg_datahist_%s'%POI[i]).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		
		if ch == 'el':
			plotmin = 1e-2
			plotmax = 50
			if signal_category == 'WZ':
				plotmin = 3e-4
		elif ch == 'mu':
			plotmin = 1e-2
			plotmax = 1e2
			if signal_category == 'WZ':
				plotmin = 1e-3
				plotmax = 50
		plots[i].GetYaxis().SetRangeUser(plotmin,plotmax)
		can[i].cd()
		can[i].SetLogy()	
		plots[i].SetTitle('')
		plots[i].Draw()
		can[i].Update()
		can[i].SaveAs('docuplots/%s_neg_%s.pdf'%(POI[i],channel))
		can[i].SaveAs('docuplots/%s_neg_%s.png'%(POI[i],channel))

		
		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		wtmp.data('SMdatahist').plotOn(plots2[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		wtmp.data('pos_datahist_%s'%POI[i]).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())

		normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())
		wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],\
					      RooFit.LineColor(kBlack),\
					      RooFit.Normalization(normval.getVal(), RooAbsReal.NumEvent))
		wtmp.var(POI[i]).setVal(par_max[POI[i]])
		normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())

		normval.setVal(wtmp.function('normfactor_3d_%s'%channel).getVal() * wtmp.data('SMdatahist').sumEntries())
        	wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],\
					      RooFit.LineColor(kBlue),\
					      RooFit.Normalization(normval.getVal(), RooAbsReal.NumEvent))
		tmp='''
		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		for j in range(15):
			wtmp.var(POI[i]).setVal(par_max[POI[i]]/10 * j+2)
			normval.setVal(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries())
			wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],\
						      RooFit.LineColor(kGray+1),\
						      RooFit.LineWidth(1),\
						      RooFit.Normalization(normval.getVal(),RooAbsReal.NumEvent))
		'''

		for j in range(3):
			wtmp.var(POI[j]).setVal(0)
		if options.eband:
        		normval.setVal(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries())
			#draw_error_band(wtmp,normval,fitres[0],'a_SM_%s'%channel,plots[i],channel,POI[i])
			draw_error_band(wtmp,wtmp.pdf('aTGC_model_%s'%channel),wtmp.var('rrv_mass_lvj'),normval,[fitres[0]],[['a_SM_4fit_%s'%channel]],plots[i],channel,POI[i])
		wtmp.var(POI[i]).setVal(par_max[POI[i]])
		if options.eband:
			#draw_error_band(wtmp,wtmp.pdf('aTGC_model_%s_%s'%(cat,ch)),wtmp.var('rrv_mass_lvj'),normval,fitres[i+1],['a_quad_%s_%s'%(POI[i],channel)],plots[i],channel,POI[i])
        		normval.setVal(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal() * wtmp.data('SMdatahist').sumEntries())
			if options.lin1 or options.inter or options.std:
				float_pars = [['a_SM_4fit_%s'%channel],['a_quad_4fit_%s_%s'%(POI[i],channel)]]
			if options.lin2 or options.linter:
				float_pars = [['a_SM_4fit_%s'%channel],['a_lin_4fit_%s_%s'%(POI[i],channel),'a_quad_4fit_%s_%s'%(POI[i],channel)]]
			draw_error_band(wtmp,wtmp.pdf('aTGC_model_%s'%channel),wtmp.var('rrv_mass_lvj'),normval,[fitres[0],fitres[i+1]],float_pars,plots[i],channel,POI[i])


		wtmp.data('SMdatahist').plotOn(plots2[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		wtmp.data('pos_datahist_%s'%POI[i]).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
		plots2[i].GetYaxis().SetRangeUser(plotmin,plotmax)
		can2[i].cd()
		can2[i].SetLogy()	
		plots2[i].SetTitle('')
		plots2[i].Draw()
		can2[i].Update()
		can[i].SaveAs('docuplots/%s_pos_%s.pdf'%(POI[i],channel))
		can[i].SaveAs('docuplots/%s_pos_%s.png'%(POI[i],channel))


	raw_input('plots plotted')



def make_input(ch = 'el',binlo=900,binhi=3500):

	if not options.lin1 and not options.lin2 and not options.inter and not options.linter and not options.std:
		raise RuntimeError('no options chosen!')

	cat		= signal_category
	channel		= cat+"_"+ch
	nbins4fit	= (binhi-binlo)/100
	WS		= RooWorkspace("WS")

	#read data
	fileInData	= TFile.Open('Input/treeEDBR_data_xww_%s.root'%(ch))
	tree_tmp	= fileInData.Get('tree')  


	#read or make ATGC and new tree
	if options.newtrees:
		make_ATGCtree(ch)
	fileInATGC	= TFile.Open('Output/ATGC-Tree_%s_%s.root'%(sigreg,ch))
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


	#make and fill SM histogram, SM fit
	SMhist		= TH1F('SMhist','SMhist',nbins4fit,binlo,binhi)
	SMhist.Sumw2(kTRUE)
	##read workspace containing background pdfs
	fileInWs	= TFile.Open('Input/wwlvj_%s_HP%s_workspace.root'%(ch[:2],cat[1]))
	w		= fileInWs.Get('workspace4limit_')
	rrv_mass_lvj	= w.var('rrv_mass_lvj')
	rrv_mass_lvj.SetTitle('M_{WV}')
	rrv_mass_lvj.setRange(binlo,binhi)

	SMPdf		= RooExponential('SMPdf_%s'%channel,'SMPdf_%s'%channel,rrv_mass_lvj,a1)
	for i in range(treeATGC.GetEntries()):
		treeATGC.GetEntry(i)
		if treeATGC.c_wwwl == 0 and treeATGC.c_wl == 0 and treeATGC.c_bl == 0:
	    		SMhist.Fill(treeATGC.MWW,treeATGC.totEventWeight)
	SMdatahist	= RooDataHist('SMdatahist','SMdatahist',RooArgList(rrv_mass_lvj),SMhist)
	##actual fit to determine SM shape parameter
	fitresSM	= SMPdf.fitTo(SMdatahist, RooFit.SumW2Error(kTRUE), RooFit.Save(kTRUE))
	fitresults.append(fitresSM)
	a1_4fit.setConstant(kTRUE)
	N_SM		= RooRealVar('N_SM_%s'%channel,'N_SM_%s'%channel,SMdatahist.sumEntries())
	N_SM.setConstant(kTRUE)

	getattr(wtmp,'import')(cwww);
	getattr(wtmp,'import')(ccw);
	getattr(wtmp,'import')(cb);
	getattr(wtmp,'import')(eps4cbWZ)
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
		getattr(wtmp,'import')(pos_datahist)
		getattr(wtmp,'import')(neg_datahist)
#
		hist4fit.SetBinContent(1,neg_datahist.sumEntries()/N_SM.getVal())
		hist4fit.SetBinContent(2,1)
		hist4fit.SetBinContent(3,pos_datahist.sumEntries()/N_SM.getVal())
		#fit parabel
		gROOT.SetBatch(True)
		cc1 = TCanvas()
		cc1.cd()
		hist4fit.Fit('pol2')
		hist4fit.GetXaxis().SetTitle(par_titles[para]+' (TeV^{-2})')
		hist4fit.GetYaxis().SetTitle('N_{events}^{SM+%s} / N_{events}^{SM}'%par_titles[para])
		hist4fit.GetYaxis().SetTitleSize(0.04)
		hist4fit.Draw()
		cc1.Update()
		cc1.SaveAs("tmp_plots/%s_%s_%s_%s.png"%(para,channel,binlo,binhi))
		cc1.SaveAs('docuplots/yields_%s_%s.pdf'%(para,channel))
		gROOT.SetBatch(False)
		print str(pos_datahist.sumEntries()) +"/"+str(SMdatahist.sumEntries())+"/"+ str(neg_datahist.sumEntries())
		cc1.Close()
		fitfunc		= hist4fit.GetFunction('pol2')
		par0		= RooRealVar('par0_%s_%s'%(para,channel),'par0_%s_%s'%(para,channel),fitfunc.GetParameter(0)); 		par0.setConstant(kTRUE);
		par1		= RooRealVar('par1_%s_%s'%(para,channel),'par1_%s_%s'%(para,channel),fitfunc.GetParameter(1)); 		par1.setConstant(kTRUE);
		par2		= RooRealVar('par2_%s_%s'%(para,channel),'par2_%s_%s'%(para,channel),fitfunc.GetParameter(2)); 		par2.setConstant(kTRUE);
#
		N_quad		= RooRealVar('N_quad_%s_%s'%(para,channel),'N_quad_%s_%s'%(para,channel), ((pos_datahist.sumEntries()+neg_datahist.sumEntries())/2)-N_SM.getVal() )
		N_lin		= RooRealVar('N_lin_%s_%s'%(para,channel),'N_lin_%s_%s'%(para,channel), (pos_datahist.sumEntries()-neg_datahist.sumEntries())/2 )

		#scaleshape is the relative change to SM		
		scaleshape	= RooFormulaVar('scaleshape_%s_%s'%(para,channel),'scaleshape_%s_%s'%(para,channel),\
						'(@0+@1*@3+@2*@3**2)-1',\
						RooArgList(par0,par1,par2,wtmp.var(para)))			
		
		a2_4fit		= RooRealVar('a_quad_4fit_%s_%s'%(para,channel),'a_quad_4fit_%s_%s'%(para,channel),-0.001,-0.01,0.)
		a3_4fit		= RooRealVar('a_lin_4fit_%s_%s'%(para,channel),'a_lin_4fit_%s_%s'%(para,channel),-0.001,-0.01,0.)
		a2_4fit.setConstant(kTRUE)
		a3_4fit.setConstant(kTRUE)
		##bigger uncertainty for cb in WZ-category
		if cat=='WZ' and para=='cb':
			a2		= RooFormulaVar('a_quad_nuis_%s_%s'%(para,channel),'a_quad_nuis_%s_%s'%(para,channel),'@0*@1',RooArgList(a2_4fit,eps4cbWZ))
			a3		= RooFormulaVar('a_lin_nuis_%s_%s'%(para,channel),'a_lin_nuis_%s_%s'%(para,channel),'@0*@1',RooArgList(a3_4fit,eps4cbWZ))
		else:
			a2		= RooFormulaVar('a_quad_nuis_%s_%s'%(para,channel),'a_quad_nuis_%s_%s'%(para,channel),'@0*@1',RooArgList(a2_4fit,eps))
			a3		= RooFormulaVar('a_lin_nuis_%s_%s'%(para,channel),'a_lin_nuis_%s_%s'%(para,channel),'@0*@1',RooArgList(a3_4fit,eps))
		cPdf_quad	= RooExponential('Pdf_quad_%s_%s'%(para,channel),'Pdf_quad_%s_%s'%(para,channel),rrv_mass_lvj,a2)
		cPdf_lin	= RooExponential('Pdf_lin_%s_%s'%(para,channel),'Pdf_lin_%s_%s'%(para,channel),rrv_mass_lvj,a3)
		if cat=='WZ' and para=='cb':
			a2.Print()
			a3.Print()
			cPdf_quad.Print()
			cPdf_lin.Print()


		getattr(wtmp,'import')(cPdf_quad,RooFit.RecycleConflictNodes())
		getattr(wtmp,'import')(cPdf_lin,RooFit.RecycleConflictNodes())
		getattr(wtmp,'import')(N_quad)
		getattr(wtmp,'import')(N_lin)
		getattr(wtmp,'import')(scaleshape)


		wtmp.Print()

	#make model
	paralist	= RooArgList(N_SM)

	#include SM-interference
	if options.lin2:
		paralist.add(RooArgList(wtmp.function('N_quad_%s_%s'%(POI[0],channel)),wtmp.function('N_lin_%s_%s_%s'%(POI[0],cat,ch)),wtmp.var('cwww'),\
					wtmp.function('N_quad_%s_%s'%(POI[1],channel)),wtmp.function('N_lin_%s_%s_%s'%(POI[1],cat,ch)),wtmp.var('ccw'),\
					wtmp.function('N_quad_%s_%s'%(POI[2],channel)),wtmp.function('N_lin_%s_%s_%s'%(POI[2],cat,ch)),wtmp.var('cb')))
		Pdf_norm	= RooFormulaVar( 'Pdf_norm_%s'%channel,'Pdf_norm_%s'%channel,\
							'@0+@1*(@3/12)**2+@2*(@3/12)+@4*(@6/20)**2+@5*(@6/20)+@7*(@9/60)**2+@8*(@9/60)',paralist )
		paralistN1	= RooArgList(Pdf_norm,N_SM) 
		paralistN2	= RooArgList(Pdf_norm,paralist.at(1),paralist.at(3))
		paralistN3	= RooArgList(Pdf_norm,paralist.at(2),paralist.at(3)) 
		paralistN4	= RooArgList(Pdf_norm,paralist.at(4),paralist.at(6))
		paralistN5	= RooArgList(Pdf_norm,paralist.at(5),paralist.at(6))
		paralistN6	= RooArgList(Pdf_norm,paralist.at(7),paralist.at(9))
		paralistN7	= RooArgList(Pdf_norm,paralist.at(8),paralist.at(9))

		N1		= RooFormulaVar( 'N1_%s'%channel,'N1_%s'%channel,'@1/@0',paralistN1 )
		N2		= RooFormulaVar( 'N2_%s'%channel,'N2_%s'%channel,'(@1*(@2/12)**2)/@0',paralistN2 )
		#N3		= RooFormulaVar( 'N3_%s'%channel,'N3_%s'%channel,'(@1*(@2/12))/@0',paralistN3 )
		N3		= RooRealVar( 'N3_%s'%channel,'N3_%s'%channel,0 )
		N4		= RooFormulaVar( 'N4_%s'%channel,'N4_%s'%channel,'(@1*(@2/20)**2)/@0',paralistN4 )
		N5		= RooFormulaVar( 'N5_%s'%channel,'N5_%s'%channel,'(@1*(@2/20))/@0',paralistN5 )
		N6		= RooFormulaVar( 'N6_%s'%channel,'N6_%s'%channel,'(@1*(@2/60)**2)/@0',paralistN6 )
		N7		= RooFormulaVar( 'N7_%s'%channel,'N7_%s'%channel,'(@1*(@2/60))/@0',paralistN7 )
		#N7		= RooRealVar( 'N7_%s'%channel,'N7_%s'%channel,0 )


		N_list		= RooArgList(N1,N2,N3,N4,N5,N6,N7)
		Pdf_list	= RooArgList(SMPdf,
						wtmp.pdf('Pdf_quad_%s_%s'%(POI[0],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[0],channel)),\
						wtmp.pdf('Pdf_quad_%s_%s'%(POI[1],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[1],channel)),\
						wtmp.pdf('Pdf_quad_%s_%s'%(POI[2],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[2],channel)))
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
			wtmp.var('a_quad_4fit_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			wtmp.var('a_lin_4fit_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			fitres		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
			fitresults.append(fitres)
			wtmp.var('a_quad_4fit_%s_%s'%(POI[i],channel)).setConstant(kTRUE)
			wtmp.var('a_lin_4fit_%s_%s'%(POI[i],channel)).setConstant(kTRUE)

		print len(fitresults)
		for i in range(4):
			fitresults[i].Print()
	
		model.Print()	
		getattr(wtmp,'import')(model)

	#include SM- and aTGC-interference
	if options.linter:
		#get ratio WW/WZ in WW- and WZ-category for each aTGC-parameter
		WWinWW,WZinWW,WWinWZ,WZinWZ = 0,0,0,0
		fileInWWinWW	= TFile.Open('Output/ATGC-Tree_WW_lo_%s.root'%ch)
		treeInWWinWW	= fileInWWinWW.Get('BasicTree')
		for i in range(treeInWWinWW.GetEntries()):
			treeInWWinWW.GetEntry(i)
			WWinWW		+= treeInWWinWW.totEventWeight
		fileInWWinWZ	= TFile.Open('Output/ATGC-Tree_WW_hi_%s.root'%ch)
		treeInWWinWZ	= fileInWWinWZ.Get('BasicTree')
		for i in range(treeInWWinWZ.GetEntries()):
			treeInWWinWZ.GetEntry(i)
			WWinWZ		+= treeInWWinWZ.totEventWeight
		fileInWZinWW	= TFile.Open('Output/ATGC-Tree_WZ_lo_%s.root'%ch)
		treeInWZinWW	= fileInWZinWW.Get('BasicTree')
		for i in range(treeInWZinWW.GetEntries()):
			treeInWZinWW.GetEntry(i)
			WZinWW		+= treeInWZinWW.totEventWeight
		fileInWZinWZ	= TFile.Open('Output/ATGC-Tree_WZ_hi_%s.root'%ch)
		treeInWZinWZ	= fileInWZinWZ.Get('BasicTree')
		for i in range(treeInWZinWZ.GetEntries()):
			treeInWZinWZ.GetEntry(i)
			WZinWZ		+= treeInWZinWZ.totEventWeight
		ratio_WW_reg	= (WWinWW)/(WWinWW+WZinWW)
		ratio_WZ_reg	= (WWinWZ)/(WWinWZ+WZinWZ)

		WWWZ_r		= {'WW' : ratio_WW_reg, 'WZ' : ratio_WZ_reg}
		print WWWZ_r


		##workaround for WW el channel -> neglected anyway
		fileInterWWel	= TFile.Open('Input/genlevel_WW_mu.root')
		w2WW		= fileInterWWel.Get('w2')
		fileInterWZ	= TFile.Open('Input/genlevel_WZ_%s.root'%ch)
		w2WZ		= fileInterWZ.Get('w2')
		a5WW		= w2WW.var('a5').getVal()
		a5WZ		= w2WZ.var('a5').getVal()
		a5val		= a5WZ
		a5_tmp		= RooRealVar('a_cwww_ccw_%s'%channel,'a_cwww_ccw_%s'%channel,w2WW.var('a5').getVal())
		##
		#read results of generator level study
		fileInterWW	= TFile.Open('Input/genlevel_WW_%s.root'%ch)
		w2WW		= fileInterWW.Get('w2')
		fileInterWZ	= TFile.Open('Input/genlevel_WZ_%s.root'%ch)
		w2WZ		= fileInterWZ.Get('w2')
		

		#define slopes and coefficients
		#a5WW		= w2WW.var('a5').getVal()
		#a5WZ		= w2WZ.var('a5').getVal()
		#a5val		= WWWZ_r[cat]*a5WW + (1-WWWZ_r[cat])*a5WZ
		#a5_tmp		= RooRealVar('a_cwww_ccw_%s'%channel,'a_cwww_ccw_%s'%channel,a5val)
		a6WW		= w2WW.var('a6').getVal()
		a6WZ		= w2WZ.var('a6').getVal()
		a6val		= WWWZ_r[cat]*a6WW + (1-WWWZ_r[cat])*a6WZ
		a6_tmp		= RooRealVar('a_cwww_cb_%s'%channel,'a_cwww_cb_%s'%channel,a6val)
		a7WW		= w2WW.var('a7').getVal()
		a7WZ		= w2WZ.var('a7').getVal()
		a7val		= WWWZ_r[cat]*a7WW + (1-WWWZ_r[cat])*a7WZ
		a7_tmp		= RooRealVar('a_ccw_cb_%s'%channel,'a_ccw_cb_%s'%channel,a7val)
		a5		= RooFormulaVar('a_cwww_ccw_nuis_%s'%channel,'a_cwww_ccw_nuis_%s'%channel,'@0*@1',RooArgList(a5_tmp,eps))
		a5_tmp.setConstant(kTRUE)
		a6_tmp.setConstant(kTRUE)
		a7_tmp.setConstant(kTRUE)
		if cat=='WZ':
			a6		= RooFormulaVar('a_cwww_cb_nuis_%s'%channel,'a_cwww_cb_nuis_%s'%channel,'@0*@1',RooArgList(a6_tmp,eps4cbWZ))
			a7		= RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,eps4cbWZ))
		else:
			a6		= RooFormulaVar('a_cwww_cb_nuis_%s'%channel,'a_cwww_cb_nuis_%s'%channel,'@0*@1',RooArgList(a6_tmp,eps))
			a7		= RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,eps))

		NSMWW		= w2WW.var('N_SM').getVal()
		NSMWZ		= w2WZ.var('N_SM').getVal()
		NSM		= WWWZ_r[cat]*NSMWW + (1-WWWZ_r[cat])*NSMWZ

		Pdf_cwww_ccw	= RooExponential('Pdf_cwww_ccw_%s'%channel,'Pdf_cwww_ccw_%s'%channel,rrv_mass_lvj,a5)
		Pdf_cwww_cb	= RooExponential('Pdf_cwww_cb_%s'%channel,'Pdf_cwww_cb_%s'%channel,rrv_mass_lvj,a6)
		Pdf_ccw_cb	= RooExponential('Pdf_ccw_cb_%s'%channel,'Pdf_ccw_cb_%s'%channel,rrv_mass_lvj,a7)

		#get factor to scale number of events to simulation level (ratio of N_events for all atgc-parameters negative)
		hist_all3	= TH1F('hist_all3','hist_all3',nbins4fit,binlo,binhi)
		hist_all3.Sumw2(kTRUE)
		for i in range(treeATGC.GetEntries()):
			treeATGC.GetEntry(i)
		    	if (treeATGC.c_wwwl == -12 and treeATGC.c_wl == -20 and treeATGC.c_bl == -60):
		      		hist_all3.Fill(treeATGC.MWW,treeATGC.totEventWeight)
		datahist_all3	= RooDataHist('datahist_all3','datahist_all3',RooArgList(rrv_mass_lvj),hist_all3)
		N_4norm		= WWWZ_r[cat]*w2WW.var('N_4norm').getVal() + (1-WWWZ_r[cat])*w2WZ.var('N_4norm').getVal()
		corr_factor	= datahist_all3.sumEntries() / N_4norm

		
		#define more coefficients
		N1220WW		= w2WW.var('N_cwww_ccw__12__20').getVal()
		N1260WW		= w2WW.var('N_cwww_cb__12__60').getVal()
		N2060WW		= w2WW.var('N_ccw_cb__20__60').getVal()
		N1220WZ		= w2WZ.var('N_cwww_ccw__12__20').getVal()
		N1260WZ		= w2WZ.var('N_cwww_cb__12__60').getVal()
		N2060WZ		= w2WZ.var('N_ccw_cb__20__60').getVal()
		#scale coefficients, take into account that there are no WZ events for cb
		N1220		= WWWZ_r[cat]*N1220WW + (1-WWWZ_r[cat])*N1220WZ
		N1260		= WWWZ_r[cat]*N1260WW + (1-WWWZ_r[cat])*N1260WZ
		N2060		= WWWZ_r[cat]*N2060WW + (1-WWWZ_r[cat])*N2060WZ
	
		NSMval= wtmp.data('SMdatahist').sumEntries()
		N12val 	= wtmp.data('pos_datahist_cwww').sumEntries()
		N12_val	= wtmp.data('neg_datahist_cwww').sumEntries()
		N20val	= wtmp.data('pos_datahist_ccw').sumEntries()
		N20_val	= wtmp.data('neg_datahist_ccw').sumEntries()
		N60val	= wtmp.data('pos_datahist_cb').sumEntries()
		N60_val	= wtmp.data('neg_datahist_cb').sumEntries()

		N12	= WWWZ_r[cat]*w2WW.var('N_cwww_12').getVal() + (1-WWWZ_r[cat])*w2WZ.var('N_cwww_12').getVal()
		N12_	= WWWZ_r[cat]*w2WW.var('N_cwww__12').getVal() + (1-WWWZ_r[cat])*w2WZ.var('N_cwww__12').getVal()
		N20	= WWWZ_r[cat]*w2WW.var('N_ccw_20').getVal() + (1-WWWZ_r[cat])*w2WZ.var('N_ccw_20').getVal()
		N20_	= WWWZ_r[cat]*w2WW.var('N_ccw__20').getVal() + (1-WWWZ_r[cat])*w2WZ.var('N_ccw__20').getVal()
		N60	= WWWZ_r[cat]*w2WW.var('N_cb_60').getVal() + (1-WWWZ_r[cat])*w2WZ.var('N_cb_60').getVal()
		N60_	= WWWZ_r[cat]*w2WW.var('N_cb__60').getVal() + (1-WWWZ_r[cat])*w2WZ.var('N_cb__60').getVal()

		print WWWZ_r

		corr_factor_SM 		= NSMval/NSM
		corr_factor_cwww	= N12val/N12
		corr_factor_ccw		= N20val/N20
		corr_factor_cb		= N60val/N60

		#N1220 	= N1220*corr_factor
		#N1260 	= N1260*corr_factor
		#N2060 	= N2060*corr_factor
		print str(N1220) + ',' + str(NSM) + ',' + str(N12_) + ',' + str(N20)
		print str(N12val/N12) + ' , ' + str(N12_val/N12_)
		print str(N20val/N20) + ' , ' + str(N20_val/N20_)
		print str(N60val/N60) + ' , ' + str(N60_val/N60_)

	
		cf_cwww_cb	= (corr_factor_cwww+corr_factor_cb)/2
		cf_cwww_ccw	= (corr_factor_cwww+corr_factor_ccw)/2
		cf_ccw_cb	= (corr_factor_ccw+corr_factor_cb)/2
		#cf		= (corr_factor_cwww+corr_factor_ccw+corr_factor_cb)/3
		cf 		= corr_factor
		#N_cwww_ccw	= RooRealVar('N_cwww_ccw_%s'%channel,'N_cwww_ccw_%s'%channel,\
		#				((N1220*cf+NSMval)-(N12_val+N20val)))
		#N_cwww_cb	= RooRealVar('N_cwww_cb_%s'%channel,'N_cwww_cb_%s'%channel,\
		#				((N1260*cf+NSMval)-(N12_val+N60val)))
		#N_ccw_cb	= RooRealVar('N_ccw_cb_%s'%channel,'N_ccw_cb_%s'%channel,\
		#				((N2060*cf+NSMval)-(N20_val+N60val)))

		##no cwww-ccw interference
		#N_cwww_ccw	= RooRealVar('N_cwww_ccw_%s'%channel,'N_cwww_ccw_%s'%channel,\
		#				(1-WWWZ_r[cat])*cf*((N1220+NSM)-(N12+N20_)))
		N_cwww_ccw	= RooRealVar('N_cwww_ccw_%s'%channel,'N_cwww_ccw_%s'%channel,0)
		N_cwww_cb	= RooRealVar('N_cwww_cb_%s'%channel,'N_cwww_cb_%s'%channel,\
						cf*((N1260+NSM)-(N12+N60_)))
		N_ccw_cb	= RooRealVar('N_ccw_cb_%s'%channel,'N_ccw_cb_%s'%channel,\
						cf*((N2060+NSM)-(N20+N60_)))

		#print str((N1220+NSM)-(N12_+N20))+'/'+str(N_cwww_ccw.getVal())+' : '+str(((N1220+NSM)-(N12_+N20))/(N_cwww_ccw.getVal()))
		print str((N1260+NSM)-(N12_+N60))+'/'+str(N_cwww_cb.getVal())+' : '+str(((N1260+NSM)-(N12_+N60))/(N_cwww_cb.getVal()))
		print str((N2060+NSM)-(N20_+N60))+'/'+str(N_ccw_cb.getVal())+' : '+str(((N2060+NSM)-(N20_+N60))/(N_ccw_cb.getVal()))
		print str(N_4norm)+'/'+str(N_4norm*cf)+'/'+str(datahist_all3.sumEntries())
		print str(N1220)+'/'+str(N1220*cf)
		print str(N1260)+'/'+str(N1260*cf)
		print str(N2060)+'/'+str(N2060*cf)


		paralist.add(RooArgList(wtmp.function('N_quad_%s_%s'%(POI[0],channel)),wtmp.function('N_lin_%s_%s'%(POI[0],channel)),wtmp.var('cwww'),\
					wtmp.function('N_quad_%s_%s'%(POI[1],channel)),wtmp.function('N_lin_%s_%s'%(POI[1],channel)),wtmp.var('ccw'),\
					wtmp.function('N_quad_%s_%s'%(POI[2],channel)),wtmp.function('N_lin_%s_%s'%(POI[2],channel)),wtmp.var('cb')))
		paralist.add(RooArgList(N_cwww_ccw,N_cwww_cb,N_ccw_cb))
		#neglect SM interference for cwww
		#cwww_s		= '+@1*(@3/12)**2+@2*(@3/12)'
		cwww_s		= '+@1*(@3/12)**2'
		ccw_s		= '+@4*(@6/20)**2+@5*(@6/20)'
		cb_s		= '+@7*(@9/60)**2+@8*(@9/60)'
		cwww_ccw_s	= ''#'+@10*(@3/12)*(@6/-20)'
		cwww_cb_s	= '+@11*(@3/12)*(@9/-60)'
		ccw_cb_s	= '+@12*(@6/20)*(@9/-60)'
		Pdf_norm	= RooFormulaVar( 'Pdf_norm_%s'%channel,'Pdf_norm_%s'%channel,\
							'@0'+cwww_s+ccw_s+cb_s+cwww_ccw_s+cwww_cb_s+ccw_cb_s,paralist)
		paralistN	= RooArgList()
		for i in range(13):
			paralistN.add(RooArgList(paralist.at(i)))
		paralistN.add(RooArgList(Pdf_norm))

		N1		= RooFormulaVar( 'N1_%s'%channel,'N1_%s'%channel,'@0/@13',paralistN )
		N2		= RooFormulaVar( 'N2_%s'%channel,'N2_%s'%channel,'(@1*(@3/12)**2)/@13',paralistN )
		#no SM-interference for c_WWW
		#N3		= RooFormulaVar( 'N3_%s'%channel,'N3_%s'%channel,'(@2*(@3/12))/@13',paralistN )
		N3		= RooRealVar( 'N3_%s'%channel,'N3_%s'%channel,0 )
		N4		= RooFormulaVar( 'N4_%s'%channel,'N4_%s'%channel,'(@4*(@6/20)**2)/@13',paralistN )
		N5		= RooFormulaVar( 'N5_%s'%channel,'N5_%s'%channel,'(@5*(@6/20))/@13',paralistN )
		N6		= RooFormulaVar( 'N6_%s'%channel,'N6_%s'%channel,'(@7*(@9/60)**2)/@13',paralistN )
		N7		= RooFormulaVar( 'N7_%s'%channel,'N7_%s'%channel,'(@8*(@9/60))/@13',paralistN )
		#no aTGC-interference for c_WWW/c_W
		#N8		= RooFormulaVar( 'N8_%s'%channel,'N8_%s'%channel,'(@10*(@3/12)*(@6/-20))/@13',paralistN )
		N8		= RooRealVar('N8_%s'%channel,'N8_%s'%channel,0)
		N9		= RooFormulaVar( 'N9_%s'%channel,'N9_%s'%channel,'(@11*(@3/12)*(@9/-60))/@13',paralistN )
		N10		= RooFormulaVar( 'N10_%s'%channel,'N10_%s'%channel,'(@12*(@6/20)*(@9/-60))/@13',paralistN )

		N_list		= RooArgList(N1,N2,N3,N4,N5,N6,N7)
		N_list.add(RooArgList(N8,N9,N10))
		Pdf_list	= RooArgList(SMPdf)
		Pdf_list.add(RooArgList(wtmp.pdf('Pdf_quad_%s_%s'%(POI[0],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[0],channel)),\
					wtmp.pdf('Pdf_quad_%s_%s'%(POI[1],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[1],channel)),\
					wtmp.pdf('Pdf_quad_%s_%s'%(POI[2],channel)),wtmp.pdf('Pdf_lin_%s_%s'%(POI[2],channel))))
		Pdf_list.add(RooArgList(Pdf_cwww_ccw,Pdf_cwww_cb,Pdf_ccw_cb))
		N_list.Print()
		Pdf_list.Print()
		model		= RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list)
		model.Print()

		scale_list	= RooArgList(wtmp.function('scaleshape_%s_%s'%(POI[0],channel)),\
						wtmp.function('scaleshape_%s_%s'%(POI[1],channel)),\
						wtmp.function('scaleshape_%s_%s'%(POI[2],channel)))
		normfactor_3d	= RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list)
		wtmp.Print()

		#fit 3 pdfs
		for i in range(3):
			for j in range(3):
				wtmp.var(POI[j]).setVal(0)
			wtmp.var(POI[i]).setVal(par_max[POI[i]])
			wtmp.var('a_quad_4fit_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			wtmp.var('a_lin_4fit_%s_%s'%(POI[i],channel)).setConstant(kFALSE)
			fitres		= model.fitTo(wtmp.data('pos_datahist_%s'%POI[i]),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
			fitresults.append(fitres)
			wtmp.var('a_quad_4fit_%s_%s'%(POI[i],channel)).setConstant(kTRUE)
			wtmp.var('a_lin_4fit_%s_%s'%(POI[i],channel)).setConstant(kTRUE)



		model.Print()
		getattr(WS,'import')(normfactor_3d)	
		getattr(wtmp,'import')(normfactor_3d)
		getattr(wtmp,'import')(model)
		
		#print coefficients to see contribution for all atgc-parameter positive
		for i in range(3):
			wtmp.var(POI[i]).setVal(par_max[POI[i]])
		for i in range(13):
			if i!=3 and i!=6 and i!=9:
				print paralist.at(i).getVal()

		#print fitresults
		for i in range(10):
			print 'N%s_%s'%(i+1,channel) + ' : ' + str(wtmp.function('N%s_%s'%(i+1,channel)).getVal())
		for i in range(4):
			fitresults[i].Print()
		a5.Print()
		a6.Print()
		a7.Print()

	#no interference
	if options.std:
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
		getattr(wtmp,'import')(model)

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
	output 	= TFile('%s/%s.root'%(path,channel),'recreate')

	data_obs = TH1F('data_obs','data_obs',nbins4fit,binlo,binhi)
	#data_obs	= TTree('data_obs','data_obs')
	#MWW_data	= array('f',[1])
	#data_obs.Branch('observable_%s'%channel,MWW_data,'observable_%s'%channel)
	for i in range(tree_tmp.GetEntries()):
	  	tree_tmp.GetEntry(i)
	    	#add cuts to data!
		#data still blinded!
		data_obs.Fill(random.random())
		#MWW_data[0] = random.random()
		#data_obs.Fill()
	data_obs.Write()

	WS.SetName('w')
	WS.Write();
	output.Close()
	print 'Write to file ' + output.GetName()


	#make plots
	if options.make_plots:
		make_plots(rrv_mass_lvj,wtmp,ch,cat,channel,fitresults)
	if options.make_plots2:
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
		#wtmp.data('SMdatahist').plotOn(p4,RooFit.MarkerColor(kBlack),RooFit.MarkerSize(0.75))

		datahist_all3.plotOn(p4,RooFit.LineColor(kBlue),RooFit.LineWidth(2),RooFit.DrawOption('E'))
		#for i in range(3):
		#		wtmp.var(POI[i]).setVal(0)
		#model.plotOn(p4,RooFit.LineColor(kBlack),RooFit.Normalization(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*wtmp.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))
		for i in range(3):
				wtmp.var(POI[i]).setVal(par_max[POI[i]])
		normval 	= RooRealVar('normval','normval',wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*wtmp.data('SMdatahist').sumEntries())
		model.plotOn(p4,RooFit.LineColor(kBlue),RooFit.Normalization(normval.getVal(),RooAbsReal.NumEvent))
		if options.eband:
			fitparas4plot	=['a_SM_4fit_%s'%channel,'a_quad_4fit_cwww_%s'%channel,'a_lin_4fit_cwww_%s'%channel,'a_quad_4fit_ccw_%s'%channel,'a_lin_4fit_ccw_%s'%channel,'a_quad_4fit_cb_%s'%channel,'a_lin_4fit_cb_%s'%channel]
			draw_error_band(wtmp,wtmp.pdf('aTGC_model_%s'%channel),wtmp.var('rrv_mass_lvj'),normval,fitresults,fitparas4plot,p4,channel,wtmp.var('cwww'))
			for i in range(3):
					wtmp.var(POI[i]).setVal(par_max[POI[i]])
			model.plotOn(p4,RooFit.LineColor(kRed),RooFit.LineStyle(kDashed),RooFit.Normalization(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*wtmp.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))

		#wtmp.data('SMdatahist').plotOn(p4,RooFit.MarkerColor(kBlack),RooFit.MarkerSize(0.75))
		datahist_all3.plotOn(p4,RooFit.LineColor(kBlue),RooFit.LineWidth(2),RooFit.DrawOption('E'))

		for i in range(3):
			wtmp.var(POI[i]).setVal(-par_max[POI[i]])
		normval = wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*wtmp.data('SMdatahist').sumEntries()
		print str(normval) + ' / ' + str(datahist_all3.sumEntries())
		print str((normval/datahist_all3.sumEntries() -1)*100) + ' %'

		for i in range(3):
			wtmp.var(POI[i]).setVal(0)
		for i in range(10):
			for j in range(3):
				wtmp.var(POI[j]).setVal((i+1)*par_max[POI[j]]/10)
			model.plotOn(p4,RooFit.LineColor(kMagenta),RooFit.LineWidth(1),RooFit.LineStyle(kDashed),RooFit.Normalization(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*wtmp.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))
			for j in range(3):
				wtmp.var(POI[j]).setVal((i+1)*(-par_max[POI[j]])/10)
			model.plotOn(p4,RooFit.LineColor(kOrange),RooFit.LineWidth(1),RooFit.LineStyle(kDashed),RooFit.Normalization(wtmp.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*wtmp.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))
	
		#eps.setVal(1.1)
		#model.plotOn(p4,RooFit.LineColor(kMagenta),RooFit.Normalization(normval,RooAbsReal.NumEvent))
		#eps.setVal(0.9)
		#model.plotOn(p4,RooFit.LineColor(kMagenta),RooFit.Normalization(normval,RooAbsReal.NumEvent))

		canvas.cd()
		plotmin = 0.5
		#if signal_category == 'WZ':
		#	plotmin = 1e-3
		p4.GetYaxis().SetRangeUser(plotmin,p4.GetMaximum()*3)
		canvas.SetLogy()
		p4.Draw()
		canvas.SetLogy
		canvas.Update()
		canvas.SaveAs('docuplots/atgc3_%s.pdf'%channel)
		canvas.SaveAs('docuplots/atgc3_%s.png'%channel)
			
		
		getattr(wtmp,'import')(datahist_all3)
		fileOut4plot = TFile.Open('docuplots/make_atgc3_plot/%s.root'%channel,'RECREATE')
		wtmp.Write()
		fileOut4plot.Close()

		wtmp.Print()
		raw_input('cross check plot plotted')
		canvas.Close()


#gROOT.SetBatch() 
make_input('el')
make_input('mu')
#i = 900
#binsize = 2600
#while i<3500:
#	make_input('el',i,i+binsize)
#	i = i + binsize

