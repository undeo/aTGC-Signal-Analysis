from ROOT import *
import os
from optparse import OptionParser
import CMS_lumi
from array import array

gSystem.Load('%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so'%os.environ['CMSSW_BASE'])

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

gSystem.Load("Util_cxx.so")

from ROOT import draw_error_band

parser	= OptionParser()
parser.add_option('--cat', dest='cat', default='WW', help='category, WW or WZ, defines signal region')
parser.add_option('--ch', dest='ch', default='el', help='channel, el or mu')
parser.add_option('-P', '--POI', dest='POI', default='cwww')
(options,args) = parser.parse_args()

ch		= options.ch
cat		= options.cat
channel		= cat+'_'+ch
latex_par	= {'cwww' : 'c_{WWW}/#Lambda^{2}=12 TeV^{-2}', 'ccw' : '#frac{c_{W}}{#Lambda^{2}}=20 TeV^{-2}', 'cb' : 'c_{B}/#Lambda^{2}=60 TeV^{-2}'}

fileInWS	= TFile.Open('../Input/wwlvj_%s_HP%s_workspace.root'%(ch,cat[1]))
w		= fileInWS.Get('workspace4limit_')
#fileInWS.Close()
fileInsig	= TFile.Open('../../../CombinedEWKAnalysis/CommonTools/data/anomalousCoupling/%s_%s.root'%(cat,ch))
w_sig		= fileInsig.Get('w')
#fileInsig.Close()
#fileInSM	= TFile.Open('../docuplots/make_plots/%s.root'%channel)
fileInSM	= TFile.Open('../Input/wwlvj_%s_HP%s_workspace.root'%(ch,cat[1]))
w_SM		= fileInSM.Get('workspace4limit_')
#w_SM		= fileInSM.Get('wtmp')
fileInSM.Close()
fileInData	= TFile.Open('../../%s_ws.root'%channel)
w_data		= fileInData.Get('proc_%s'%channel)
fileInData.Close()
w_err		= RooWorkspace('w_err')

w.Print()
w_sig.Print()
w_SM.Print()
w_err.Print()


#getattr(w_err,'import')(w.pdf('VV_xww_%s_HP%s'%(ch,cat[1])))
getattr(w_err,'import')(w_sig.pdf('STop'))
getattr(w_err,'import')(w_sig.pdf('TTbar'))
getattr(w_err,'import')(w_sig.pdf('WJets'))
getattr(w_err,'import')(w.var('rrv_number_VV_xww_%s_mj_prefit_signal_region'%ch))
getattr(w_err,'import')(w.var('rrv_number_STop_xww_%s_mj_postfit_signal_region'%ch))
getattr(w_err,'import')(w.var('rrv_number_TTbar_xww_%s_mj_postfit_signal_region'%ch))
getattr(w_err,'import')(w.var('rrv_number_WJets0_xww_%s_mj_postfit_signal_region'%ch))
getattr(w_err,'import')(w_sig.function('normfactor_3d_%s'%channel))
getattr(w_err,'import')(w_sig.pdf('aTGC_model_%s_%s'%(cat,ch)))



rrv_x		= w_err.var('rrv_mass_lvj')
rrv_x.setVal(1500)
rrv_x.setRange(900,3500)
rrv_x.SetName("observable")
bins	= RooBinning(26,900,3500)
rrv_x.setBinning(bins)


data_obs2	= w.data("data_obs_xww_%s_HP%s"%(ch,cat[1]))
data_obs	= w_data.data('data_obs')



#data_obs	= w_data.data("data_obs")
#data_obs	= RooDataSet('data_obs','data_obs',fileInsig.Get("data_obs_tree"),RooArgSet(rrv_x))


#VV_pdf		= w_err.pdf('VV_xww_%s_HP%s'%(ch,cat[1]))
STop_pdf	= w_err.pdf('STop')
TTbar_pdf	= w_err.pdf('TTbar')
WJets_pdf	= w_err.pdf('WJets')
Sig_pdf		= w_err.pdf('aTGC_model_%s_%s'%(cat,ch))
Sig_pdf.SetName('Sig')
VV_norm		= w_err.var('rrv_number_VV_xww_%s_mj_prefit_signal_region'%ch)
STop_norm	= w_err.var('rrv_number_STop_xww_%s_mj_postfit_signal_region'%ch)
TTbar_norm	= w_err.var('rrv_number_TTbar_xww_%s_mj_postfit_signal_region'%ch)
WJets_norm	= w_err.var('rrv_number_WJets0_xww_%s_mj_postfit_signal_region'%ch)
allbkg_pdf	= RooAddPdf('allbkg','allbkg',RooArgList(Sig_pdf,STop_pdf,TTbar_pdf,WJets_pdf),RooArgList(VV_norm,STop_norm,TTbar_norm,WJets_norm))

VV_norm_val	= VV_norm.getVal()
STop_norm_val	= STop_norm.getVal()
TTbar_norm_val	= TTbar_norm.getVal()
WJets_norm_val	= WJets_norm.getVal()

allbkg_norm	= RooRealVar('allbkg_norm','allbkg_norm',VV_norm_val+STop_norm_val+TTbar_norm_val+WJets_norm_val)
#allbkg_norm.setError(TMath.Sqrt(VV_norm.getError()*VV_norm.getError()+STop_norm.getError()*STop_norm.getError()+TTbar_norm.getError()*TTbar_norm.getError()+WJets_norm.getError()*WJets_norm.getError()))
allbkg_norm_val	= allbkg_norm.getVal()
bkg_norm_error	= {'WW_el' : 0.107, 'WW_mu' : 0.091, 'WZ_el' : 0.1167, 'WZ_mu' : 0.101, 'WV_el' : 0.2, 'WV_mu' : 0.2}
allbkg_norm.setError(allbkg_norm_val*bkg_norm_error[channel])
allbkg_norm_err_val = allbkg_norm.getError()



w_err.var('cwww').setVal(0)
w_err.var('ccw').setVal(0)
w_err.var('cb').setVal(0)
Sig_norm	= RooRealVar('Sig_norm','Sig_norm',w_err.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*VV_norm_val)


allbkgsig_pdf	= RooAddPdf('allbkgsig','allbkgsig',RooArgList(allbkg_pdf,Sig_pdf),RooArgList(allbkg_norm,Sig_norm))



floatparams	= RooArgList(w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig0'%(ch,cat[1])),w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig1'%(ch,cat[1])),w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig2'%(ch,cat[1])),w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig3'%(ch,cat[1])),w_err.var('Deco_TTbar_xww_signal_region_%s_HP%s_mlvj_13TeV_eig0'%(ch,cat[1])),w_err.var('Deco_TTbar_xww_signal_region_%s_HP%s_mlvj_13TeV_eig1'%(ch,cat[1])))

c		= TCanvas('c','c',1)
c.cd()
p		= rrv_x.frame(900,3500)
c.SetLogy()
pad1=TPad("pad1","pad1",0.,0. ,1,0.30); #pad1 - pull
pad2=TPad("pad2","pad2",0.,0.3,1.,1. ); #pad0
pad2.SetRightMargin(0.1);
pad2.SetTopMargin(0.1);
pad2.SetBottomMargin(0.0001);
pad1.SetRightMargin(0.1)
pad1.SetTopMargin(0)
pad1.SetBottomMargin(0.4)   
pad1.Draw();
pad2.Draw();
pad2.cd()
pad2.SetLogy()

draw_error_band(allbkg_pdf, rrv_x.GetName(), allbkg_norm, floatparams , w_err, p, kBlack, "F", 3013)
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val,RooAbsReal.NumEvent), RooFit.Components('Sig,STop,TTbar,WJets'), RooFit.DrawOption('F'), RooFit.FillColor(kGreen+1), RooFit.LineColor(kBlack), RooFit.LineWidth(1),RooFit.Name('WJets'))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val,RooAbsReal.NumEvent), RooFit.Components('Sig,STop,TTbar'), RooFit.DrawOption('F'), RooFit.FillColor(kOrange), RooFit.LineColor(kBlack), RooFit.LineWidth(1),RooFit.Name('TTbar'))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val,RooAbsReal.NumEvent), RooFit.Components('Sig,STop'), RooFit.DrawOption('F'), RooFit.FillColor(kRed), RooFit.LineColor(kBlack), RooFit.LineWidth(1),RooFit.Name('VV'))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val,RooAbsReal.NumEvent), RooFit.Components('STop'), RooFit.DrawOption('F'), RooFit.FillColor(kBlue), RooFit.LineColor(kBlack), RooFit.LineWidth(1),RooFit.Name('STop'))

allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop,TTbar,WJets'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop,TTbar'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))


draw_error_band(allbkg_pdf, rrv_x.GetName(), allbkg_norm, floatparams , w_err, p, kBlack, "F", 3013)



if options.POI=='cwww':
	w_err.var('cwww').setVal(12)
	w_err.var('ccw').setVal(0)
	w_err.var('cb').setVal(0)
elif options.POI=='ccw':
	w_err.var('cwww').setVal(0)
	w_err.var('ccw').setVal(20)
	w_err.var('cb').setVal(0)
elif options.POI=='cb':
	w_err.var('cwww').setVal(0)
	w_err.var('ccw').setVal(0)
	w_err.var('cb').setVal(60)
Sig_norm.setVal(w_err.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*VV_norm_val-VV_norm_val)
norm4error	= RooRealVar('norm4error','norm4error',allbkg_norm_val + Sig_norm.getVal())
allbkgsig_pdf.plotOn(p,RooFit.Normalization(norm4error.getVal(),RooAbsReal.NumEvent), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.LineStyle(2), RooFit.LineColor(kViolet), RooFit.Name('sig'))#, RooFit.LineStyle(kDashed))
w_err.var("slope_nuis").setError(0.05)
#sig_floatparams	= RooArgList(w_err.var('a_SM_4fit_%s'%channel),w_err.var('a_quad_4fit_cwww_%s'%channel),w_err.var('a_quad_4fit_ccw_%s'%channel),w_err.var('a_lin_4fit_ccw_%s'%channel),w_err.var('a_quad_4fit_cb_%s'%channel),w_err.var('a_lin_4fit_cb_%s'%channel),w_err.var('slope_nuis'))
#norm4error.setError(norm4error.getVal()*0.14)
#draw_error_band(allbkgsig_pdf, rrv_x.GetName(), norm4error, sig_floatparams , w_err, p, kMagenta, "F", 3001)
#rrv_x.setRange('hi',900,3500)
w_sig.var('slope_nuis').setVal(1)
rrv_x.Print()

data_histo	= data_obs.binnedClone("data_obs_binned","data_obs_binned").createHistogram("data_histo",rrv_x)
data_plot	= RooHist(data_histo,100)
data_plot.SetMarkerStyle(20)
data_plot.SetMarkerSize(1)
alpha		= 1-0.6827
for iPoint in range(data_plot.GetN()):
	N = data_plot.GetY()[iPoint]
	if N==0 :
		L = 0
	else:
		L = ROOT.Math.gamma_quantile(alpha/2.,N,1.)
		U = ROOT.Math.gamma_quantile_c(alpha/2,N+1,1.)
		data_plot.SetPointEYlow(iPoint, N-L)
		data_plot.SetPointEYhigh(iPoint, U-N)
		data_plot.SetPointEXlow(iPoint,0)
		data_plot.SetPointEXhigh(iPoint,0)
p.addPlotable(data_plot,"PE")




#hpull = p.pullHist('data_histo__rrv_mass_lvj','allbkg_Norm[rrv_mass_lvj]_Comp[Sig,STop,TTbar,WJets]');
hpull = p.pullHist('data_histo__observable','allbkg_Norm[observable]_Comp[Sig,STop,TTbar,WJets]');
x = Double(0.); y = Double(0.) ;
for ipoint in range(0,hpull.GetN()):
  hpull.GetPoint(ipoint,x,y);
  if(y == 0):
   hpull.SetPoint(ipoint,x,10)
gt = TH1F("gt","gt",26,900,3500);
gt.SetMinimum(-3.999);
gt.SetMaximum(3.999);
gt.SetDirectory(0);
gt.SetStats(0);
gt.SetLineStyle(0);
gt.SetMarkerStyle(20);
gt.GetXaxis().SetTitle("M_{WV} (GeV)");
gt.GetXaxis().SetLabelFont(42);
gt.GetXaxis().SetLabelOffset(0.02);
gt.GetXaxis().SetLabelSize(0.15);
gt.GetXaxis().SetTitleSize(0.15);
gt.GetXaxis().SetTitleOffset(1.3);
gt.GetXaxis().SetTitleFont(42);
gt.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}");
gt.GetYaxis().CenterTitle(True);
gt.GetYaxis().SetNdivisions(205);
gt.GetYaxis().SetLabelFont(42);
gt.GetYaxis().SetLabelOffset(0.007);
gt.GetYaxis().SetLabelSize(0.15);
gt.GetYaxis().SetTitleSize(0.14);
gt.GetYaxis().SetTitleOffset(1);
gt.GetYaxis().SetTitleFont(42);
#gt.GetXaxis().SetNdivisions(505)
hpull.SetHistogram(gt)

leg	= TLegend(0.6,0.35,0.89,0.85)
if ch == 'el':
	leg.SetHeader('e#nu,%s-category'%cat)
	datalatex	= 'Data W#rightarrowe#nu'
else:
	leg.SetHeader('#mu#nu,%s-category'%cat)
	datalatex	= 'Data W#rightarrow#mu#nu'
leg.AddEntry(p.findObject(data_plot.GetName()),datalatex,'pe')
leg.AddEntry(p.findObject('sig'),'signal %s'%latex_par[options.POI],'l')
leg.AddEntry(p.findObject('WJets'),'W+jets','F')
leg.AddEntry(p.findObject('TTbar'),'t#bar{t}','F')
leg.AddEntry(p.findObject('VV'),'WW/WZ','F')
leg.AddEntry(p.findObject('STop'),'Single Top','F')
leg.AddEntry(p.getObject(0),'Background uncertainty','F')
leg.SetBorderSize(0)
leg.SetFillColor(0)
#leg.SetTextFont(43)


p.GetYaxis().SetRangeUser(0.07,5e3)
p.GetXaxis().SetRangeUser(900,3500)
p.GetYaxis().SetTitleSize(0.05)
p.GetYaxis().SetLabelSize(0.06)
p.GetYaxis().SetTitleOffset(0.8)
p.GetYaxis().SetTitle('Events / (100 GeV)')
p.SetTitle('')
p.Draw()
leg.Draw("SAME")
CMS_lumi.lumi_13TeV 	= '2.3 fb^{-1}'
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText	= 'preliminary'
CMS_lumi.CMS_lumi(pad2,4,11)
pad1.cd()
hpull.SetTitle('')
hpull.GetYaxis().SetTitleOffset(0.25)
hpull.Draw("AP")
medianLine = TLine(900,0.,3500,0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
medianLine.Draw()
hpull.Draw("Psame");
c.Draw()
c.Update()
c.SaveAs('../docuplots/%s_%s_%s.png'%(cat,ch,options.POI))
c.SaveAs('../docuplots/%s_%s_%s.pdf'%(cat,ch,options.POI))

c2=TCanvas('c2','c2',1)
c2.cd()
pad3=TPad("pad3","pad3",0.,0. ,1,0.30); #pad1 - pull
pad4=TPad("pad4","pad4",0.,0.3,1.,1. ); #pad0
pad4.SetRightMargin(0.1);
pad4.SetTopMargin(0.1);
pad4.SetBottomMargin(0.0001);
pad3.SetRightMargin(0.1)
pad3.SetTopMargin(0)
pad3.SetBottomMargin(0.4)   
pad3.Draw();
pad4.Draw();
pad4.cd()
if ch=='el':
	plothi = 120
elif ch=='mu':
	plothi = 175
p.GetYaxis().SetRangeUser(0,plothi)
p.GetXaxis().SetRangeUser(900,1500)
p.GetYaxis().SetLabelSize(0.04)
pad4.cd()
p.Draw()
leg.Draw("SAME")
CMS_lumi.CMS_lumi(pad4,4,11)
pad3.cd()
hpull.GetXaxis().SetRangeUser(900,1500)
hpull.Draw("AP")
medianLine2 = TLine(900,0.,1500,0); medianLine2.SetLineWidth(2); medianLine2.SetLineColor(kRed);
medianLine2.Draw()

hpull.Draw("Psame")
c2.Draw()
c2.Update()
c2.SaveAs('../docuplots/%s_%s_%s_lin.png'%(cat,ch,options.POI))
c2.SaveAs('../docuplots/%s_%s_%s_lin.pdf'%(cat,ch,options.POI))

print VV_norm_val
print STop_norm_val
print TTbar_norm_val
print WJets_norm_val
print allbkg_norm_val
print data_obs.sumEntries()
print data_obs2.sumEntries()
raw_input('.,')

