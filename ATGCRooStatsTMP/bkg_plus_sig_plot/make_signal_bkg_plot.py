from ROOT import *
import os
from optparse import OptionParser
import CMS_lumi

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
latex_par	= {'cwww' : '#frac{c_{WWW}}{#Lambda^{2}}=12 TeV^{-2}', 'ccw' : '#frac{c_{W}}{#Lambda^{2}}=20 TeV^{-2}', 'cb' : '#frac{c_{B}}{#Lambda^{2}}=60 TeV^{-2}'}

fileInWS	= TFile.Open('../Input/wwlvj_%s_HP%s_workspace.root'%(ch,cat[1]))
w		= fileInWS.Get('workspace4limit_')
fileInWS.Close()
fileInsig	= TFile.Open('../../../CombinedEWKAnalysis/CommonTools/data/anomalousCoupling/%s_%s.root'%(cat,ch))
w_sig		= fileInsig.Get('w')
fileInsig.Close()
fileInSM	= TFile.Open('../docuplots/make_plots/%s.root'%channel)
w_SM		= fileInSM.Get('wtmp')
fileInSM.Close()
w_err		= RooWorkspace('w_err')

w.var('rrv_mass_lvj').setRange(900,3500)

getattr(w_err,'import')(w.pdf('VV_xww_%s_HP%s'%(ch,cat[1])))
#getattr(w_err,'import')(w.pdf('STop_xww_%s_HP%s'%(ch,cat[1])))
#getattr(w_err,'import')(w.pdf('TTbar_xww_%s_HP%s'%(ch,cat[1])))
#getattr(w_err,'import')(w.pdf('WJets_xww_%s_HP%s'%(ch,cat[1])))
getattr(w_err,'import')(w_sig.pdf('STop'))
getattr(w_err,'import')(w_sig.pdf('TTbar'))
getattr(w_err,'import')(w_sig.pdf('WJets'))
getattr(w_err,'import')(w.var('rate_VV_xww_for_unbin'))
getattr(w_err,'import')(w.var('rate_STop_xww_for_unbin'))
getattr(w_err,'import')(w.var('rate_TTbar_xww_for_unbin'))
getattr(w_err,'import')(w.var('rate_TTbar_xww_for_unbin'))
getattr(w_err,'import')(w.var('rate_WJets_xww_for_unbin'))
getattr(w_err,'import')(w_sig.function('normfactor_3d_%s'%channel))
getattr(w_err,'import')(w_sig.pdf('aTGC_model_%s_%s'%(cat,ch)))



rrv_x		= w_err.var('rrv_mass_lvj')
rrv_x.setRange(900,3500)

VV_pdf		= w_err.pdf('VV_xww_%s_HP%s'%(ch,cat[1]))
#STop_pdf	= w_err.pdf('STop_xww_%s_HP%s'%(ch,cat[1]))
#TTbar_pdf	= w_err.pdf('TTbar_xww_%s_HP%s'%(ch,cat[1]))
#WJets_pdf	= w_err.pdf('WJets_xww_%s_HP%s'%(ch,cat[1]))
STop_pdf	= w_err.pdf('STop')
TTbar_pdf	= w_err.pdf('TTbar')
WJets_pdf	= w_err.pdf('WJets')
Sig_pdf		= w_err.pdf('aTGC_model_%s_%s'%(cat,ch))
Sig_pdf.SetName('Sig')
VV_norm		= w_err.var('rate_VV_xww_for_unbin')
STop_norm	= w_err.var('rate_STop_xww_for_unbin')
TTbar_norm	= w_err.var('rate_TTbar_xww_for_unbin')
WJets_norm	= w_err.var('rate_WJets_xww_for_unbin')
#allbkg_pdf	= RooAddPdf('allbkg','allbkg',RooArgList(VV_pdf,STop_pdf,TTbar_pdf,WJets_pdf),RooArgList(VV_norm,STop_norm,TTbar_norm,WJets_norm))
allbkg_pdf	= RooAddPdf('allbkg','allbkg',RooArgList(Sig_pdf,STop_pdf,TTbar_pdf,WJets_pdf),RooArgList(VV_norm,STop_norm,TTbar_norm,WJets_norm))



VV_norm_val	= VV_norm.getVal()
STop_norm_val	= STop_norm.getVal()
TTbar_norm_val	= TTbar_norm.getVal()
WJets_norm_val	= WJets_norm.getVal()
#allbkg_norm	= RooRealVar('allbkg_norm','allbkg_norm',VV_norm_val+STop_norm_val+TTbar_norm_val+WJets_norm_val)
#allbkg_norm.setError(TMath.Sqrt(VV_norm.getError()*VV_norm.getError()+STop_norm.getError()*STop_norm.getError()+TTbar_norm.getError()*TTbar_norm.getError()+WJets_norm.getError()*WJets_norm.getError()))
allbkg_norm	= RooRealVar('allbkg_norm','allbkg_norm',STop_norm_val+TTbar_norm_val+WJets_norm_val)
allbkg_norm.setError(TMath.Sqrt(VV_norm.getError()*VV_norm.getError()+STop_norm.getError()*STop_norm.getError()+TTbar_norm.getError()*TTbar_norm.getError()+WJets_norm.getError()*WJets_norm.getError()))
allbkg_norm_val	= allbkg_norm.getVal()
allbkg_norm_err_val = allbkg_norm.getError()






w_err.var('cwww').setVal(0)
w_err.var('ccw').setVal(0)
w_err.var('cb').setVal(0)
Sig_norm	= RooRealVar('Sig_norm','Sig_norm',w_err.function('normfactor_3d_%s_%s'%(cat,ch)).getVal()*VV_norm_val)
#Sig_norm	= RooRealVar('Sig_norm','Sig_norm',w_SM.function('normfactor_3d_4fit_%s'%channel).getVal()*w_SM.data('SMdatahist').sumEntries())




allbkgsig_pdf	= RooAddPdf('allbkgsig','allbkgsig',RooArgList(allbkg_pdf,Sig_pdf),RooArgList(allbkg_norm,Sig_norm))



floatparams	= RooArgList(w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig0'%(ch,cat[1])),w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig1'%(ch,cat[1])),w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig2'%(ch,cat[1])),w_err.var('Deco_WJets0_xww_sim_%s_HP%s_mlvj_13TeV_eig3'%(ch,cat[1])),w_err.var('Deco_TTbar_xww_signal_region_%s_HP%s_mlvj_13TeV_eig0'%(ch,cat[1])),w_err.var('Deco_TTbar_xww_signal_region_%s_HP%s_mlvj_13TeV_eig1'%(ch,cat[1])))

c		= TCanvas('c','c',1)
c.cd()
p		= rrv_x.frame(900,3500)


c.SetLogy()
#allbkgsig_pdf.plotOn(p, RooFit.Normalization(norm4error.getVal()), RooFit.DrawOption('F'), RooFit.FillColor(kCyan+1))
#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s,TTbar_xww_%s_HP%s,WJets_xww_%s_HP%s'%(ch,cat[1],ch,cat[1],ch,cat[1])), RooFit.DrawOption('F'), RooFit.FillColor(kGreen+1), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s,TTbar_xww_%s_HP%s'%(ch,cat[1],ch,cat[1])), RooFit.DrawOption('F'), RooFit.FillColor(kOrange), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s'%(ch,cat[1])), RooFit.DrawOption('F'), RooFit.FillColor(kRed), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s'%(ch,cat[1])), RooFit.DrawOption('F'), RooFit.FillColor(kBlue), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop,TTbar,WJets'), RooFit.DrawOption('F'), RooFit.FillColor(kGreen+1), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop,TTbar'), RooFit.DrawOption('F'), RooFit.FillColor(kOrange), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop'), RooFit.DrawOption('F'), RooFit.FillColor(kRed), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop'), RooFit.DrawOption('F'), RooFit.FillColor(kBlue), RooFit.LineColor(kBlack), RooFit.LineWidth(1))


#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s,TTbar_xww_%s_HP%s,WJets_xww_%s_HP%s'%(ch,cat[1],ch,cat[1],ch,cat[1])), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s,TTbar_xww_%s_HP%s'%(ch,cat[1],ch,cat[1])), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s'%(ch,cat[1])), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
#allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop_xww_%s_HP%s'%(ch,cat[1])), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.LineStyle(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop,TTbar,WJets'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop,TTbar'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('Sig,STop'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))
allbkg_pdf.plotOn(p, RooFit.Normalization(allbkg_norm_val), RooFit.Components('STop'), RooFit.LineColor(kBlack), RooFit.LineWidth(1))

draw_error_band(allbkg_pdf, rrv_x.GetName(), allbkg_norm, floatparams , w_err, p, kBlack, "F", 3013)
#allbkg_norm.setError(0)
#draw_error_band(allbkg_pdf, rrv_x.GetName(), allbkg_norm, floatparams , w_err, p, kGreen+1, "F", 3001)
#allbkg_norm.setError(allbkg_norm_err_val)
slope_nuis	= w_err.var('slope_nuis')
slope_nuis.setError(0.05)



#draw_error_band(allbkgsig_pdf, rrv_x.GetName(), norm4error, RooArgList(slope_nuis), w_err, p, kWhite, "F", 3144)
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
allbkgsig_pdf.plotOn(p,RooFit.Normalization(norm4error.getVal(),RooAbsReal.NumEvent), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.LineStyle(2), RooFit.LineColor(kViolet))#, RooFit.LineStyle(kDashed))
rrv_x.setRange('hi',900,3500)
#Sig_pdf.plotOn(p,RooFit.Normalization(Sig_norm.getVal(),RooAbsReal.NumEvent), RooFit.LineColor(kViolet))
#VV_pdf.plotOn(p,RooFit.Normalization(VV_norm_val,RooAbsReal.NumEvent), RooFit.LineColor(kOrange))
#VV_pdf.plotOn(p,RooFit.Normalization(VV_norm_val,RooAbsReal.NumEvent), RooFit.LineColor(kOrange))
w_sig.var('slope_nuis').setVal(1)


leg	= TLegend(0.65,0.5,0.89,0.89)
#leg.AddEntry(p.getObject(9),'bkg + signal %s'%latex_par[options.POI],'l')
leg.AddEntry(p.getObject(9),'signal %s'%latex_par[options.POI],'l')
leg.AddEntry(p.getObject(0),'W+Jets','F')
leg.AddEntry(p.getObject(1),'t#bar{t}','F')
leg.AddEntry(p.getObject(2),'VV','F')
leg.AddEntry(p.getObject(3),'single top','F')
leg.AddEntry(p.getObject(8),'background uncertainty','F')
leg.SetBorderSize(0)
leg.SetFillColor(0)
if ch == 'el':
	leg.SetHeader('e#nu,%s-category'%cat)
else:
	leg.SetHeader('#mu#nu,%s-category'%cat)




p.GetYaxis().SetRangeUser(0.002,5e3)
p.GetYaxis().SetTitle('events / 100 GeV')
p.SetTitle('')
p.Draw()
leg.Draw("SAME")
CMS_lumi.lumi_13TeV 	= '2.3 fb^{-1}'
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText	= 'preliminary'
CMS_lumi.CMS_lumi(c,4,11)
c.Draw()
c.Update()
#c.SaveAs('%s_%s_%s.png'%(cat,ch,options.POI))

c2=TCanvas('c2','c2',1)
c2.cd()
if ch=='el':
	plothi = 80
elif ch=='mu':
	plothi = 110
p.GetYaxis().SetRangeUser(0,plothi)
p.GetXaxis().SetRangeUser(900,1500)
p.Draw()
leg.Draw("SAME")
CMS_lumi.CMS_lumi(c2,4,11)
c2.Draw()
#c2.SaveAs('%s_%s_%s_lin.png'%(cat,ch,options.POI))

print allbkg_norm.getError()/WJets_norm.getVal()
print Sig_norm.getVal(),' / ',VV_norm_val
print VV_norm_val
print STop_norm_val
print TTbar_norm_val
print WJets_norm_val
raw_input('.,')

