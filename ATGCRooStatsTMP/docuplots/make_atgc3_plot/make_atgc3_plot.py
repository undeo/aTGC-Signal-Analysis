from ROOT import *





fileInWWel	= TFile.Open('WW_el.root')
ws_WWel		= fileInWWel.Get('wtmp')
fileInWWmu	= TFile.Open('WW_mu.root')
ws_WWmu		= fileInWWmu.Get('wtmp')
fileInWZel	= TFile.Open('WZ_el.root')
ws_WZel		= fileInWZel.Get('wtmp')
fileInWZmu	= TFile.Open('WZ_mu.root')
ws_WZmu		= fileInWZmu.Get('wtmp')



c1	= TCanvas('c1','c1',1)
c1.SetLogy()

p1	= ws_WWel.var('rrv_mass_lvj').frame()
p1.SetTitle('')

ws_WWel.var('cwww').setVal(-12)
ws_WWel.var('ccw').setVal(-20)
ws_WWel.var('cb').setVal(-60)
ws_WWel.data('datahist_all3').plotOn(p1,RooFit.LineColor(kBlue),RooFit.LineWidth(1),RooFit.DrawOption('E'))
ws_WWel.pdf('aTGC_model_WW_el').plotOn(p1,RooFit.LineColor(kBlue),RooFit.Normalization(ws_WWel.function('normfactor_3d_WW_el').getVal()*ws_WWel.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))

ws_WZel.var('cwww').setVal(-12)
ws_WZel.var('ccw').setVal(-20)
ws_WZel.var('cb').setVal(-60)
ws_WZel.data('datahist_all3').plotOn(p1,RooFit.LineColor(kRed),RooFit.LineWidth(1),RooFit.DrawOption('E'))
ws_WZel.pdf('aTGC_model_WZ_el').plotOn(p1,RooFit.LineColor(kRed),RooFit.Normalization(ws_WZel.function('normfactor_3d_WZ_el').getVal()*ws_WZel.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))

l1	= TLegend(0.6,0.6,0.89,0.89)
l1.SetFillColor(0)
l1.SetBorderSize(0)
l1.SetHeader('electron channel')
l1.AddEntry(p1.getObject(0),'MC WW-category','le')
l1.AddEntry(p1.getObject(1),'signal model WW-category','l')
l1.AddEntry(p1.getObject(2),'MC WZ-category','le')
l1.AddEntry(p1.getObject(3),'signal model WZ-category','l')

p1.GetYaxis().SetRangeUser(1,100)
p1.Draw()
c1.Draw()
l1.Draw()
c1.Update()
c1.SaveAs('../atgc3_el.pdf')
c1.SaveAs('../atgc3_el.png')


c2	= TCanvas('c2','c2',1)
c2.SetLogy()

p2	= ws_WWmu.var('rrv_mass_lvj').frame()
p2.SetTitle('')

ws_WWmu.var('cwww').setVal(-12)
ws_WWmu.var('ccw').setVal(-20)
ws_WWmu.var('cb').setVal(-60)
ws_WWmu.data('datahist_all3').plotOn(p2,RooFit.LineColor(kBlue),RooFit.LineWidth(1),RooFit.DrawOption('E'))
ws_WWmu.pdf('aTGC_model_WW_mu').plotOn(p2,RooFit.LineColor(kBlue),RooFit.Normalization(ws_WWmu.function('normfactor_3d_WW_mu').getVal()*ws_WWmu.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))

ws_WZmu.var('cwww').setVal(-12)
ws_WZmu.var('ccw').setVal(-20)
ws_WZmu.var('cb').setVal(-60)
ws_WZmu.data('datahist_all3').plotOn(p2,RooFit.LineColor(kRed),RooFit.LineWidth(1),RooFit.DrawOption('E'))
ws_WZmu.pdf('aTGC_model_WZ_mu').plotOn(p2,RooFit.LineColor(kRed),RooFit.Normalization(ws_WZmu.function('normfactor_3d_WZ_mu').getVal()*ws_WZmu.data('SMdatahist').sumEntries(),RooAbsReal.NumEvent))

l2	= TLegend(0.6,0.6,0.89,0.89)
l2.SetFillColor(0)
l2.SetBorderSize(0)
l2.SetHeader('muon channel')
l2.AddEntry(p2.getObject(0),'MC WW-category','le')
l2.AddEntry(p2.getObject(1),'signal model WW-category','l')
l2.AddEntry(p2.getObject(2),'MC WZ-category','le')
l2.AddEntry(p2.getObject(3),'signal model WZ-category','l')

p2.GetYaxis().SetRangeUser(1,100)
p2.Draw()
c2.Draw()
l2.Draw()
c2.Update()
c2.SaveAs('../atgc3_mu.pdf')
c2.SaveAs('../atgc3_mu.png')
raw_input("...")



















