#! /usr/bin/env python

# If you want to edit paddings or text in the plot
# go to the bottom of this text file.

#root crap
import ROOT
from ROOT import TGraph2D, TGraph, TGraphPolar, TFile, TCanvas, TMath, TLegend, TPaveText, TGraphAsymmErrors, TLine

#python modules
from glob import glob
from array import array
from math import pi, fabs
from ctypes import c_double
import string, re, random
from optparse import OptionParser

ROOT.gROOT.ProcessLine('.X CMSStyle.C')

parser = OptionParser(description="%prog : Creates 2D contour plots from Higgs Combined Limit output.",
                      usage="make_hcl_atgc_contour.py --input=dir/*.root")
parser.add_option('--input',dest='input',help='The input set of root files "some/dir/*.root"')
parser.add_option('--par1',dest='par1',help='The name of the x-axis parameter')
parser.add_option('--par2',dest='par2',help='The name of the y-axis parameter')
parser.add_option('--par1Latex',dest='par1Latex',default="NotSet",help='The name of the x-axis parameter')
parser.add_option('--par2Latex',dest='par2Latex',default="NotSet",help='The name of the y-axis parameter')
parser.add_option('--flavorText',dest='flavorText',default="",help='Flavor text for the final plot')
(options,args) = parser.parse_args()

if options.input is None:
    print 'You must specify --input to define the input scanning grid.'
    exit(1)

whichParm = None

if options.par1 is None:
    print 'You must specify --par1'
    exit(1)
    
if options.par2 is None:
    print 'You must specify --par2'
    exit(1)
        
def extractParValue(par,file):
    pattern = r'%s\-?[0-9]*\.?[0-9]*e?\-?[0-9]*'%par
    match = re.search(pattern,file)
    val = 0.0

    if match is not None:
        print file
        print match.group(0),match.group(0).split(par)[-1]
        
        temp = match.group(0)
        val = float(temp.split(par)[-1])
        #print val
    else:
        val = -1
        print "Could not find a valid value for par %s in filename: %s !"%(par,file)        

    return val

ROOT.gStyle.SetPalette(1)

par1 = options.par1
par2 = options.par2

par1Latex = options.par1Latex
par2Latex = options.par2Latex

limits = {'-2s' : [TGraph2D(),0],
          '-1s' : [TGraph2D(),1],
          'mean': [TGraph2D(),2],
          '+1s' : [TGraph2D(),3],
          '+2s' : [TGraph2D(),4],
          'obs' : [TGraph2D(),5]}

n = 0
for file in glob(options.input):
    if '.root' not in file:
        print "Skipped non-root file: %s"%file
        continue
    
    par1val = extractParValue(par1,file)
    par2val = extractParValue(par2,file)
    
    if par1val == -1 or par2val == -2:
        continue
    
    pwd = ROOT.gDirectory.GetPath()
    f = TFile.Open(file)
    ROOT.gDirectory.cd(pwd)
    if f.Get("limit"):
        lTree = f.Get("limit").CloneTree()
        f.Close()
    else:
        f.Close()
        continue     
    
    for it in limits:
        graph = limits[it][0]
        idx = limits[it][1]
        lTree.GetEntry(idx)
        
        rval = lTree.GetLeaf("limit").GetValue()        
        graph.SetPoint(n,par1val,par2val,rval)
        
    n+=1

canv2 = TCanvas('two','two',500,500)
canv2.cd()

limits['+2s'][0].Draw('TRI')


canv = TCanvas('dummy','dummy',500,500)
limits['+2s'][0].Draw("cont 5z list")    

limconts = {'-2s' : 0,
            '-1s' : 0,
            'mean': 0,
            '+1s' : 0,
            '+2s' : 0,
            'obs' : 0}

#process TGraph2Ds into contours at r == 1
for it in limits:
    limits[it][0].Draw("cont 2z list")    
    limits[it][0].GetHistogram().SetContour(1,array('d',[1.0]))    
    canv.Update()    
        
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    print it, contours.At(0).First()
    limconts[it] = contours.At(0).First().Clone()
    
del canv

#from here were build the two-dimensional aTGC limit

finalPlot = TCanvas('final','limits',500,500)
finalPlot.cd()

limconts['+2s'].SetLineColor(ROOT.kGreen)
limconts['+2s'].SetFillColor(ROOT.kGreen)
limconts['+2s'].Draw("ACF")

limconts['+1s'].SetLineColor(ROOT.kYellow)
limconts['+1s'].SetFillColor(ROOT.kYellow)
limconts['+1s'].Draw("SAME CF")

limconts['-1s'].SetLineColor(ROOT.kGreen)
limconts['-1s'].SetFillColor(ROOT.kGreen)
limconts['-1s'].Draw("SAME CF")

limconts['-2s'].SetFillColor(ROOT.kWhite)
limconts['-2s'].SetLineColor(ROOT.kGreen)
limconts['-2s'].Draw("SAME CF")

limconts['mean'].SetLineColor(ROOT.kBlack)
limconts['mean'].SetLineWidth(2)
limconts['mean'].SetLineStyle(2)
limconts['mean'].Draw("SAME C")

limconts['obs'].SetLineColor(ROOT.kBlack)
limconts['obs'].SetLineWidth(2)
limconts['obs'].Draw("SAME C")

limconts['+2s'].GetYaxis().SetRangeUser(-1.25*limconts['+2s'].GetYaxis().GetXmax(),
                                        +2.0*limconts['+2s'].GetYaxis().GetXmax())

limconts['+2s'].SetTitle()
limconts['+2s'].GetXaxis().SetTitle(par1Latex)
limconts['+2s'].GetXaxis().SetTitleFont(132)
limconts['+2s'].GetYaxis().SetTitle(par2Latex)
limconts['+2s'].GetYaxis().SetTitleFont(132)
limconts['+2s'].GetYaxis().SetTitleOffset(1.20)

limconts['+2s'].GetXaxis().SetRangeUser(-0.1,0.1)
limconts['+2s'].GetYaxis().SetRangeUser(-0.3,0.5)

SMpoint = TGraph(1)
SMpoint.SetPoint(1,0,0)
SMpoint.Draw("SAME Po")

smLabel = TPaveText(0,
                    limconts['-2s'].GetYaxis().GetXmax()/8,
                    limconts['-2s'].GetXaxis().GetXmax()/3.5,
                    -limconts['-2s'].GetYaxis().GetXmax()/8)
smLabel.SetFillStyle(0)
smLabel.SetBorderSize(0)
smLabel.AddText(' SM')
smLabel.Draw()

legend = TLegend(0.212,0.686,0.554,0.917,"","NDC")
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetHeader("#bf{#it{CMS}} Preliminary")

legend.AddEntry(limconts['obs'],"Observed","L")
legend.AddEntry(limconts['mean'],"Expected","L")
legend.AddEntry(limconts['+1s'],"#pm 1#sigma","F")
legend.AddEntry(limconts['+2s'],"#pm 2#sigma","F")

legend.Draw()

text = TPaveText(0.516,0.720,0.915,0.951,"NDC")
text.SetFillStyle(0)
text.SetBorderSize(0)
text.AddText("95% CLs Limit on "+"#bf{%s} and #bf{%s}"%(par1Latex,par2Latex))
text.AddText(0,0.35,"#intL dt= 4.7 fb^{-1}, #sqrt{s} = 7 TeV")
text.Draw()

text2 = TPaveText(0.155,0.199,0.974,0.244,"NDC")
text2.SetFillStyle(0)
text2.SetBorderSize(0)
text2.AddText("aTGC values outside contour excluded")
text2.Draw()

text3 = TPaveText(0.506,0.699,0.905,0.758,"NDC")
text3.SetFillStyle(0)
text3.SetBorderSize(0)
text3.AddText(options.flavorText)
text3.Draw()    

finalPlot.RedrawAxis()
finalPlot.ResetAttPad()
finalPlot.Update()

finalPlot.Draw()
finalPlot.Print("%s_%s_2dlimit.pdf"%(options.par1,options.par2))
finalPlot.Print("%s_%s_2dlimit.eps"%(options.par1,options.par2))
finalPlot.Print("%s_%s_2dlimit.png"%(options.par1,options.par2))

def make1DLimit(limits2d, par, parlatex, parmin, parmax, samples, boundScale, isX):
    limitplot = TCanvas("%slimit"%par,"%s limit"%par,500,500)

    limits1d = {'2s' :  [TGraphAsymmErrors(),0],
              '1s' :  [TGraphAsymmErrors(),1],
              'mean': [TGraphAsymmErrors(),2],
              'obs' : [TGraphAsymmErrors(),3]}

    parSize = parmax - parmin

    lowerLimit = 0.0
    lowerLimitFound = False
    upperLimit = 0.0
    upperLimitFound = False

    cutoff = parmin*boundScale
    bound = 0.0
    if isX:
        bound = limits2d['-2s'][0].Interpolate(cutoff,0.0)
    else:
        bound = limits2d['-2s'][0].Interpolate(0.0,cutoff)
    
    for i in range(samples): 
        parval = parmin + i*parSize/samples
        
        if isX:
            obs  = limits2d['obs'][0].Interpolate(parval,0.0)
            mean = limits2d['mean'][0].Interpolate(parval,0.0)
            p1s = limits2d['+1s'][0].Interpolate(parval,0.0)
            m1s = limits2d['-1s'][0].Interpolate(parval,0.0)
            p2s = limits2d['+2s'][0].Interpolate(parval,0.0)
            m2s = limits2d['-2s'][0].Interpolate(parval,0.0)
        else:
            obs  = limits2d['obs'][0].Interpolate(0.0,parval)
            mean = limits2d['mean'][0].Interpolate(0.0,parval)
            p1s = limits2d['+1s'][0].Interpolate(0.0,parval)
            m1s = limits2d['-1s'][0].Interpolate(0.0,parval)
            p2s = limits2d['+2s'][0].Interpolate(0.0,parval)
            m2s = limits2d['-2s'][0].Interpolate(0.0,parval)


        #print m2s, m1s, mean, p1s, p2s
        
        if obs > 1 and not lowerLimitFound and parval < 0:
            lowerLimit = parval
            lowerLimitFound = True

        if obs < 1 and not upperLimitFound and parval > 0:
            upperLimit = parval
            upperLimitFound = True

        if m2s < bound  and fabs(parval) > fabs(cutoff):
            #par1 observed limit
            limits1d['obs'][0].SetPoint(i,parval,obs)
            #mean and one sigma expected
            limits1d['mean'][0].SetPoint(i,parval,mean)
            limits1d['1s'][0].SetPoint(i,parval,mean)
            limits1d['1s'][0].SetPointError(i,0,0,mean-m1s,p1s-mean)
            #two sigma expected    
            limits1d['2s'][0].SetPoint(i,parval,mean)
            limits1d['2s'][0].SetPointError(i,0,0,mean-m2s,p2s-mean)
        else:
            limits1d['obs'][0].SetPoint(i,parval,bound+0.1)
            limits1d['mean'][0].SetPoint(i,parval,bound+0.1)
            limits1d['1s'][0].SetPoint(i,parval,bound+0.1)
            limits1d['1s'][0].SetPointError(i,0,0,0,0)
            limits1d['2s'][0].SetPoint(i,parval,bound+0.1)
            limits1d['2s'][0].SetPointError(i,0,0,0,0)

    print "95% CL on"+" %s = [%.3g,%.3g]"%(par,lowerLimit,upperLimit)

    limitplot.cd()
    limitplot.SetLogy()

    limits1d['2s'][0].SetFillColor(ROOT.kGreen)
    limits1d['2s'][0].Draw("A E3")
    
    limits1d['1s'][0].SetFillColor(ROOT.kYellow)
    limits1d['1s'][0].Draw("SAME E3")
    
    limits1d['mean'][0].SetLineStyle(2)
    limits1d['mean'][0].SetLineWidth(2)
    limits1d['mean'][0].Draw("SAME C")
    
    limits1d['obs'][0].SetLineWidth(2)
    limits1d['obs'][0].Draw("SAME C")
    
    #titles
    limits1d['2s'][0].GetYaxis().SetTitle("95% CL limit on #sigma/#sigma_{aTGC}")
    limits1d['2s'][0].GetYaxis().SetTitleFont(132)
    limits1d['2s'][0].GetXaxis().SetTitle(parlatex)
    limits1d['2s'][0].GetXaxis().SetTitleFont(132)
    
    limits1d['2s'][0].GetYaxis().SetRangeUser(limits1d['2s'][0].GetYaxis().GetXmin()*0.75,
                                              bound)
    limits1d['2s'][0].GetXaxis().SetRangeUser(parmin*0.985,parmax*0.96)
    
    legend.SetX1NDC(0.183)
    legend.SetY1NDC(0.699)
    legend.SetX2NDC(0.524)
    legend.SetY2NDC(0.930)
    legend.Draw()
    
    text1d = TPaveText(0.359,0.265,0.758,0.496,"NDC")
    text1d.SetFillStyle(0)
    text1d.SetBorderSize(0)
    text1d.AddText("95% CLs Limit on "+"#bf{%s}"%(parlatex))
    text1d.AddText(0,0.35,"#intL dt= 4.7 fb^{-1}, #sqrt{s} = 7 TeV")
    text1d.Draw()
    
    text3.SetX1NDC(0.357)
    text3.SetY1NDC(0.246)
    text3.SetX2NDC(0.756)
    text3.SetY2NDC(0.305)
    text3.Draw()
    
    textlim = TPaveText(0.357,0.197,0.754,0.246,"NDC")
    textlim.SetFillStyle(0)
    textlim.SetBorderSize(0)
    textlim.AddText("%.2g < %s  < %.2g"%(lowerLimit,parlatex,upperLimit))
    textlim.Draw()
    
    #lowLimitLine = TLine(lowerLimit,limits1d['2s'][0].GetYaxis().GetXmin()*0.75,
    #                     lowerLimit,1)
    #lowLimitLine.SetLineColor(14)
    #lowLimitLine.SetLineWidth(2)
    #lowLimitLine.Draw()
    #upLimitLine = TLine(upperLimit,limits1d['2s'][0].GetYaxis().GetXmin()*0.75,
    #                    upperLimit,1)
    #upLimitLine.SetLineColor(14)
    #upLimitLine.SetLineWidth(2)
    #upLimitLine.Draw()

    oneLine = TLine(parmin*0.985,1,parmax*0.96,1)
    oneLine.SetLineStyle(9)
    oneLine.SetLineColor(14)
    oneLine.Draw()
    
    limitplot.Draw()
    
    return limitplot.Clone()

# now, do the 1-d limit for par1 and par2

par1Max = limconts['+2s'].GetXaxis().GetXmax()
par1Min = limconts['+2s'].GetXaxis().GetXmin()
par2Max = limconts['+2s'].GetYaxis().GetXmax()
par2Min = limconts['+2s'].GetYaxis().GetXmin()
samples = 1000

par1Plot = make1DLimit(limits, par1, par1Latex, par1Min, par1Max, samples, 0.17, True)
par1Plot.Draw()
par1Plot.Print("%s_1dlimit.pdf"%par1)
par1Plot.Print("%s_1dlimit.eps"%par1)
par1Plot.Print("%s_1dlimit.png"%par1)

par2Plot = make1DLimit(limits, par2, par2Latex, par2Min, par2Max, samples, 0.26, False)
par2Plot.Draw()
par2Plot.Print("%s_1dlimit.pdf"%par2)
par2Plot.Print("%s_1dlimit.eps"%par2)
par2Plot.Print("%s_1dlimit.png"%par2)
