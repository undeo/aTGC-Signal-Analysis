import os
cmssw_base = os.environ['CMSSW_BASE']

import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex, TClass
import subprocess

def cmsLabel(canvas, lumi, prelim = False, lumiLabel = 'fb', s = 8.):
    l = TLatex();
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextAlign(31);
    l.SetTextSize(0.045);

    canvas.cd()
    prelimText = ''
    if prelim:
        prelimText = ' preliminary'
    l.DrawLatex(1. - canvas.GetRightMargin(), 1. - canvas.GetTopMargin() + 0.01,
                'CMS%s, #scale[0.5]{#lower[-0.15]{#it{#int}}}#it{L} dt = %0.1f#kern[0.2]{%s}^{-1}, #sqrt{#it{s}} = %.0f#kern[0.1]{TeV}' % \
                (prelimText, lumi, lumiLabel, s)
                )
    canvas.Update()

def cmsPrelim(canvas, lumi):
    cmsLabel(canvas, lumi, True)

def RooFitInclude():
    print "adding RooFit ...",
    scramCmd = ['scram','tool','info','roofitcore']
    grepCmd = ['grep', 'INCLUDE']
    pscram = subprocess.Popen(scramCmd, stdout = subprocess.PIPE)
    pgrep = subprocess.Popen(grepCmd, stdin=pscram.stdout,
                             stdout=subprocess.PIPE)
    pscram.stdout.close()
    output = pgrep.communicate()[0]
    if (pgrep.returncode == 0):
        print "done"
        roofitinc = output.split("=")[1].rstrip()
        # print roofitinc
        gROOT.GetInterpreter().AddIncludePath(roofitinc)
        roofitinc = '-I"' + roofitinc + '"'
        # print roofitinc
        gSystem.AddIncludePath(roofitinc)
        return True
    else:
        print "failed"
        print 'scram returned:',pscram.returncode,'grep:',pgrep.returncode
        return False

import atexit
import readline
import rlcompleter

historyPath = os.path.expanduser("~/.pyhistory")

def save_history(historyPath=historyPath):
    import readline
    readline.set_history_length(1000)
    readline.write_history_file(historyPath)

if os.path.exists(historyPath):
    readline.read_history_file(historyPath)

atexit.register(save_history)
## macroPath = gROOT.GetMacroPath()
## macroPath += os.environ['CMSSW_BASE'] + '/src/ElectroWeakAnalysis/VPlusJets/test:'
## gROOT.SetMacroPath(macroPath)
del atexit, readline, rlcompleter, save_history, historyPath

gROOT.SetStyle('Plain')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetOptStat("iouRMe")
gStyle.SetPalette(1)
gStyle.SetOptFit(1112)
gStyle.SetOptTitle(0)

gStyle.SetCanvasDefH(600) ## Height of canvas
gStyle.SetCanvasDefW(600) ## Width of canvas
gStyle.SetErrorX(0.)

gStyle.SetMarkerStyle(20)

## For the fit/function:
gStyle.SetFuncColor(2)
gStyle.SetFuncStyle(1)
gStyle.SetFuncWidth(1)

##  Margins:
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadLeftMargin(0.20) ## was 0.16
gStyle.SetPadRightMargin(0.06)## was 0.02

gStyle.SetTitleColor(1, "XYZ")
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.07, "XYZ")
gStyle.SetTitleXOffset(0.9)
gStyle.SetTitleYOffset(1.4) ## was 1.25

##  For the axis labels:
gStyle.SetLabelColor(1, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelOffset(0.007, "XYZ")
gStyle.SetLabelSize(0.06, "XYZ")
gStyle.SetNdivisions(505, "XYZ")

if (gSystem.DynamicPathName("libFWCoreFWLite.so",True)):
    gSystem.Load('libHiggsAnalysisCombinedLimit')
    gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMMozerpowhegweight.so")
    res = gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    gROOT.GetInterpreter().AddIncludePath(cmssw_base + '/src')
    gSystem.AddIncludePath('-I"' + cmssw_base + '/src"')
    if not RooFitInclude():
        workingdir = os.getcwd()
        print 'changing to', cmssw_base, 'directory'
        os.chdir(cmssw_base)
        RooFitInclude()
        print 'returning to working directory', workingdir
        os.chdir(workingdir)
    gROOT.ProcessLine('.L EffTableReader.cc+')
    gROOT.ProcessLine('.L EffTableLoader.cc+')
    gROOT.ProcessLine('.L CPWeighter.cc+')
    
    if not TClass.GetClass('RooPowerLaw'):
        # print 'importing RooFit PDF classes'
        gROOT.ProcessLine('.L RooPowerLaw.cc+')
    if not TClass.GetClass('RooPowerExpPdf'):
        gROOT.ProcessLine('.L RooPowerExpPdf.cxx+')
    if not TClass.GetClass('RooErfExpPdf'):
        gROOT.ProcessLine('.L RooErfExpPdf.cxx+')
    if not TClass.GetClass('RooErfPdf'):
        gROOT.ProcessLine('.L RooErfPdf.cxx+')
    if not TClass.GetClass('RooTH1DPdf'):
        gROOT.ProcessLine('.L RooTH1DPdf.cxx+')
    if not TClass.GetClass('RooChebyshevPdf'):
        gROOT.ProcessLine('.L RooChebyshevPDF.cc+')
    if not TClass.GetClass('alphaFunction'):
        gROOT.ProcessLine('.L alphaFunction.cxx+')        
    
print 'end of pyroot_logon'

if __name__ == '__main__':
    from ROOT import *
