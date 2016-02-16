import pyroot_logon
import limits
import os
import sys

from array import *

from ROOT import *
from optparse import OptionParser
from ConfigParser import SafeConfigParser

parser = OptionParser(description="%prog : A RooStats Implementation of Anomalous Triple Gauge Coupling Analysis.",
                      usage="buildWZworkspace --config=example_config.cfg")
cfgparse = SafeConfigParser()

parser.add_option("--config",dest="config",help="The name of the input configuration file.")
(options,args) = parser.parse_args()

miss_options = False

if options.config is None:
    print 'Need to specify --config'
    miss_options=True
    
if miss_options:
    exit(1)
        
cfgparse.read(options.config)
options.config = cfgparse # put the parsed config file into our options

cfg = options.config

norm_sig_sm = -1
norm_sig_sm_up = -1
norm_sig_sm_down = -1
norm_bkg = -1
norm_obs = -1

fit_sections = cfg.sections()

signalModel ='-1'
signalModel = cfg.get('Global','model')     
print 'modelName= ',signalModel

fit_sections.remove('Global') #don't need to iterate over the global configuration

limtype = -1
Ndim = -1

if ( signalModel == 'par1_TH1'):
    limtype = 0
    Ndim = 1
elif ( signalModel == 'par1_TF1'):
    limtype = 1
    Ndim = 1
elif ( signalModel == 'par1par2_TF2par'):
    limtype = 0
    Ndim = 2
elif ( signalModel == 'par1par2_TH2'):
    limtype = 1
    Ndim = 2
elif ( signalModel == 'par1par2_TF2'):
    limtype = 2
    Ndim = 2
elif ( signalModel == 'par1par2par3_TH3'):
    limtype = 0
    Ndim = 3
elif ( signalModel == 'par1par2par3_TF3'):
    limtype = 1
    Ndim = 3
else:
    raise RuntimeError('InvalidCouplingChoice %s'%signalModel,
                       'We can only use 1D (par1_TH1, par1_TF1)  2D (par1par2_TF2par, par1par2_TH2, par1par2_TF2) or 3D (par1par2par3_TH3, par1par2par3_TF3) models right now!')

par1name = 'par1'
par1low = None
par1high = None

par2name ='par2'
par2low = None
par2high = None

par3name = 'par3'
par3low = None
par3high = None

NlnN = 0

cfg_items=cfg.items('Global')
for cfg_item in cfg_items:
    if 'par1name' in cfg_item:
        par1name = cfg.get('Global','par1name')     
    if 'par2name' in cfg_item:
        par2name = cfg.get('Global','par2name')     
    if 'par3name' in cfg_item:
        par3name = cfg.get('Global','par3name')     
    if 'nlnn' in cfg_item:
        NlnN = int(cfg.get('Global','nlnn'))     

lnN_name = []
for i in range(1,NlnN+1):
    lnN_name.append(cfg.get('Global','lnN%i_name'%i))
lnN_value = []
for i in range(0,NlnN):
    lnN_value.append([])
    for value in cfg.get('Global','lnN%i_value'%(i+1)).split(','):
        lnN_value[i].append(value)


lnN_for = []
for i in range(0,NlnN):
    lnN_for.append([])
    for name in cfg.get('Global','lnN%i_for'%(i+1)).split(','):
        lnN_for[i].append(name)

print '\n\t\t=============================================> lnN: ',lnN_name
print 'lnN value: ',lnN_value
print 'lnN for: ',lnN_for



if Ndim == 1:
    par1low  = float(cfg.get('Global', 'par1Low'))
    par1high = float(cfg.get('Global', 'par1High'))

elif Ndim == 2:
    par1low  = float(cfg.get('Global', 'par1Low'))
    par1high = float(cfg.get('Global', 'par1High'))
    par2low  = float(cfg.get('Global', 'par2Low'))
    par2high = float(cfg.get('Global', 'par2High'))

elif Ndim == 3:
    par1low  = float(cfg.get('Global', 'par1Low'))
    par1high = float(cfg.get('Global', 'par1High'))
    par2low  = float(cfg.get('Global', 'par2Low'))
    par2high = float(cfg.get('Global', 'par2High'))
    par3low  = float(cfg.get('Global', 'par3Low'))
    par3high = float(cfg.get('Global', 'par3High'))

else:
    print 'Only dimensions 1, 2 and 3 implemented... this thing will crash.'

NSigBkg_corr_unc_int=0

basepath = '%s/src/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'%os.environ['CMSSW_BASE']


for section in fit_sections:
    codename = section
    lType = codename
    f = TFile('%s/%s.root'%(basepath,codename))

    Nbkg = cfg.get(codename,'Nbkg')
    print "Nbkg= ",Nbkg
    Nbkg_int=int(Nbkg)

    bkg_name = []
    for i in range(1,Nbkg_int+1):
        bkg_name.append(cfg.get(codename,'bkg%i_name'%i))

    background = []

    for i in range(0,Nbkg_int):
        background.append(f.Get(bkg_name[i]))


    print 'backgrounds= ',background
    background_shapeSyst = []
    for i in range(0,Nbkg_int):
        background_shapeSyst.append([])
        cfg_items=cfg.items(codename)
        for cfg_item in cfg_items:
            if ('bkg%i_shape_syst'%(i+1) in cfg_item):
                for name in cfg.get(codename,'bkg%i_shape_syst'%(i+1)).split(','):
                    background_shapeSyst[i].append(name)


    background_backshapeUp = []
    background_backshapeDown = []

    for j in range(0,Nbkg_int):
        background_backshapeUp.append([])
        background_backshapeDown.append([])
        for i in range(0,len(background_shapeSyst[j])):
            background_backshapeUp[j].append(f.Get('%sUp'%background_shapeSyst[j][i]))
            background_backshapeDown[j].append(f.Get('%sDown'%background_shapeSyst[j][i]))


    data_obs = f.Get('data_obs')
    diboson = f.Get('diboson')

    doSignalShape_unc=False
    doSignalBkg_corr_unc=False
    cfg_items=cfg.items(codename)
    for cfg_item in cfg_items:
        if 'signal_shape_syst' in cfg_item:
            doSignalShape_unc = True
        if 'nsigbkg_corr_unc' in cfg_item:
            doSignalBkg_corr_unc = True


    if (doSignalShape_unc):
        diboson_up = {}
        diboson_down = {}
        norm_sig_sm_up = {}
        norm_sig_sm_down = {}
        signal_shapeSyst = [string(i) for i in cfg.get(codename,'signal_shape_syst').split(',')]
        for i in range(0,len(signal_shapeSyst)):
            diboson_up[i] = f.Get('%sUp'%signal_shapeSyst[i])
            diboson_down[i] = f.Get('%sDown'%signal_shapeSyst[i])
            norm_sig_sm_up[i] = diboson_up[i].Integral()
            norm_sig_sm_down[i] = diboson_down[i].Integral()

    if (doSignalBkg_corr_unc):
        NSigBkg_corr_unc = cfg.get(codename,'NSigBkg_corr_unc')
        NSigBkg_corr_unc_int=int(NSigBkg_corr_unc)

        SignalBkg_corr_name_ws = []
        for i in range(1,NSigBkg_corr_unc_int+1):
            SignalBkg_corr_name_ws.append(cfg.get(codename,'correlated_SigBkg_unc%i_name'%i))


        SignalBkg_corr_name = {}
        for i in range(1,NSigBkg_corr_unc_int+1):
            SignalBkg_corr_name[i-1]=[string(j) for j in cfg.get(codename,'correlated_SigBkg_unc%s'%i).split(',')]

# check if shape uncertainty name is one of those where Signal and bkg uncertainties are correlated
    def isItCorrelated(name):
        isItCorr=False
        for i in range(0,NSigBkg_corr_unc_int):
            if (name in SignalBkg_corr_name[i]):
                isItCorr=True
        return isItCorr

    def isItCorrelated_name(name):
        name_out=name
        for i in range(0,NSigBkg_corr_unc_int):
            if (name in SignalBkg_corr_name[i]):
                name_out=SignalBkg_corr_name_ws[i]
        return name_out

    
    norm_sig_sm = diboson.Integral()
    norm_bkg = []
    for i in range(0,Nbkg_int):
        norm_bkg.append(background[i].Integral())
    norm_obs = data_obs.Integral()
    
    theWS = RooWorkspace('proc_%s'%codename, 'proc_%s'%codename)
    
    wpt = theWS.factory('observable_%s[%f,%f]' % (codename,data_obs.GetBinLowEdge(1), 
                                         data_obs.GetBinLowEdge(data_obs.GetNbinsX())+data_obs.GetBinWidth(data_obs.GetNbinsX())))

    binning=array('d',[])

    for i in range(1, data_obs.GetNbinsX()+1):
        binning.append(data_obs.GetBinLowEdge(i))
    binning.append(data_obs.GetBinLowEdge(data_obs.GetNbinsX()+1))


    bins=RooBinning(len(binning)-1, binning)

    wpt.setBinning(bins)

    if Ndim == 1:
        par1 = theWS.factory('%s[0., %f, %f]' % (par1name, par1low, par1high))
    elif Ndim == 2:
        par1 = theWS.factory('%s[0., %f, %f]' % (par1name, par1low, par1high))
        par2 = theWS.factory('%s[0., %f, %f]' % (par2name, par2low, par2high))
    elif Ndim == 3:
        par1 = theWS.factory('%s[0., %f, %f]' % (par1name, par1low, par1high))
        par2 = theWS.factory('%s[0., %f, %f]' % (par2name, par2low, par2high))
        par3 = theWS.factory('%s[0., %f, %f]' % (par3name, par3low, par2high))
    else:
        raise RuntimeError('Ndim = %s'%Ndim,
                           'We can only use 1D (par1_TH1, par1_TF1)  2D (par1par2_TF2par, par1par2_TH2, par1par2_TF2) or 3D (par1par2par3_TH3, par1par2par3_TF3) models right now!')

    vars = RooArgList(wpt)
    varSet = RooArgSet(wpt)

    data = RooDataHist('data_obs', 'data_obs_proc_%s'%codename, vars, data_obs)
    print 'data integral: ',data.Print()
    bkgHist = {}
    for i in range(0,Nbkg_int):
        bkgHist[i] = RooDataHist('anomalousCoupling_bkg%i_%s'%(i+1,codename),
                              'anomalousCoupling_bkg%i_%s'%(i+1,codename),
                              vars,
                              background[i])
    
    bkgHist_systUp = []
    bkgHist_systDown = []

    for j in range(0,Nbkg_int):
        bkgHist_systUp.append([])
        bkgHist_systDown.append([])
        for i in range(0,len(background_shapeSyst[j])):
            name_forCorr=background_shapeSyst[j][i]
            if (isItCorrelated(background_shapeSyst[j][i])):
                name_forCorr=isItCorrelated_name(background_shapeSyst[j][i])

            bkgHist_systUp[j].append(RooDataHist('anomalousCoupling_bkg%i_%s_%sUp'%(j+1,codename,name_forCorr),
                                                 'anomalousCoupling_bkg%i_%s_%sUp'%(j+1,codename,name_forCorr),
                                                 vars,
                                                 background_backshapeUp[j][i]))
            bkgHist_systDown[j].append(RooDataHist('anomalousCoupling_bkg%i_%s_%sDown'%(j+1,codename,name_forCorr),
                                                   'anomalousCoupling_bkg%i_%s_%sDown'%(j+1,codename,name_forCorr),
                                                   vars,
                                                   background_backshapeDown[j][i]))
   
    dibosonHist = RooDataHist('anomalousCoupling_SM_%s_rawshape'%codename,
                              'anomalousCoupling_SM_%s_rawshape'%codename,
                              vars,
                              diboson)

    if (doSignalShape_unc):
        dibosonHist_up = {}
        dibosonHist_down = {}
        for i in range(0,len(signal_shapeSyst)):
            name_forCorr=signal_shapeSyst[i]
            if (isItCorrelated(signal_shapeSyst[i])):
                name_forCorr=isItCorrelated_name(signal_shapeSyst[i])

            
            dibosonHist_up[i] = RooDataHist('anomalousCoupling_SM_%s_rawshape_%sUp'%(codename,name_forCorr),
                                         'anomalousCoupling_SM_%s_rawshape_%sUp'%(codename,name_forCorr),
                                         vars,
                                         diboson_up[i])
            dibosonHist_down[i] = RooDataHist('anomalousCoupling_SM_%s_rawshape_%sDown'%(codename,name_forCorr),
                                           'anomalousCoupling_SM_%s_rawshape_%sDown'%(codename,name_forCorr),
                                           vars,
                                           diboson_down[i])

    dibosonPdf = RooHistFunc('anomalousCoupling_SM_%s_shape'%codename,
                             'anomalousCoupling_SM_%s_shape'%codename,
                             varSet,
                             dibosonHist)

    if (doSignalShape_unc):
        dibosonPdf_up = {}
        dibosonPdf_down = {}
        for i in range(0,len(signal_shapeSyst)):


            name_forCorr=signal_shapeSyst[i]
            if (isItCorrelated(signal_shapeSyst[i])):
                name_forCorr=isItCorrelated_name(signal_shapeSyst[i])

            dibosonPdf_up[i] = RooHistFunc('anomalousCoupling_SM_%s_shape_%sUp'%(codename,name_forCorr),
                                        'anomalousCoupling_SM_%s_shape_%sUp'%(codename,name_forCorr),
                                        varSet,
                                        dibosonHist_up[i])
            dibosonPdf_down[i] = RooHistFunc('anomalousCoupling_SM_%s_shape_%sDown'%(codename,name_forCorr),
                                          'anomalousCoupling_SM_%s_shape_%sDown'%(codename,name_forCorr),
                                          varSet,
                                          dibosonHist_down[i])
    

    print '$$$$$$$$$$$$$ \t $$$$$$$$$  setting up for %s plane!'%signalModel,' limitType: ', limtype, ' Ndim: %s ',Ndim
   

    if (doSignalShape_unc):
        kappaLow = {}
        kappaHigh = {}
        aTGCPdf_norm = {}
        theta = {}
        kappaLow_sum_d = 1.
        kappaHigh_sum_d = 1.
        
        for i in range(0,len(signal_shapeSyst)):
            kappaLow[i] = RooRealVar("kappaL_%s_%s"%(i+1,codename),"kappaL_%s_%s"%(i+1,codename),norm_sig_sm_down[i]/norm_sig_sm)
            kappaLow[i].setConstant(True)
            kappaHigh[i] = RooRealVar("kappaH_%s_%s"%(i+1,codename),"kappaH_%s_%s"%(i+1,codename),norm_sig_sm_up[i]/norm_sig_sm)
            kappaHigh[i].setConstant(True)
            kappaLow_sum_d = kappaLow_sum_d*norm_sig_sm_down[i]/norm_sig_sm
            kappaHigh_sum_d = kappaHigh_sum_d*norm_sig_sm_up[i]/norm_sig_sm

            name_forCorr=signal_shapeSyst[i]
            name_forCorr=isItCorrelated_name(signal_shapeSyst[i])
            
            theWS.factory("%s[-7,7]"%name_forCorr)
            theta[i] = theWS.var("%s"%name_forCorr)
            
            aTGCPdf_norm[i] = AsymPow('ATGCPdf_anoCoupl_process_%s_integral%s'%(codename,i+1),
                                      'ATGCPdf_proc_%s_integral%s'%(codename,i+1),
                                      kappaLow[i],
                                      kappaHigh[i],
                                      theta[i])

        if (len(signal_shapeSyst)==1):
            aTGCPdf_norm_sum = aTGCPdf_norm[0]
        else:
            for i in range(0,len(signal_shapeSyst)):
                if (i==0): prodset=RooArgList(aTGCPdf_norm[i])
                else: prodset.add(RooArgList(aTGCPdf_norm[i]))
            aTGCPdf_norm_sum = RooProduct("aTGCPdf_norm_sum","aTGCPdf_norm_sum",prodset)

        kappaLow_sum = RooRealVar("kappaLow_sum","kappaLow_sum",kappaLow_sum_d)
        kappaHigh_sum = RooRealVar("kappaHigh_sum","kappaHigh_sum",kappaHigh_sum_d)

        aTGCPdf_norm_sum.SetNameTitle('ATGCPdf_anoCoupl_process_%s_norm'%codename,
                                      'ATGCPdf_proc_%s_norm'%codename)
        

    if (Ndim==1):
        aTGCPdf = RooACSemiAnalyticPdf_1D('ATGCPdf_anoCoupl_process_%s'%codename,
                                          'ATGCPdf_proc_%s'%codename,
                                          wpt,
                                          par1,
                                          dibosonPdf,
                                          '%s/signal_proc_%s.root'%(basepath,codename),
                                          limtype
                                          )
    elif (Ndim==2):
        aTGCPdf = RooACSemiAnalyticPdf_2D('ATGCPdf_anoCoupl_process_%s'%codename,
                                          'ATGCPdf_proc_%s'%codename,
                                          wpt,
                                          par1,
                                          par2,                                 
                                          dibosonPdf,
                                          '%s/signal_proc_%s.root'%(basepath,codename),
                                          limtype
                                          )
    elif (Ndim==3):
        aTGCPdf = RooACSemiAnalyticPdf_3D('ATGCPdf_anoCoupl_process_%s'%codename,
                                          'ATGCPdf_proc_%s'%codename,
                                          wpt,
                                          par1,
                                          par2,                                 
                                          par3,
                                          dibosonPdf,
                                          '%s/signal_proc_%s.root'%(basepath,codename),
                                          limtype
                                          )
    else:
        raise RuntimeError('Ndim = %s'%Ndim,
                           'We can only use 1D (par1_TH1, par1_TF1)  2D (par1par2_TF2par, par1par2_TH2, par1par2_TF2) or 3D (par1par2par3_TH3, par1par2par3_TF3) models right now!')

    if (doSignalShape_unc):
        aTGCPdf_up = {}
        aTGCPdf_down = {}
        for i in range(0,len(signal_shapeSyst)):

            name_forCorr=signal_shapeSyst[i]
            name_forCorr=isItCorrelated_name(signal_shapeSyst[i])

            if (Ndim==1):
                aTGCPdf_up[i] = RooACSemiAnalyticPdf_1D('ATGCPdf_anoCoupl_process_%s_%sUp'%(codename,name_forCorr),
                                                        'ATGCPdf_proc_%s'%codename,
                                                        wpt,
                                                        par1,
                                                        dibosonPdf_up[i],
                                                        '%s/signal_proc_%s.root'%(basepath,codename),
                                                        limtype
                                                        )
            elif (Ndim==2):
                aTGCPdf_up[i] = RooACSemiAnalyticPdf_2D('ATGCPdf_anoCoupl_process_%s_%sUp'%(codename,name_forCorr),
                                                        'ATGCPdf_proc_%s'%codename,
                                                        wpt,
                                                        par1,
                                                        par2,                                 
                                                        dibosonPdf_up[i],
                                                        '%s/signal_proc_%s.root'%(basepath,codename),
                                                        limtype
                                                        )
            elif (Ndim==3):
                aTGCPdf_up[i] = RooACSemiAnalyticPdf_3D('ATGCPdf_anoCoupl_process_%s_%sUp'%(codename,name_forCorr),
                                                        'ATGCPdf_proc_%s'%codename,
                                                        wpt,
                                                        par1,
                                                        par2,                        
                                                        par3,
                                                        dibosonPdf_up[i],
                                                        '%s/signal_proc_%s.root'%(basepath,codename),
                                                        limtype
                                                        )
            else:
                raise RuntimeError('Ndim = %s'%Ndim,
                                   'We can only use 2D (par1par2_TF2par, par1par2_TH2, par1par2_TF2) or 3D (par1par2par3_TH3, par1par2par3_TF3) models right now!')
            

            if (Ndim==1):
                aTGCPdf_down[i] = RooACSemiAnalyticPdf_1D('ATGCPdf_anoCoupl_process_%s_%sDown'%(codename,name_forCorr),
                                                          'ATGCPdf_proc_%s'%codename,
                                                          wpt,
                                                          par1,
                                                          dibosonPdf_down[i],
                                                          '%s/signal_proc_%s.root'%(basepath,codename),
                                                          limtype
                                                          )
            elif (Ndim==2):
                aTGCPdf_down[i] = RooACSemiAnalyticPdf_2D('ATGCPdf_anoCoupl_process_%s_%sDown'%(codename,name_forCorr),
                                                          'ATGCPdf_proc_%s'%codename,
                                                          wpt,
                                                          par1,
                                                          par2,                                 
                                                          dibosonPdf_down[i],
                                                          '%s/signal_proc_%s.root'%(basepath,codename),
                                                          limtype
                                                          )
            elif (Ndim==3):
                aTGCPdf_down[i] = RooACSemiAnalyticPdf_3D('ATGCPdf_anoCoupl_process_%s_%sDown'%(codename,name_forCorr),
                                                          'ATGCPdf_proc_%s'%codename,
                                                          wpt,
                                                          par1,
                                                          par2,
                                                          par3,
                                                          dibosonPdf_down[i],
                                                          '%s/signal_proc_%s.root'%(basepath,codename),
                                                          limtype
                                                          )
            else:
                raise RuntimeError('Ndim = %s'%Ndim,
                                   'We can only use 2D (par1par2_TF2par, par1par2_TH2, par1par2_TF2) or 3D (par1par2par3_TH3, par1par2par3_TF3) models right now!' )
  
    
    getattr(theWS, 'import')(data)
    for i in range(0,Nbkg_int):
        getattr(theWS, 'import')(bkgHist[i])
    for j in range(0,Nbkg_int):
        for i in range(0,len(background_shapeSyst[j])):
            print '\t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adding bkgUp in ws: ',bkgHist_systUp[j][i]
 
            getattr(theWS, 'import')(bkgHist_systUp[j][i])
            getattr(theWS, 'import')(bkgHist_systDown[j][i])

    getattr(theWS, 'import')(aTGCPdf)
    if (doSignalShape_unc):
        for i in range(0,len(signal_shapeSyst)):
            getattr(theWS, 'import')(aTGCPdf_up[i])
            getattr(theWS, 'import')(aTGCPdf_down[i])
#            getattr(theWS, 'import')(aTGCPdf_norm[i])
        getattr(theWS, 'import')(aTGCPdf_norm_sum)
    theWS.Print()
    
    fout = TFile('%s_ws.root'%(codename), 'recreate')
    theWS.Write()
    fout.Close()

### make the card for this channel and plane ID
    card = """
# Simple counting experiment, with one signal and a few background processes 
imax 1  number of channels
jmax {Nbkg_int}  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,Nbkg_int=Nbkg_int)
    for i in range(0,Nbkg_int):
        card += """
shapes anomalousCoupling_bkg{Nbkg_int}_{codename}  {codename} ./{codename}_ws.root proc_{codename}:$PROCESS proc_{codename}:$PROCESS_$SYSTEMATIC""".format(Nbkg_int=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
    card += """
shapes data_obs                {codename} ./{codename}_ws.root proc_{codename}:$PROCESS """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,Nbkg_int=Nbkg_int)
    if (doSignalShape_unc):
        card += """   
shapes anoCoupl_process_{codename} {codename} ./{codename}_ws.root proc_{codename}:ATGCPdf_$PROCESS proc_{codename}:ATGCPdf_$PROCESS_$SYSTEMATIC """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
    else:
        card += """   
shapes anoCoupl_process_{codename} {codename} ./{codename}_ws.root proc_{codename}:ATGCPdf_$PROCESS
""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        
    card += """   
------------
bin {codename} 
observation {norm_obs}
------------
bin                         {codename}\t\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """\t\t\t{codename}""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    card += """       
process\t\t\t    anoCoupl_process_{codename}    """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """\tanomalousCoupling_bkg{Nbkg_int}_{codename}""".format(Nbkg_int=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    card += """       
process                     0	  	      """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """ \t\t\t\t{i}""".format(i=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)
        
    card += """       
rate                        {norm_sig_sm}\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """ \t\t\t{norm_bkg}""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

	
    card += """           
------------
"""
    for i in range(0,NlnN):

# if found signal or codename in the list of names that are affected by lnN unc:
        if (('%s_signal'%codename in lnN_for[i]) or (any(codename in s for s in lnN_for[i]))):
            card+="""
{lnN_name}         lnN """.format(lnN_name=lnN_name[i])
            if ('%s_signal'%codename in lnN_for[i]):
                # if lnN syst affects signal:
                index_s=lnN_for[i].index('%s_signal'%codename)
                card+=""" {lnN_value}      """.format(lnN_value=lnN_value[i][index_s])
            else:
                card+=""" -      """
            for j in range(0,Nbkg_int):
                name_for_lnN=codename
                name_for_lnN+='_'
                name_for_lnN+=bkg_name[j]
                if (name_for_lnN in lnN_for[i]):
                    index=lnN_for[i].index(name_for_lnN)
                    card+="""\t\t\t\t{lnN_value}""".format(lnN_value=lnN_value[i][index])
                else:
                    card+="""\t\t\t\t-"""

################ bkg shape syst:

    for j in range(0,Nbkg_int):
        for i in range(0,len(background_shapeSyst[j])):
            # write out only those bkg shape unc that are not correlated with signal shape syst:
            if not(isItCorrelated(background_shapeSyst[j][i])):
                card += """
{background_shapeSyst} shape1  """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])
                for k in range(0,j+1):
                    card += """-\t\t\t\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])
                card += """1.0 """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])
                for k in range(1,Nbkg_int-j):
                    card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])

    if (doSignalShape_unc):
        for i in range(0,len(signal_shapeSyst)):
            name_forCorr=signal_shapeSyst[i]
            if (isItCorrelated(signal_shapeSyst[i])):
                name_forCorr=isItCorrelated_name(signal_shapeSyst[i])
            card += """
{signal_shapeSyst}        shape1  1.0          """.format(signal_shapeSyst=name_forCorr)
    
            for j in range(0,Nbkg_int):
                if (isItCorrelated(signal_shapeSyst[i])):
                    isitcorr=false
                    for k in range(0,len(background_shapeSyst[j])):
                        if (isItCorrelated(background_shapeSyst[j][k])):
                            isitcorr=true
                    if (isitcorr):
                        card += """\t\t\t\t1.0""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
                    if not(isitcorr):
                        card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
                else:
                    card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    print card

    cardfile = open('aC_%s.txt'%(codename),'w')
    cardfile.write(card)
    cardfile.close
