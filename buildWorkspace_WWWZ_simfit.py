import pyroot_logon
import limits
import os
import sys

from array import *

from ROOT import *
from optparse import OptionParser
from ConfigParser import SafeConfigParser
#@#import new pdfs
gSystem.Load("/afs/cern.ch/work/c/crenner/CMSSW_7_1_5/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so")
#gSystem.Load('ATGCRooStatsTMP/PDFs/PdfDiagonalizer_cc.so')
#gSystem.Load('ATGCRooStatsTMP/PDFs/Util_cxx.so')
#gSystem.Load('ATGCRooStatsTMP/PDFs/hyperg_2F1_c.so')
#gSystem.Load('ATGCRooStatsTMP/PDFs/HWWLVJRooPdfs_cxx.so')


from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf, RooErfExpDecoPdf
#@#

parser = OptionParser(description="%prog : A RooStats Implementation of Anomalous Triple Gauge Coupling Analysis.",
                      usage="buildWZworkspace --config=example_config.cfg")
cfgparse = SafeConfigParser()

parser.add_option("--config",dest="config",help="The name of the input configuration file.")
parser.add_option("-o","--outputcard",dest="ocard",help="Name of the generated combined datacard")
parser.add_option("--binned",action='store_true',dest="binned",default=False,help="use RooDataHist for data -> binned fit for combine") 
(options,args) = parser.parse_args()



if options.config is None:
    print 'Need to specify --config'
    exit(1)
        
if options.binned:
	raw_input("Will result in a binned fit! Press any key to continue")
cfgparse.read(options.config)
options.config = cfgparse # put the parsed config file into our options

cfg = options.config

norm_sig_sm = -1
norm_sig_sm_up = -1
norm_sig_sm_down = -1
norm_bkg = -1
norm_obs = -1
cardfilenames = []

fit_sections = cfg.sections()

signalModel ='-1'
signalModel = cfg.get('Global','model')     
print 'modelName= ',signalModel

fit_sections.remove('Global') #don't need to iterate over the global configuration

limtype = -1
Ndim = -1

if ( signalModel == 'par1par2par3_TH3'):
    limtype = 0
    Ndim = 3
elif ( signalModel == 'par1par2par3_TF3'):
    limtype = 1
    Ndim = 3
else:
    raise RuntimeError('InvalidCouplingChoice %s'%signalModel,
                       'We can only use 3D (par1par2par3_TH3, par1par2par3_TF3) models right now!')

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

par1low  = float(cfg.get('Global', 'par1Low'))
par1high = float(cfg.get('Global', 'par1High'))
par2low  = float(cfg.get('Global', 'par2Low'))
par2high = float(cfg.get('Global', 'par2High'))
par3low  = float(cfg.get('Global', 'par3Low'))
par3high = float(cfg.get('Global', 'par3High'))


NSigBkg_corr_unc_int=0

basepath = '%s/src/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'%os.environ['CMSSW_BASE']

for section in fit_sections:
  
    codename = section
    #@#
    if 'el' in codename:
      ch = 'el'
    elif 'mu' in codename:
      ch = 'mu'
    else:
      raise RuntimeError('wrong channel name: ' + section)
    if 'WW' in codename:
      cat = 'WW'
    elif 'WZ' in codename:
      cat = 'WZ'
    #@#




    f = TFile('%s/%s.root'%(basepath,codename))

    #@# get old workspace 
    oldWS = f.Get("w")
    old_obs = oldWS.var("rrv_mass_lvj")
    old_obs.setVal(2000)
    old_obs.setRange(900,3500)
    #old_obs.SetName("m_lvj")
    old_obs2 = oldWS.var('rrv_mass_j')
    #old_obs2.SetName('m_j')
    oldWS.Print()

    theWS = RooWorkspace('proc_%s'%codename, 'proc_%s'%codename)

    Nbkg = cfg.get(codename,'Nbkg')
    print "Nbkg= ",Nbkg
    Nbkg_int=int(Nbkg)

    bkg_name 	= []
    for i in range(1,Nbkg_int+1):
        bkg_name.append(cfg.get(codename,'bkg%i_name'%i))

    background 	= []
    bkg_norms	= []
    norm_bkg = []
    bkgFits = {}
    oldWS.allFunctions().Print("V")


    for i in range(0,Nbkg_int):
        background.append(oldWS.pdf(bkg_name[i]))
        bkg_norms.append(oldWS.function(bkg_name[i]+'_norm'))
        norm_bkg.append(1)
        bkgFits[i] = oldWS.pdf(bkg_name[i])
        bkg_norms[i].Print()


    aTGCPdf = oldWS.pdf("ATGCPdf_%s"%codename)
    normfactor_atgc	= oldWS.function("ATGCPdf_%s_norm"%codename)



    #get normalisation of SM VV
    #norm_sig_sm = oldWS.var("rate_VV").getVal()
    norm_sig_sm	= 1


    getattr(theWS,'import')(normfactor_atgc,RooFit.RecycleConflictNodes())
    #normfactor_atgc	= oldWS.function("normfactor_3d_%s"%channel)
    #normfactor_atgc.SetName("ATGCPdf_%s_norm"%codename)
    #normfactor_atgc.SetName("normfactor_3d_%s"%codename)
    #normfactor.Print()
    #raw_input('.')

    cfg_items=cfg.items(codename)

    #@#


    par1 = theWS.factory('%s[0., %f, %f]' % (par1name, par1low, par1high))
    par2 = theWS.factory('%s[0., %f, %f]' % (par2name, par2low, par2high))
    par3 = theWS.factory('%s[0., %f, %f]' % (par3name, par3low, par2high))
    
    wpt = old_obs
    wpt2 = old_obs2

    getattr(theWS,'import')(wpt)
    getattr(theWS,'import')(wpt2)

    data	= oldWS.data('dataset_mj_%s'%codename)
    data.SetName('data_obs')
    data.Print()

	
    getattr(theWS, 'import')(data)
    for i in range(0,Nbkg_int):
        getattr(theWS, 'import')(bkgFits[i],RooFit.RecycleConflictNodes())
        getattr(theWS, 'import')(bkg_norms[i],RooFit.RecycleConflictNodes())
        theWS.allVars().Print("V")
    getattr(theWS, 'import')(aTGCPdf)
    theWS.allFunctions().Print("V")
   
    theWS.Print()
    
    fout = TFile('%s_ws.root'%(codename), 'recreate')
    theWS.Write()
    fout.Close()
    f.Close()
### make the card for this channel and plane ID
    card = """
imax 1  number of channels
jmax {Nbkg_int}  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,Nbkg_int=Nbkg_int)
    for i in range(0,Nbkg_int):
        card += """
shapes {bkg_name}\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename,bkg_name=bkg_name[i])
    card += """
shapes data_obs\t\t\t\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename)    
    card += """   
shapes ATGCPdf_{codename}\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS
""".format(codename=codename)
        
    card += """   
------------
bin {codename} 
observation -1
------------
bin                         {codename}\t\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg)

    for i in range(0,Nbkg_int):
        card += """\t\t\t{codename}""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i])

    card += """       
process\t\t\t    ATGCPdf_{codename}    """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i])

    for i in range(0,Nbkg_int):
        card += """\t\t{bkg_name}""".format(Nbkg_int=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],bkg_name=bkg_name[i])

    card += """       
process                     0	  	      """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i])

    for i in range(0,Nbkg_int):
        card += """ \t\t\t\t{i}""".format(i=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i])
        
    card += """       
rate                        {norm_sig_sm}\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i])

    for i in range(0,Nbkg_int):
        card += """ \t\t\t{norm_bkg}""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i])

	
    card += """           
------------
"""
    for i in range(0,NlnN):

# if found signal or codename in the list of names that are affected by lnN unc:
        if (('%s_signal'%codename in lnN_for[i]) or (any(codename in s for s in lnN_for[i]))):
            card+="""
{lnN_name} \t\t\t lnN """.format(lnN_name=lnN_name[i])
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

# add pdf parameter uncertainties
    card += """	
Deco_WJets0_xww_sb_lo_from_fitting_{ch}_HP_mlvj_13TeV_eig0 param 0.0 1.4
Deco_WJets0_xww_sb_lo_from_fitting_{ch}_HP_mlvj_13TeV_eig1 param 0.0 1.4
Deco_WJets0_xww_sb_lo_from_fitting_{ch}_HP_mlvj_13TeV_eig2 param 0.0 1.4 
Deco_WJets0_xww_sim_{ch}_HP{cat}_mlvj_13TeV_eig0 param 0.0 1.4 
Deco_WJets0_xww_sim_{ch}_HP{cat}_mlvj_13TeV_eig1 param 0.0 1.4 
Deco_WJets0_xww_sim_{ch}_HP{cat}_mlvj_13TeV_eig2 param 0.0 1.4 
Deco_WJets0_xww_sim_{ch}_HP{cat}_mlvj_13TeV_eig3 param 0.0 1.4
Deco_TTbar_xww_signal_region_{ch}_HP{cat}_mlvj_13TeV_eig0 param 0.0 2.0
Deco_TTbar_xww_signal_region_{ch}_HP{cat}_mlvj_13TeV_eig1 param 0.0 2.0
slope_nuis param 1.0 0.05
rrv_c_ExpN_WJets0_xww_sb_lo_from_fitting_el flatParam
rrv_c_ExpN_WJets0_xww_sb_lo_from_fitting_mu flatParam
rrv_n_ExpN_WJets0_xww_sb_lo_from_fitting_el flatParam
rrv_n_ExpN_WJets0_xww_sb_lo_from_fitting_mu flatParam
normvar_WJets_el flatParam
normvar_WJets_mu flatParam""".format(ch=ch[:2],cat=cat[1])
    print card

    cardfile = open('aC_%s.txt'%(codename),'w')
    cardfile.write(card)
    cardfile.close()
    cardfilenames.append('aC_%s.txt'%(codename))

print ''
if options.binned:
	print 'data_obs added as RooDataHist -> combine will do a binned fit!'
else:
	print 'data_obs added as RooDataSet -> combine will do an unbinned fit!'
print ''

print '### combining Cards ###'
cmd = 'combineCards.py aC_WW_sig_el.txt aC_WZ_sig_el.txt aC_sb_lo_el.txt aC_sb_hi_el.txt aC_WW_sig_mu.txt aC_WZ_sig_mu.txt aC_sb_lo_mu.txt aC_sb_hi_mu.txt > %s'%(options.ocard)
print cmd
os.system(cmd)
print 'generated Card : %s'%options.ocard

