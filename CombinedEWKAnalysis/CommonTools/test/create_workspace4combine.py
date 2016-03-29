import os
from ROOT import *
from optparse import OptionParser


parser = OptionParser()
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='make new input trees')
parser.add_option('-p', '--plots', action='store_true', dest='do_plots', default=False, help='make plots')
(options,args) = parser.parse_args()


signalmodelcommand = 'python ATGCRooStatsTMP/make_PDF_input_all.py '
if options.newtrees:
	signalmodelcommand += '-n '
if options.do_plots:
	signalmodelcommand += '-p '

datacardname	= 'aC_WWWZ.txt'
model 		= 'CombinedEWKAnalysis.CommonTools.ACModel:par1par2par3_TF3_Model'
channels	= 'ch_WW_ele,ch_WW_mu,ch_WZ_ele,ch_WZ_mu'
wsout		= 'workspace_WWWZ.root'
r_cwww		= 'range_cwww=-20,20'
r_ccw		= 'range_ccw=-30,30'
r_cb		= 'range_cb=-70,70'
pois		= 'poi=cwww,ccw,cb'

print '''################################
### creating singal model WW ###
################################'''
print signalmodelcommand + '--cat WW'
os.system(signalmodelcommand + '--cat WW')
print '''################################
### creating singal model WZ ###
################################'''
print signalmodelcommand + '--cat WZ'
os.system(signalmodelcommand + '--cat WZ')

print '''################################################
### creating workspace for text2workspace.py ###
################################################'''
print 'python buildWorkspace_WWWZ.py --config=config_ATGC2 -o %s'%datacardname
os.system('python buildWorkspace_WWWZ.py --config=config_ATGC2 -o %s'%datacardname)

print '''################################
### creating final workspace ###
################################'''
print 'text2workspace.py %s -o %s -P %s --PO channels=%s --PO %s --PO %s --PO %s --PO %s'%(datacardname,wsout,model,channels,pois,r_cwww,r_ccw,r_cb)
os.system('text2workspace.py %s -o %s -P %s --PO channels=%s --PO %s --PO %s --PO %s --PO %s'%(datacardname,wsout,model,channels,pois,r_cwww,r_ccw,r_cb))
