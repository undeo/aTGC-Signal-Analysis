from ROOT import *
import os
from optparse import OptionParser

parser	= OptionParser()
parser.add_option('-P', '--POI', dest='POI', default='cwww')
parser.add_option('--pos', action='store_true', default=False)
parser.add_option('--neg', action='store_true', default=False)
(options,args) = parser.parse_args()


if options.neg:
	startvals	= {'cwww' : -8.28, 'ccw' : -11.4, 'cb' : -53}
	endvals		= {'cwww' : -8.58, 'ccw' : -11.8, 'cb' : -60}
	intvals		= {'cwww' : -0.01, 'ccw' : -0.01, 'cb' : -0.1}
elif options.pos:
	startvals	= {'cwww' : 8.4, 'ccw' : 10.5, 'cb' : 51.5}
	endvals		= {'cwww' : 8.9, 'ccw' : 10.9, 'cb' : 62.5}
	intvals		= {'cwww' : 0.01, 'ccw' : 0.01, 'cb' : 0.1}
else:
	raise RuntimeError('run with either --pos or --neg!')
iters	= {'cwww' : abs(int((endvals['cwww']-startvals['cwww'])/intvals['cwww'])), 'ccw' : abs(int((endvals['ccw']-startvals['ccw'])/intvals['ccw'])), 'cb' : abs(int((endvals['cb']-startvals['cb'])/intvals['cb']))}


POI	= options.POI 
if POI=='cwww':
	noPOI	= ['ccw','cb']
	vals	= {'cwww' : startvals[POI], 'ccw' : 0,'cb' : 0}
elif POI=='ccw':
	noPOI	= ['cwww','cb']
	vals	= {'cwww' : 0,'ccw' : startvals[POI], 'cb' : 0}
elif POI=='cb':
	noPOI	= ['cwww','ccw']
	vals	= {'cwww': 0,'ccw' : 0, 'cb' : startvals[POI]}



#cwww
for i in range(iters[POI]):
	vals[POI] = startvals[POI] + intvals[POI]*i
	comstring = 'combine workspace_WWWZ_linter1.root -M MultiDimFit --floatOtherPOIs=0 --algo grid -t -1 --expectSignal=1 --points 3 --freezeNuisances %s,%s --setPhysicsModelParameters cwww=%s,ccw=%s,cb=%s -n _sens%s -P %s --setPhysicsModelParameterRange %s=-2,2'%(noPOI[0],noPOI[1],vals['cwww'],vals['ccw'],vals['cb'],POI,POI,POI)
	print comstring
	os.system(comstring)

	fileIn = TFile.Open('higgsCombine_sens%s.MultiDimFit.mH120.root'%POI)
	tree = fileIn.Get('limit')
	for j in range(3):
		tree.GetEntry(1+j)
		if 2*tree.deltaNLL>3.84:
			print "stop, %s - %s"%(startvals[POI]+intvals[POI]*(i-1),startvals[POI]+intvals[POI]*i)
			fileIn.Close()
			exit(0)
	fileIn.Close()



