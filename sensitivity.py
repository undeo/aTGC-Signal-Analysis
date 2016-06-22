from ROOT import *
import os
from optparse import OptionParser

parser	= OptionParser()
parser.add_option('-P', dest='POI', default='cwww')
(options,args) = parser.parse_args()

startvalcwww = 10.6
endvalcwww = 10.7
intvalcwww = 0.01
iterscwww = abs(int((endvalcwww-startvalcwww)/intvalcwww))
startvalccw = 13.68
endvalccw = 13.8
intvalccw = 0.01
itersccw = abs(int((endvalccw-startvalccw)/intvalccw))
startvalcb = 60.
endvalcb = 63
intvalcb = 0.1
iterscb = abs(int((endvalcb-startvalcb)/intvalcb))

cwww=False
ccw=False
cb=False
if options.POI=='cwww':
	cwww=True
elif options.POI=='ccw':
	ccw=True
elif options.POI=='cb':
	cb=True

#cwww
if cwww:
	for i in range(iterscwww):
		comstring = 'combine workspace_WWWZ_linter1.root -M MultiDimFit --floatOtherPOIs=0 --algo grid -t -1 --expectSignal=1 --points 3 --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=%s,ccw=0,cb=0 -n _senscwww -P cwww --setPhysicsModelParameterRange cwww=-2,2'%(startvalcwww+intvalcwww*i)
		print comstring
		os.system(comstring)

		fileIn = TFile.Open('higgsCombine_senscwww.MultiDimFit.mH120.root')
		tree = fileIn.Get('limit')
		for j in range(3):
			tree.GetEntry(1+j)
			if 2*tree.deltaNLL>3.84:
				print "stop, %s - %s"%(startvalcwww+intvalcwww*(i-1),startvalcwww+intvalcwww*i)
				fileIn.Close()
				exit(0)
		fileIn.Close()


#ccw
if ccw:	
	for i in range(itersccw):
		comstring = 'combine workspace_WWWZ_linter2.root -M MultiDimFit --floatOtherPOIs=0 --algo grid -t -1 --expectSignal=1 --points 3 --freezeNuisances cwww,cb --setPhysicsModelParameters cwww=0,ccw=%s,cb=0 -n _sensccw -P ccw --setPhysicsModelParameterRange ccw=-2,2'%(startvalccw+intvalccw*i)
		print comstring
		os.system(comstring)

		fileIn = TFile.Open('higgsCombine_sensccw.MultiDimFit.mH120.root')
		tree = fileIn.Get('limit')
		for j in range(3):
			tree.GetEntry(1+j)
			if 2*tree.deltaNLL>3.84:
				print "stop, %s - %s"%(startvalccw+intvalccw*(i-1),startvalccw+intvalccw*i)
				fileIn.Close()
				exit(0)
		fileIn.Close()

#cb
if cb:	
	for i in range(iterscb):
		comstring = 'combine workspace_WWWZ_linter3.root -M MultiDimFit --floatOtherPOIs=0 --algo grid -t -1 --expectSignal=1 --points 3 --freezeNuisances cwww,ccw --setPhysicsModelParameters cwww=0,ccw=0,cb=%s -n _senscb -P cb --setPhysicsModelParameterRange cb=-2,2'%(startvalcb+intvalcb*i)
		print comstring
		os.system(comstring)

		fileIn = TFile.Open('higgsCombine_senscb.MultiDimFit.mH120.root')
		tree = fileIn.Get('limit')
		for j in range(3):
			tree.GetEntry(1+j)
			if 2*tree.deltaNLL>3.84:
				print "stop, %s - %s"%(startvalcb+intvalcb*(i-1),startvalcb+intvalcb*i)
				fileIn.Close()
				exit(0)
		fileIn.Close()
