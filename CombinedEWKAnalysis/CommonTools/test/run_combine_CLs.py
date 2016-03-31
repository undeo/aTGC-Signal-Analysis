import os

negative = False

start1 = 8.5
end1   = 10.5

start2 = -9
end2 = -8

prec = 0.05

name = 'CLS'
outname = 'CLSLHC_limit.root'
outnames = ''
poi = 'cwww'
filename = 'CLSLHC_upperlimit'

i = start1
n = 0
while i < end1:
	n += 1
	print 'combine workspace_WWWZ.root --verbose=-1 -M  HybridNew -n %s --clsAcc 0 -i 1 -T 100 --frequentist --testStat LHC --importantContours=0.95 --fork 4 --redefineSignalPOI %s --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --freezeNuisances cwww,ccw,cb --seed -1 --singlePoint %s --saveToys --saveHybridResult --expectSignal=1 -t -1'%(name,poi,i)
	os.system('combine workspace_WWWZ.root --verbose=-1 -M HybridNew -n %s --clsAcc 0 -i 1 -T 100 --frequentist --testStat LHC --importantContours=0.95 --fork 4 --redefineSignalPOI %s --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --freezeNuisances cwww,ccw,cb --seed -1 --singlePoint %s --saveToys --saveHybridResult --expectSignal=1 -t -1'%(name,poi,i))
	i += prec
	if n%2 == 0:
		print '###############################################################################################'
		print 'hadd -f -k file_2tmp%s.root higgsCombine%s*'%(n,name)
		os.system('hadd -f -k file_2tmp%s.root higgsCombine%s*'%(n,name))
		print 'rm higgsCombine%s*'%name
		os.system('rm higgsCombine%s*'%name)
		print '###############################################################################################'
	if n%4 == 0:
		print '###############################################################################################'
		print 'hadd -f -k file_4tmp%s.root file_2tmp*.root'%n
		os.system('hadd -f -k file_4tmp%s.root file_2tmp*.root'%n)
		print 'rm file_2tmp*.root'
		os.system('rm file_2tmp*.root')
		print '###############################################################################################'
	if n%8 == 0 or n == end1:
		outnames += '%s%s '%(n,outname)
		print '###############################################################################################'
		print 'hadd -f -k %s%s file_4tmp*.root'%(n,outname)
		os.system('hadd -f -k %s%s file_4tmp*.root'%(n,outname))
		print 'rm higgsCombine%s*'%name
		os.system('rm higgsCombine%s*'%name)
		print outnames
		print '###############################################################################################'
print '###############################################################################################'
print 'hadd -f %s_%s%s-%s.root %s'%(filename,poi,start1,end1,outnames)
os.system('hadd -f %s_%s%s-%s.root %s'%(filename,poi,start1,end1,outnames))
print 'mv %s tmp/'%outnames
os.system('mv %s tmp/'%outnames)
print '###############################################################################################'





if negative:
	outnames = ''
	i = start2
	n = 0
	while i < end2:
		n += 1
		print 'combine workspace_WWWZ.root --verbose=-1 -M  HybridNew -n %s --frequentist --clsAcc 0 --testStat PL --importantContours=0.95 --fork 4 --redefineSignalPOI %s --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --freezeNuisances cwww,ccw,cb --seed -1 --singlePoint %s --saveToys --saveHybridResult --expectSignal=1 -t -1'%(name,poi,i)
		os.system('combine workspace_WWWZ.root --verbose=-1 -M HybridNew -n %s --frequentist --clsAcc 0 --testStat PL --importantContours=0.95 --fork 4 --redefineSignalPOI %s --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --freezeNuisances cwww,ccw,cb --seed -1 --singlePoint %s --saveToys --saveHybridResult --expectSignal=1 -t -1'%(name,poi,i))
		i += prec
		if n %5 == 0:
			outnames += '%s%s '%(n,outname)
			print '###############################################################################################'
			print 'hadd -f %s%s higgsCombine%s*'%(n,outname,name)
			os.system('hadd -f %s%s higgsCombine%s*'%(n,outname,name))
			print 'rm higgsCombine%s*'%name
			os.system('rm higgsCombine%s*'%name)
			print outnames
			print '###############################################################################################'
	print '###############################################################################################'
	print 'hadd -f CLS_lowerlimit_%s%s-%s.root %s'%(poi,outnames,start2,end2)
	os.system('hadd -f CLS_lowerlimit_%s%s-%s.root %s'%(poi,outnames,start2,end2))
	print 'rm %s'%outnames
	os.system('rm %s'%outnames)
	print '###############################################################################################'


	#print '###############################################################################################'
	#print 'hadd -f CLS_limit_%s.root CLS_*limit_%s.root'%(poi,poi)
	#os.system('hadd -f CLS_limit_%s.root CLS_*limit_%s.root'%(poi,poi))
	#print '###############################################################################################'
