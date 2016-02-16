from ROOT import *
from array import array



def make_tree(ch = "ele", iscwwwl=1, iscwl=0, iscbl=0):
	
	fileIn	= TFile.Open("Input/WW-aTGC-%s.root"%ch)
	treeIn	= fileIn.Get("treeDumper/BasicTree")
	fileOut = TFile("Output/WW-aTGC-%s-inputfile.root"%ch, "recreate")
	treeIn.SetBranchStatus("*",0);
	treeOut = treeIn.CloneTree(0);
	treeIn.SetBranchStatus("aTGCWeights",1);
	treeIn.SetBranchStatus("m_lvj",1);
	treeIn.SetBranchStatus("PUweight",1);

	weight		= array('f',[1])
	branch_weight	= treeOut.Branch("weight",weight,"weight")
	m_lvj		= array('f',[1])
	branch_m_lvj	= treeOut.Branch("m_lvj",m_lvj,"m_lvj")
	c_wwwl		= array('f',[1])
	c_wl		= array('f',[1])	
	c_bl		= array('f',[1])
	if iscwwwl:
		branch_c_wwwl	= treeOut.Branch("c_wwwl_grid",c_wwwl,"c_wwwl")
	if iscwl:
		branch_c_wl	= treeOut.Branch("c_wl_grid",c_wl,"c_wl")	
	if iscbl:
		branch_c_bl	= treeOut.Branch("c_bl_grid",c_bl,"c_bl")


	lumi_tmp 	= 2093.917403402
	NEntries	= treeIn.GetEntries()
	for i in range(NEntries):
		if i%10000==0:
			print i
		treeIn.GetEntry(i)
		m_lvj[0]	= treeIn.m_lvj
		weight_part	= 1/20. * lumi_tmp * treeIn.PUweight

		c_wwwl[0]	= 0
		c_wl[0]		= 0
		c_bl[0]		= 0
		weight[0]	= treeIn.aTGCWeights[7] * weight_part
		treeOut.Fill()

		if iscwwwl:
			c_wwwl[0]	= 12
			c_wl[0]		= 0
			c_bl[0]		= 0
			weight[0]	= treeIn.aTGCWeights[0] * weight_part
			treeOut.Fill()

			c_wwwl[0]	= -12
			c_wl[0]		= 0
			c_bl[0]		= 0
			weight[0]	= treeIn.aTGCWeights[1] * weight_part
			treeOut.Fill()
		if iscwl:
			c_wwwl[0]	= 0
			c_wl[0]		= 20
			c_bl[0]		= 0
			weight[0]	= treeIn.aTGCWeights[2] * weight_part
			treeOut.Fill()

			c_wwwl[0]	= 0
			c_wl[0]		= -20
			c_bl[0]		= 0
			weight[0]	= treeIn.aTGCWeights[3] * weight_part
			treeOut.Fill()
		if iscbl:
			c_wwwl[0]	= 0
			c_wl[0]		= 0
			c_bl[0]		= 60
			weight[0]	= treeIn.aTGCWeights[4] * weight_part
			treeOut.Fill()

			c_wwwl[0]	= 0
			c_wl[0]		= 0
			c_bl[0]		= -60
			weight[0]	= treeIn.aTGCWeights[5] * weight_part
			treeOut.Fill()

		
	print NEntries
	print "--> " + str(NEntries) + "*" + str(treeOut.GetEntries()/NEntries) + " = " + str(treeOut.GetEntries()) + " total Entries" 
	print "Write to file " + fileOut.GetName()
	treeOut.Write()
	fileOut.Close()


make_tree("ele",1,0,0)
make_tree("mu",1,0,0)
