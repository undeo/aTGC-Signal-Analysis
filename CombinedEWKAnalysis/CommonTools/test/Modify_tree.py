from ROOT import *
from array import array
import math



def merge_trees(list_):    

  mergedtree = TChain("BasicTree")
 
  for file_ in list_:    
    mergedtree.Add(file_)

  return mergedtree
  

def prepare_trees_with_cuts(filenameOut = "Output.root", list_of_trees = [""], is_data = False, is_ttbar=False, is_WJets=False):
  
  if is_data:
    oldtree = TChain("treeDumper/BasicTree")
    for file_ in list_of_trees:    
      oldtree.Add(file_)
  else:
    oldtree = merge_trees(list_of_trees)
  
  
  nEntries = oldtree.GetEntries()
  newtree = oldtree.CloneTree(0)  
  
  ## create new branches
  var_puweight = array('f',[1])
  branch_puweight = newtree.Branch("puweight", var_puweight, "puweight")
  var_genweight = array('f',[1])
  branch_genweight = newtree.Branch("genweight", var_genweight, "genweight")
  var_lumiweight = array('f',[1])
  branch_lumiweight = newtree.Branch("lumiweight", var_lumiweight, "lumiweight")
  var_issignal = array('f',[1])
  branch_issignal = newtree.Branch("issignal", var_issignal, "issignal")
  var_MWW = array('f',[1])
  branch_MWW = newtree.Branch("MWW", var_MWW, "MWW")
  var_Mjpruned = array('f',[1]) 
  branch_Mjpruned = newtree.Branch("Mjpruned", var_Mjpruned, "Mjpruned")
  if is_data:
    var_weight = array('f',[1])
    branch_weight = newtree.Branch("weight", var_weight, "weight")
  
  ##fill new branches  
  used_events = 0
  print list_of_trees
  print "Number of Entries: " + str(nEntries) + " , " "Number of used events: "
  for i in range(nEntries):    
    use_event = False
    oldtree.GetEntry(i)
    ##apply cuts
    ##additional cuts for ttbar and wjets
#    if is_ttbar:
#	if oldtree.nbtag>0:
#	    event = True
#	else:
#	    event = False
#   elif is_WJets:
#	if oldtree.jet_mass_pruned < 65 or oldtree.jet_mass_pruned > 95 and oldtree.nbtag == 0:
#	    event = True
#	else:
#	    event = False
#    else:
#        event = True
    event = True
    #other cuts
    if oldtree.jet_pt>200 and oldtree.jet_tau2tau1<0.5 and oldtree.jet_mass_pruned<150 and oldtree.jet_mass_pruned>40 and oldtree.W_pt>200 and abs(oldtree.deltaR_LeptonWJet)>math.pi/2 and abs(oldtree.deltaPhi_WJetMet)>2 and abs(oldtree.deltaPhi_WJetWlep)>2 and event:
      use_event = True
      used_events += 1
      if used_events%1000 == 0:
	print str(used_events) + " / " + str(i) + " / " + str(nEntries)
      
    if use_event==True:      
      var_issignal[0] = 1
      var_MWW[0] = oldtree.m_lvj
      var_Mjpruned[0] = oldtree.jet_mass_pruned
      if is_data:
	var_weight[0] = 1
	var_lumiweight[0] = 1
	var_puweight[0] = 1
	var_genweight[0] = 1
      else:
	var_lumiweight[0] = oldtree.weight
	var_puweight[0] = oldtree.PUweight
	var_genweight[0] = oldtree.genWeight      
      newtree.Fill()

  print str(used_events) + " / " + str(nEntries)
  print "____"
  
  ##set name to 'tree'
  newtree.SetName("tree")
  fileOut = TFile(filenameOut, "recreate")
  newtree.Write()
  fileOut.Close()


    

chan = ["ele", "mu"]
for ch in chan:
  ###merge and rename other trees
  data_list_of_trees = ["data_Prompt_%s.root"%(ch), "data_05Oct_%s.root"%(ch)]
  prepare_trees_with_cuts("treeEDBR_data_xww_%s.root"%(ch), data_list_of_trees, True, False, False)
  #W_list_of_trees = ["WJets_Ht100To200_%s.root"%(ch), "WJets_Ht200To400_%s.root"%(ch), "WJets_Ht400To600_%s.root"%(ch), "WJets_Ht600ToInf_%s.root"%(ch)]
  #prepare_trees_with_cuts("treeEDBR_WJets_xww_%s.root"%(ch), W_list_of_trees, False, False, True)
  #T_list_of_trees = ["s-ch_%s.root"%(ch), "t-ch_antitop_%s.root"%(ch), "t-ch_top_%s.root"%(ch), "tW-ch_%s.root"%(ch)]
  #prepare_trees_with_cuts("treeEDBR_SingleTop_xww_%s.root"%(ch), T_list_of_trees, False, False, False)
  #VV_list_of_trees = ["WW_%s.root"%(ch), "WZ_%s.root"%(ch)]
  #prepare_trees_with_cuts("treeEDBR_VV_xww_%s.root"%(ch), VV_list_of_trees, False, False, False)
  #TTBAR_list_of_trees = ["ttbar_%s.root"%(ch)]
  #prepare_trees_with_cuts("treeEDBR_TTBARpowheg_xww_%s.root"%(ch), TTBAR_list_of_trees, False, True, False)
  






