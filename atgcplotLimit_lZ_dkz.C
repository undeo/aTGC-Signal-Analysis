#include <iostream>
#include <stdlib.h>
//#define __size_t unsigned // needed for glob.h!!
#include "/usr/lib/gcc/x86_64-redhat-linux/4.4.4/include/stddef.h"
#include <glob.h>
#include <assert.h>

#include "TH2D.h"
#include "TRegexp.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGraphPolar.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TStyle.h"

#include "atgcstyle.C"

using namespace std;

const float intlumifbinv = 2.3;
const int   beamcometev  = 13;


//======================================================================

TString par2latex(const TString& parname)
{
  if (parname.EqualTo("cwww") )  return "c_{WWW}/#Lambda^{2} (TeV^{-2})";
  if (parname.EqualTo("ccw") ) return "c_{W}/#Lambda^{2} (TeV^{-2})";
  if (parname.EqualTo("cb") ) return "c_{B}/#Lambda^{2} (TeV^{-2})";
  if (parname.EqualTo("dg1z") ) return "#Deltag_{1}^{Z}";
  if (parname.EqualTo("dkz") ) return "#Delta#kappa_{Z}";
  if (parname.EqualTo("lZ") ) return "#lambda_{Z}";

  return "UNKNOWN PAR "+parname;
}

//======================================================================

float parmin(const TString& parname)
{
  if (parname.EqualTo("cwww") )  return -20;
  if (parname.EqualTo("ccw") ) return -20;
  if (parname.EqualTo("cb") ) return -105;
  if (parname.EqualTo("dg1z") ) return -0.15;
  if (parname.EqualTo("dkz") ) return -0.1;
  if (parname.EqualTo("lZ") ) return -0.075;

  return -999;
}

//======================================================================

float parmax(const TString& parname)
{
  if (parname.EqualTo("cwww") )  return 20;
  if (parname.EqualTo("ccw") ) return 35;
  if (parname.EqualTo("cb") ) return 150;
  if (parname.EqualTo("dg1z") ) return 0.2;
  if (parname.EqualTo("dkz") ) return 0.15;
  if (parname.EqualTo("lZ") ) return 0.075;


  return -999;
}

//======================================================================

float parinc(const TString& parname)
{
  if (parname.EqualTo("cwww") )  return 0.4;
  if (parname.EqualTo("ccw") ) return 0.6;
  if (parname.EqualTo("cb") ) return 1.2;
  if (parname.EqualTo("dg1z") ) return 0.0002;
  if (parname.EqualTo("dkz") ) return 0.0002;
  if (parname.EqualTo("lZ") ) return 0.0002;

  return -999;
}

//======================================================================

int parbins(const TString& parname)
{
  int parbins = 1 + ((parmax(parname)-parmin(parname))/parinc(parname));
  return parbins;
}

//======================================================================

void getFileNames(const string& fileglob,
		  vector<TString>& fnames)
{
  // File globbing pattern can select multiple files
  //
  glob_t globbuf;
  int stat = glob (fileglob.c_str(), GLOB_MARK, NULL, &globbuf);
  if (stat) {
    switch (stat) {
    case GLOB_NOMATCH: cerr << "No file matching glob pattern "; break;
    case GLOB_NOSPACE: cerr << "glob ran out of memory "; break;
    case GLOB_ABORTED: cerr << "glob read error "; break;
    default: cerr << "unknown glob error stat=" << stat << " "; break;
    }
    cerr << fileglob << endl;
    exit(-1);
  }
  cout<<globbuf.gl_pathc<<" files match the glob pattern"<<endl;
  for (size_t i=0; i<globbuf.gl_pathc; i++) {
    char *path = globbuf.gl_pathv[i];
    if (!strncmp(&path[strlen(path)-6],".root",5)) {
      cerr << "non-root file found in glob, skipping: " << path << endl;
    } else {
      fnames.push_back(TString(path));
    }
  }
  globfree(&globbuf);
}                                                        // getFileNames

//======================================================================

double extractParValue(const char *par,
		       const TString& file)
{
  TRegexp re(Form("%s_\\-?[0-9]*\\.?[0-9]*e?\\-?[0-9]*",par));

  int len = 0;
  int i = re.Index(file,&len);

  i  +=strlen(par)+1;
  len-=strlen(par)+1;

  //cout << file << "\t" << i << "\t" << len << "\t" << TString(file(i,len)) << endl;

  if (i > 0 && len > 0)
    return (TString(file(i,len))).Atof();

  cerr << "Could not find a valid value for par " <<  par << " in filename: " << file << endl;

  return -9e99;
}                                                     // extractParValue

//======================================================================
// adapted from HiggsAnalysis/CombinedLimit/test/plotting/bandUtils.cxx
//
int  getBands(TFile *file, int doSyst, int whichChannel,
	      double& obs,	 double& median,
	      double& s68hi,	 double& s68lo,
	      double& s95hi,	 double& s95lo)
{
#if 0
  printf("getBands(file=%lx,doSyst=%d,whichChannel=%d)\t",
  	 (unsigned long)file,doSyst,whichChannel);
#endif
  if (file == 0) return 0;

  TTree *t = (TTree *) file->Get("limit");

  if (t == 0) { 
    std::cerr<<"TFile "<<file->GetName()<<" does not contain the tree"<<std::endl;
    return 0; 
  }

  Double_t limit, limitErr = 0;
  Int_t syst, iChannel, iToy;

  t->SetBranchAddress("limit", &limit);
  if (t->GetBranch("limitErr"))
    t->SetBranchAddress("limitErr", &limitErr);

  t->SetBranchAddress("syst", &syst);
  t->SetBranchAddress("iChannel", &iChannel);
  t->SetBranchAddress("iToy", &iToy);

  std::vector<double>  dataset;
  std::vector<double>  errors;

  for (size_t i = 0, n = t->GetEntries(); i < n; ++i) {
    t->GetEntry(i);
    //if (!i) printf("%6d mh=%.1f iChannel=%d syst=%d limit=%8.3f +/- %8.3f toy=%5d", 
    //i, iChannel,syst,limit, limitErr, iToy);
    if (syst     != doSyst)            { /* cout<<"syst"<<endl;    */ continue; }
    if (iChannel != whichChannel)      { /* cout<<"channel"<<endl; */ continue; }
    if (limit == 0)                    { /* cout<<"limit==0"<<endl;*/ continue; }
      
    dataset.push_back(limit);
    errors.push_back(limitErr);
  }

  int nd = dataset.size();
  if (!nd) {
    cerr << "Zero entries for " << file->GetName() << endl;
    return 0;
  }

  std::sort(dataset.begin(), dataset.end());
    
  median = (dataset.size() % 2 == 0 ? 0.5*(dataset[nd/2]+dataset[nd/2+1]) : dataset[nd/2]);
  double mean = 0;
  for (int j = 0; j < nd; ++j) mean += dataset[j]; 
  mean /= nd;
  s68lo = dataset[             floor(0.5*nd*(1-0.68)+0.5)        ];
  s68hi = dataset[std::min(int(floor(0.5*nd*(1+0.68)+0.5)), nd-1)];
  s95lo = dataset[             floor(0.5*nd*(1-0.95)+0.5)        ];
  s95hi = dataset[std::min(int(floor(0.5*nd*(1+0.95)+0.5)), nd-1)];

  if (nd == 1) {
    obs = mean;
    if (errors.size() == 1) {
      s68lo = mean - errors[0];
      s68hi = mean + errors[0];
    } else {
      // could happen if limitErr is not available
      s68lo = s68hi = mean;
    }
  }
#if 0
  printf("nd=%d, mean=%.3g, median=%.3g, s68lo=%.3g, s68hi=%.3g\n",
	 nd, mean,median,s68lo,s68hi);
#endif
  return nd;
}                                                            // getBands

//======================================================================

void fillGraphsFromFiles( const TString& par1,
			  const TString& par2,
			  const vector<TString>& fnames,
			  vector<string>&  keys,
			  map<string,TGraph2D *>& m_graphs)
{
  keys.push_back("-2s");
  keys.push_back("-1s");
  keys.push_back("median");
  keys.push_back("+1s");
  keys.push_back("+2s");
  keys.push_back("obs");

  for (int i=0; i<6; i++) {
    m_graphs[keys[i]] = new TGraph2D();
    m_graphs[keys[i]]->SetName(Form("graph2D%s",keys[i].c_str()));
  }

  int nobs=0, nexp=0;
  for( size_t i=0; i<fnames.size(); i++) {
    
    double par1val = extractParValue(par1,fnames[i]);
    double par2val = extractParValue(par2,fnames[i]);

    //cout << par1val << "\t" << par2val << endl;

    if (par1val == -9e99 ||
	par2val == -9e99)
      continue;
    
    TFile *f = new TFile(fnames[i]);

    double obs,median,s68hi,s68lo,s95hi,s95lo;
    int num = getBands(f,1,0,obs,median,s68hi,s68lo,s95hi,s95lo);

    switch (num) {
    case 0: break;
    case 1:
      //cout << "SetPoint(i="<<nobs<<",par1="<<par1val<<",par2="<<par2val<<",obs="<<obs<<");"<<endl;
      m_graphs["obs"]->SetPoint(nobs,par1val,par2val,obs);
      nobs++;
      break;
    default:
      m_graphs["+2s"]->SetPoint(nexp,par1val,par2val,s95hi);
      m_graphs["-2s"]->SetPoint(nexp,par1val,par2val,s95lo);
      m_graphs["+1s"]->SetPoint(nexp,par1val,par2val,s68hi);
      m_graphs["-1s"]->SetPoint(nexp,par1val,par2val,s68lo);
      m_graphs["median"]->SetPoint(nexp,par1val,par2val,median);
      nexp++;
      break;
    }
        
    f->Close();

    delete f;

    //if (!(i%10)) cout << i << " " << std::flush;


  } // file loop

#if 0
  TGraph2D *gobs = m_graphs["obs"];
  cout << "obs has " << gobs->GetN() << " points" << endl;
  cout << "npx = " << gobs->GetNpx() << endl;
  cout << "npy = " << gobs->GetNpy() << endl;
  cout << "xmin = " << gobs->GetXmin() << endl;
  cout << "xmax = " << gobs->GetXmax() << endl;
  cout << "ymin = " << gobs->GetYmin() << endl;
  cout << "ymax = " << gobs->GetYmax() << endl;
  cout << "zmin = " << gobs->GetZmin() << endl;
  cout << "zmax = " << gobs->GetZmax() << endl;

  double *xvec = gobs->GetX();
  double *yvec = gobs->GetY();
  double *zvec = gobs->GetZ();
  for (int i=0; i<gobs->GetN(); i++)
    printf("%lf\t%lf\t%lf\n",xvec[i],yvec[i],zvec[i]);
#endif
}                                                 // fillGraphsFromFiles

//======================================================================

void fillGraphsFromTextTables( const TString& fname,
			       vector<string>&  keys,
			       map<string,TGraph2D *>& m_graphs)
{
  keys.push_back("-2s");
  keys.push_back("-1s");
  keys.push_back("median");
  keys.push_back("+1s");
  keys.push_back("+2s");
  keys.push_back("obs");

  m_graphs["obs"] = new TGraph2D(fname,"%lg %lg %lg %*lg %*lg %*lg %*lg %*lg");
  m_graphs["obs"]->SetName("graph2Dobs");

  m_graphs["-2s"] = new TGraph2D(fname,"%lg %lg %*lg %*lg %*lg %*lg %*lg %lg");
  m_graphs["-2s"]->SetName("graph2D-2s");

  m_graphs["-1s"] = new TGraph2D(fname,"%lg %lg %*lg %*lg %*lg %*lg %lg %*lg");
  m_graphs["-1s"]->SetName("graph2D-1s");

  m_graphs["median"] = new TGraph2D(fname,"%lg %lg %*lg %*lg %*lg %lg %*lg %*lg");
  m_graphs["median"]->SetName("graph2Dmedian");

  m_graphs["+1s"] = new TGraph2D(fname,"%lg %lg %*lg %*lg %lg %*lg %*lg %*lg");
  m_graphs["+1s"]->SetName("graph2D+1s");

  m_graphs["+2s"] = new TGraph2D(fname,"%lg %lg %*lg %lg %*lg %*lg %*lg %*lg");
  m_graphs["+2s"]->SetName("graph2D+2s");

  TCanvas *canv = new TCanvas("dummy2","dummy2",500,500);
}                                            // fillGraphsFromTextTables

//======================================================================

void fillGraphsFromFilesAsymp( const TString& par1,
			       const TString& par2,
			       const vector<TString>& fnames,
			       vector<string>&  keys,
			       map<string,TGraph2D *>& m_graphs)
{
  keys.push_back("-2s");
  keys.push_back("-1s");
  keys.push_back("median");
  keys.push_back("+1s");
  keys.push_back("+2s");
  keys.push_back("obs");


  for (int i=0; i<6; i++) {
    m_graphs[keys[i]] = new TGraph2D();
    m_graphs[keys[i]]->SetName(Form("graph2D%s",keys[i].c_str()));
  }

  for( size_t i=0; i<fnames.size(); i++) {
    
    double par1val = extractParValue(par1,fnames[i]);
    double par2val = extractParValue(par2,fnames[i]);

    //cout << par1val << "\t" << par2val << endl;

    if (par1val == -9e99 ||
	par2val == -9e99)
      continue;
    
    TFile *f = new TFile(fnames[i]);

    TTree *lTree = (TTree *)f->Get("limit");
    if (!lTree) {
      f->Close();
      delete f;
      cerr << "didn't find tree limit in file " << fnames[i] << endl;
      continue;
    }
    
    for (int j=0; j<6; j++) {
      lTree->GetEntry(j);
      double rval = lTree->GetLeaf("limit")->GetValue();
      m_graphs[keys[j]]->SetPoint(i,par1val,par2val,rval);
    }
        
    f->Close();

    delete f;

    //if (!(i%10)) cout << i << " " << std::flush;


  } // file loop
}                                            // fillGraphsFromFilesAsymp

//======================================================================

void fillGraphsFromFilesDeltaNLL( const TString& par1name,
				  const TString& par2name,
				  const vector<TString>& fnames,
				  vector<string>&  keys,
				  map<string,TGraph2D *>& m_graphs)
{

  std::cout << "fillGraphsFromFilesDeltaNLL 1" << std::endl;

  keys.push_back("exp68");
  keys.push_back("exp95");
  keys.push_back("exp99");

  // uncommented below to plot observed!
  keys.push_back("obs95");

  TGraph2D *grobs = new TGraph2D();
  TGraph2D *grexp = new TGraph2D();

  m_graphs["obs95"] = grobs;
  m_graphs["exp95"] = grexp;

  grobs->SetName("graph2Dobs95");
  grexp->SetName("graph2Dexp95");

  Int_t nobs=0, nexp=0;

  for( size_t i=0; i<fnames.size(); i++) {
    
    TFile *f = new TFile(fnames[i]);
    TTree *t = (TTree *) f->Get("limit");

    if (!t) { 
      std::cerr<<"TFile "<<f->GetName()<<" does not contain the tree"<<std::endl;
      return;
    }
    cout << fnames[i] << " has limit tree with " << t->GetEntries() << " entries." << endl;

    Float_t deltaNLL, par1, par2;
    Int_t iToy;

    t->SetBranchAddress("iToy", &iToy);
    t->SetBranchAddress("deltaNLL", &deltaNLL);
    t->SetBranchAddress(par1name, &par1);
    t->SetBranchAddress(par2name, &par2);

    for (size_t j = 0, n = t->GetEntries(); j < n; ++j) {
      t->GetEntry(j);
      printf ("%d\r",j);
      if( !iToy){
	//	cout <<"!iToy" << endl;
	grobs->SetPoint(nobs++,par1,par2,2*deltaNLL);
	//	cout <<"grobs->SetPoint("<<nobs++<<","<<par1<<","<< par2<< ","<< 2*deltaNLL << endl;
      }
      else if (iToy == -1) {
	//	cout <<"iToy == -1" << endl;
	grexp->SetPoint(nexp++,par1,par2,2*deltaNLL);
      }
      else {
	cerr << "Unexpected value for iToy, = " << iToy << endl;
	exit(-1);
      }
    } // tree entry loop

    f->Close();
    delete f;

  } // file loop
  cout << endl;

  m_graphs["exp68"] = (TGraph2D*)grexp->Clone("graph2Dexp68");
  m_graphs["exp99"] = (TGraph2D*)grexp->Clone("graph2Dexp99");

#if 0
  TCanvas *canv = new TCanvas("tester","tester",500,500);
  cout << grexp->GetN()<<" points. " <<endl;
  grexp->Draw("TRI"); // cont 5z list");
#endif
}                                         // fillGraphsFromFilesDeltaNLL

//======================================================================

void collectContours(map<string,TGraph2D *>& m_graphs,
		     const vector<string>&  keys,
		     map<string,double>& m_contourlevels,
		     map<string,TList *>& m_contours)
{
  cout << "CollectContours" << endl;

  TCanvas *canv = new TCanvas("dummy","dummy",100,100);
  //canv->Divide(3,2);

  std::cout << "keys.size() = " << keys.size() << std::endl;

  //process TGraph2Ds into contours at levels m_contourlevels
  for (size_t i=0; i<keys.size(); i++) {
    double clev = m_contourlevels[keys[i]];
    TGraph2D *gr2d = m_graphs[keys[i]];
    std::cout << "gr2d = " << gr2d << std::endl;
    std::cout << "gr2d->GetN() = " << gr2d->GetN() << std::endl;


    if (gr2d && (gr2d->GetN() > 0)) {
      gr2d->GetHistogram()->SetContour(1, &clev);
      //canv->cd(i+1);
      cout << "drawing... " << endl;

      gr2d->Draw("CONT LIST"); // it's stupid, but only "CONT" will generate the list
      gPad->Update();

      TObjArray *contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
      assert(contours);

      TList *newlist = 0;
      for (int ci=0; ci<contours->GetEntriesFast(); ci++) {
	TList *contLevel = (TList*)contours->At(ci);
	printf("%s: Contour %d has %d Graphs\n", keys[i].c_str(), ci, contLevel->GetSize());

	if (contLevel->GetSize()) {
	  assert(contLevel->First());
	  if (!newlist) newlist = new TList();
	  TGraph *curv = (TGraph*)(contLevel->First());

	  for (int j=0; j<contLevel->GetSize(); j++) {
	    newlist->Add((TGraph *)(curv->Clone()));
	    curv=(TGraph *)(contLevel->After(curv));
	  }
	}
      } // contour loop

      cout << "Inserting contour list for "<< keys[i] << " newlist="<<newlist<<endl;
      m_contours[keys[i]] = newlist;

    } // if (gr2d)
  } // key loop

  //delete canv;
}                                                     // collectContours

//======================================================================
// "Brazilian Flag" style

void
draw2DLimitBFstyle(map<string,TList *>& m_contours,
		     const TString& par1,
		     const TString& par2,
		     const TString& plotprefix,
		     TLegend *legend)
{

  //from here we build the two-dimensional aTGC limit

  TCanvas *finalPlot = new TCanvas("final","limits",500,500);
  finalPlot->cd();

  cout << "Drawing +2s" << endl;

  TList *contLevel = m_contours["+2s"];
  TGraph *curv;

  assert(contLevel);

  curv = (TGraph*)(contLevel->First());

  //curv->GetYaxis()->SetRangeUser(-1.25*curv->GetYaxis()->GetXmax(),
	  			   //+2.0*curv->GetYaxis()->GetXmax());
  //curv->GetYaxis()->SetRangeUser(-0.1,0.15);
  curv->GetYaxis()->SetRangeUser(parmin(par2),parmax(par2));

  curv->SetTitle();
  curv->GetXaxis()->SetTitle(par2latex(par1));
  curv->GetXaxis()->SetTitleFont(42);
  curv->GetYaxis()->SetTitle(par2latex(par2));
  curv->GetYaxis()->SetTitleFont(42);
  curv->GetYaxis()->SetTitleOffset(1.20);

  for (int i=0; i<contLevel->GetSize(); i++) {
    assert(curv);
    curv->SetLineColor(kYellow);
    curv->SetFillColor(kYellow);
    curv->GetXaxis()->SetLimits(parmin(par1),parmax(par1));
    if (!i) {
      curv->Draw("ACF");
      legend->AddEntry(curv,"#pm 2#sigma","F");
    } else 
      curv->Draw("SAME CF");
    curv=(TGraph *)(contLevel->After(curv));
  }

  cout << "Drawing +1s" << endl;
  
  contLevel = m_contours["+1s"];

  curv = (TGraph*)(contLevel->First());

  for (int i=0; i<contLevel->GetSize(); i++) {
    curv->SetLineColor(kGreen);
    curv->SetFillColor(kGreen);
    curv->Draw("SAME CF");
    if (!i) legend->AddEntry(curv,"#pm 1#sigma","F");
    curv=(TGraph *)(contLevel->After(curv));
  }

  cout << "Drawing -1s" << endl;

  contLevel = m_contours["-1s"];
  curv = (TGraph*)(contLevel->First());
  for (int i=0; i<contLevel->GetSize(); i++) {
    curv->SetLineColor(kYellow);
    curv->SetFillColor(kYellow);
    curv->Draw("SAME CF");
    curv=(TGraph *)(contLevel->After(curv));
  }

  cout << "Drawing -2s" << endl;
  
  contLevel = m_contours["-2s"];

  if (!contLevel)
    //  this can happen more often for this contour if there is insufficient
    // sensitivity close to the SM
    cerr << "No contour level for +2s, have to fill in the central region" << endl;
  else {
    curv = (TGraph*)(contLevel->First());
    for (int i=0; i<contLevel->GetSize(); i++) {
      curv->SetFillColor(kWhite);
      curv->SetLineColor(kYellow);
      curv->Draw("SAME CF");
      curv=(TGraph *)(contLevel->After(curv));
    }
  }
  cout << "Drawing median" << endl;
  
  curv = (TGraph*)(m_contours["median"]->First());
  curv->SetLineColor(kBlack);
  curv->SetLineWidth(2);
  curv->SetLineStyle(2);
  curv->Draw("SAME C");

  legend->AddEntry(curv,"Expected","L");
  
  cout << "Drawing obs" << endl;

  contLevel = m_contours["obs"];
  curv = (TGraph*)(contLevel->First());
  for (int i=0; i<contLevel->GetSize(); i++) {
    curv->SetLineColor(kBlack);
    curv->SetLineWidth(2);
    curv->Draw("SAME C");
    if (!i) legend->AddEntry(curv,"Observed","L");
    curv=(TGraph *)(contLevel->After(curv));
  }

  
  TGraph *SMpoint = new TGraph(1);
  SMpoint->SetPoint(1,0,0);
  SMpoint->Draw("SAME Po");
  
  // smLabel = TPaveText(0,
  //                     m_contours["-2s"]->GetYaxis()->GetXmax()/8,
  //                     m_contours["-2s"]->GetXaxis()->GetXmax()/3->5,
  //                     -m_contours["-2s"]->GetYaxis()->GetXmax()/8);
  // smLabel->SetFillStyle(0);
  // smLabel->SetBorderSize(0);
  // smLabel->AddText(" SM");
  // smLabel->Draw();

  legend->Draw();

  TPaveText *text = new TPaveText(0.566,0.87,0.965,1.101,"NDC");
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->AddText(Form("95%% CL Limit on %s and %s",par2latex(par1).Data(),par2latex(par2).Data()));
  text->AddText(0,0.35,Form("#intL dt= %.1f fb^{-1}, #sqrt{s} = %d TeV",intlumifbinv,beamcometev));
  text->Draw();

  // text2 = TPaveText(0.155,0.199,0.974,0.244,"NDC");
  // text2->SetFillStyle(0);
  // text2->SetBorderSize(0);
  // text2->AddText("Values outside contour excluded");
  // text2->Draw();

  //text3 = TPaveText(0.506,0.699,0.905,0.758,"NDC");
  //text3->SetFillStyle(0);
  //text3->SetBorderSize(0);
  //text3->AddText(options.flavorText);
  //text3->Draw();    
  
  finalPlot->RedrawAxis();
  finalPlot->ResetAttPad();
  finalPlot->Update();

  finalPlot->Draw();
  finalPlot->Update();
  finalPlot->Modified();
  finalPlot->Update();
  finalPlot->Print(Form("%s.pdf",plotprefix.Data()));
  finalPlot->Print(Form("%s.eps",plotprefix.Data()));
  finalPlot->Print(Form("%s.png",plotprefix.Data()));

}                                                  // draw2DlimitBFstyle

//======================================================================

void
draw2DLimitContours(map<string,TList *>& m_contours,
		    const TString& par1,
		    const TString& par2,
		    const TString& plotprefix,
		    TLegend *legend,
		    float par1_bestfit,
		    float par2_bestfit)
{

  //from here we build the two-dimensional aTGC limit

  TCanvas *finalPlot = new TCanvas("final","limits",500,500);
  finalPlot->cd();

  cout << "Drawing expected 68%" << endl;

  TList *contLevel = m_contours["exp68"];
  TGraph *curv;

  std::cout << "m_contours.size() = " << m_contours.size() << std::endl;

  for (map<string,TList *>::const_iterator iter = m_contours.begin(); iter != m_contours.end(); iter++ ){
    std::cout << "iter->first = " << iter->first << std::endl;
    std::cout << "iter->second = " << iter->second << std::endl;
  }

  std::cout << "contLevel = " << contLevel << std::endl;

  assert(contLevel);

  curv = (TGraph*)(contLevel->First());

  curv->GetXaxis()->SetLimits(parmin(par1),parmax(par1));
  curv->GetYaxis()->SetRangeUser(parmin(par2),parmax(par2));

  curv->SetTitle();
  curv->GetXaxis()->SetTitle(par2latex(par1));
  curv->GetXaxis()->SetTitleFont(42);
  curv->GetYaxis()->SetTitle(par2latex(par2));
  curv->GetYaxis()->SetTitleFont(42);
  curv->GetYaxis()->SetTitleOffset(1.20);

  legend->SetNColumns(2);

  for (int i=0; i<contLevel->GetSize(); i++) {
    assert(curv);
    curv->SetLineColor(kBlue);
    curv->SetLineWidth(2);
    curv->SetLineStyle(9);
    if (!i) {
      curv->Draw("AC");
      legend->AddEntry(curv,"Expected 68% C.L.","L");
    } else 
      curv->Draw("SAME C");
    curv=(TGraph *)(contLevel->After(curv));
  }

  cout << "Drawing expected 95%" << endl;
  
  contLevel = m_contours["exp95"];

  curv = (TGraph*)(contLevel->First());

  for (int i=0; i<contLevel->GetSize(); i++) {
    curv->SetLineColor(kGreen);
    curv->SetLineWidth(2);
    curv->SetLineStyle(9);
    curv->Draw("SAME C");
    if (!i) legend->AddEntry(curv,"Expected 95% C.L.","L");
    curv=(TGraph *)(contLevel->After(curv));
  }

  cout << "Drawing expected 99%" << endl;

  contLevel = m_contours["exp99"];
  curv = (TGraph*)(contLevel->First());
  for (int i=0; i<contLevel->GetSize(); i++) {
    curv->SetLineColor(kRed);
    curv->SetLineWidth(2);
    curv->SetLineStyle(9);
    curv->Draw("SAME C");
    if (!i) legend->AddEntry(curv,"Expected 99% C.L.","L");
    curv=(TGraph *)(contLevel->After(curv));
  }

  
  contLevel = m_contours["obs95"];

  if (contLevel) {
    cout << "Drawing obs95" << endl;

    curv = (TGraph*)(contLevel->First());

    for (int i=0; i<contLevel->GetSize(); i++) {
      curv->Draw("SAME C");
      curv->SetLineWidth(3);
      if (!i) legend->AddEntry(curv,"Observed 95% C.L.","L");
      curv=(TGraph *)(contLevel->After(curv));
    }
  }

  
  TGraph *SMpoint = new TGraph(1);
  SMpoint->SetPoint(1,0,0);
  SMpoint->Draw("SAME P");
  legend->AddEntry(SMpoint,"SM","P");
  TGraph *bestfit = new TGraph(1);
  bestfit->SetPoint(1,par1_bestfit,par2_bestfit);
  bestfit->Draw("SAME *");
  legend->AddEntry(bestfit,"Best fit value","P");
  
  //smLabel = TPaveText(0,
  //                    m_contours["-2s"]->GetYaxis()->GetXmax()/8,
  //                    m_contours["-2s"]->GetXaxis()->GetXmax()/3->5,
  //                    -m_contours["-2s"]->GetYaxis()->GetXmax()/8);
  //smLabel->SetFillStyle(0);
  //smLabel->SetBorderSize(0);
  //smLabel->AddText(" SM");
  //smLabel->Draw();

  legend->Draw();

  TPaveText *text = new TPaveText(0.566,0.89,0.965,1.13,"NDC");
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextFont(42);
  text->SetTextSize(0.05);
  //text->AddText(Form("95%% CL Limit on %s and %s",par2latex(par1).Data(),par2latex(par2).Data()));
  text->AddText(0,0.35,Form("2.3 fb^{-1} (%d TeV)",beamcometev));
  text->Draw();

  // text2 = TPaveText(0.155,0.199,0.974,0.244,"NDC");
  // text2->SetFillStyle(0);
  // text2->SetBorderSize(0);
  // text2->AddText("Values outside contour excluded");
  // text2->Draw();

  //text3 = TPaveText(0.506,0.699,0.905,0.758,"NDC");
  //text3->SetFillStyle(0);
  //text3->SetBorderSize(0);
  //text3->AddText(options.flavorText);
  //text3->Draw();    
  
  gPad->SetGrid(1,1);

  finalPlot->RedrawAxis();
  finalPlot->ResetAttPad();
  finalPlot->Update();

  finalPlot->Draw();
  finalPlot->Update();
  finalPlot->Modified();
  finalPlot->Update();
  finalPlot->Print(Form("%s.pdf",plotprefix.Data()));
  finalPlot->Print(Form("%s.eps",plotprefix.Data()));
  //finalPlot->Print(Form("%s.png",plotprefix.Data()));

}                                                 // draw2DlimitContours

//======================================================================

void
draw1DLimit(map<string,TGraph2D *> m_graphs,
	    const TString& parname,
	    const TString& plotprefix,
	    int      npts,
	    double   boundScale,   // used to exclude region closest to SM from plotting
	    double   exclusionlimit,
	    bool     isX,
	    TLegend *legend)
{
  TCanvas *c1 = new TCanvas(Form("%slimit",parname.Data()),
			    Form("%slimit",parname.Data()),
			    500,500);

  map<string,TGraphAsymmErrors *> m_limits1d;

  m_limits1d["2s"]     = new TGraphAsymmErrors();
  m_limits1d["1s"]     = new TGraphAsymmErrors();
  m_limits1d["-median"] = new TGraphAsymmErrors();
  m_limits1d["-obs"]    = new TGraphAsymmErrors();
  m_limits1d["+median"] = new TGraphAsymmErrors();
  m_limits1d["+obs"]    = new TGraphAsymmErrors();

  double parSize = parmax(parname) - parmin(parname);

  double lowerLimit=0.0,        upperLimit=0.0;
  bool   lowerLimitFound=false, upperLimitFound=false;

  double parcutoff = parmin(parname)*boundScale;

  double bound = 0.0;    // bound =  the max y value for the minus side
  if (isX)
    bound = m_graphs["-2s"]->Interpolate(parcutoff,0.0);
  else
    bound = m_graphs["-2s"]->Interpolate(0.0,parcutoff);

  printf ("par=%s, min=%f, max=%f, boundScale=%lf, bound=%lf, parcutoff=%lf\n",
	  parname.Data(),parmin(parname), parmax(parname), boundScale, bound, parcutoff);

  int nnegpts=0,npospts=0;
  for (int i=0; i<npts; i++) {
    double obs, median, p1s, m1s, p2s, m2s;
    double parval = parmin(parname) + i*parSize/npts;
        
    if (isX) {
      obs  = m_graphs["obs"]->Interpolate(parval,0.0);
      median = m_graphs["median"]->Interpolate(parval,0.0);
      p1s  = m_graphs["+1s"]->Interpolate(parval,0.0);
      m1s  = m_graphs["-1s"]->Interpolate(parval,0.0);
      p2s  = m_graphs["+2s"]->Interpolate(parval,0.0);
      m2s  = m_graphs["-2s"]->Interpolate(parval,0.0);
    } else {
      obs  = m_graphs["obs"]->Interpolate(0.0,parval);
      median = m_graphs["median"]->Interpolate(0.0,parval);
      p1s  = m_graphs["+1s"]->Interpolate(0.0,parval);
      m1s  = m_graphs["-1s"]->Interpolate(0.0,parval);
      p2s  = m_graphs["+2s"]->Interpolate(0.0,parval);
      m2s  = m_graphs["-2s"]->Interpolate(0.0,parval);
    }

    //print m2s, m1s, median, p1s, p2s
        
    if (obs > exclusionlimit && !lowerLimitFound && (parval < 0)) {
      lowerLimit = parval;
      lowerLimitFound = true;
    } else if (obs < exclusionlimit && !upperLimitFound && (parval > 0)) {
      upperLimit = parval;
      upperLimitFound = true;
    }

    if ( // (m2s < bound) &&
	 (fabs(parval) > fabs(parcutoff))
       ) {
      //par1 observed limit
      //obs and median
      if (parval < 0) {
	m_limits1d["-median"]->SetPoint(nnegpts,parval,median);
	m_limits1d["-obs"]->SetPoint(nnegpts++,parval,obs);
      } else {
	m_limits1d["+median"]->SetPoint(npospts,parval,median);
	m_limits1d["+obs"]->SetPoint(npospts++,parval,obs);
      }
      // one sigma expected
      m_limits1d["1s"]->SetPoint(i,parval,median);
      m_limits1d["1s"]->SetPointError(i,0,0,median-m1s,p1s-median);
      //two sigma expected    
      m_limits1d["2s"]->SetPoint(i,parval,median);
      m_limits1d["2s"]->SetPointError(i,0,0,median-m2s,p2s-median);
    } else {
      m_limits1d["1s"]->SetPoint(i,parval,bound+0.1);
      m_limits1d["1s"]->SetPointError(i,0,0,0,0);
      m_limits1d["2s"]->SetPoint(i,parval,bound+0.1);
      m_limits1d["2s"]->SetPointError(i,0,0,0,0);
    }
  }  // npts loop

  //print "95% CL on"+" %s = [%.3g,%.3g]"%(par,lowerLimit,upperLimit)

  c1->cd();
  c1->SetLogy();

  m_limits1d["2s"]->SetFillColor(kYellow);
  m_limits1d["2s"]->Draw("A E3");
    
  m_limits1d["1s"]->SetFillColor(kGreen);
  m_limits1d["1s"]->Draw("SAME E3");
    
  m_limits1d["-median"]->SetLineStyle(2);
  m_limits1d["-median"]->SetLineWidth(2);
  m_limits1d["-median"]->Draw("SAME C");
    
  m_limits1d["+median"]->SetLineStyle(2);
  m_limits1d["+median"]->SetLineWidth(2);
  m_limits1d["+median"]->Draw("SAME C");
    
  m_limits1d["-obs"]->SetLineWidth(2);
  m_limits1d["-obs"]->Draw("SAME C");
    
  m_limits1d["+obs"]->SetLineWidth(2);
  m_limits1d["+obs"]->Draw("SAME C");
    
  //titles
  if (exclusionlimit==1)
    m_limits1d["2s"]->GetYaxis()->SetTitle("95% CL limit on #sigma/#sigma_{aTGC}");
  else
    m_limits1d["2s"]->GetYaxis()->SetTitle(Form("p_{S+B}(%s)",par2latex(parname).Data()));

  m_limits1d["2s"]->GetYaxis()->SetTitleFont(42);
  m_limits1d["2s"]->GetXaxis()->SetTitle(par2latex(parname));
  m_limits1d["2s"]->GetXaxis()->SetTitleFont(42);
  
  m_limits1d["2s"]->GetYaxis()->SetRangeUser(m_limits1d["2s"]->GetYaxis()->GetXmin()*0.75,10); //,bound);
  m_limits1d["2s"]->GetXaxis()->SetRangeUser(parmin(parname)*0.985,parmax(parname)*0.96);

#if 0    
  legend->SetX1NDC(0.43);
  legend->SetY1NDC(0.43);
  legend->SetX2NDC(0.75);
  legend->SetY2NDC(0.65);
#endif
  legend->Draw();
  
  TPaveText *text1d = new TPaveText(0.566,0.87,0.965,1.101,"NDC");
  //TPaveText *text1d = new TPaveText(0.359,0.24,0.758,0.44,"NDC");
  text1d->SetFillStyle(0);
  text1d->SetBorderSize(0);
  text1d->AddText(Form("95%% CL Limit on #bf{%s}",par2latex(parname).Data()));
  text1d->AddText(0,0.35,Form("%f fb^{-1} (%d TeV)", intlumifbinv,beamcometev));
  text1d->Draw();
    
  // text3.SetX1NDC(0.357);
  // text3.SetY1NDC(0.246);
  // text3.SetX2NDC(0.756);
  // text3.SetY2NDC(0.305);
  // text3.Draw();
  
  TPaveText *obslimtext = new TPaveText(0.357,0.197,0.754,0.246,"NDC");
  obslimtext->SetFillStyle(0);
  obslimtext->SetBorderSize(0);
  obslimtext->AddText(Form("%.3f < %s  < %.3f",
			   lowerLimit,par2latex(parname).Data(),upperLimit));
  obslimtext->Draw();
  
  //lowLimitLine = TLine(lowerLimit,m_limits1d["2s"]->GetYaxis()->GetXmin()*0.75,
  //                     lowerLimit,1)
  //lowLimitLine->SetLineColor(14)
  //lowLimitLine->SetLineWidth(2)
  //lowLimitLine->Draw()
  //upLimitLine = TLine(upperLimit,m_limits1d["2s"]->GetYaxis()->GetXmin()*0.75,
  //                    upperLimit,1)
  //upLimitLine->SetLineColor(14)
  //upLimitLine->SetLineWidth(2)
  //upLimitLine->Draw()

  TLine *oneLine = new TLine(parmin(parname)*0.985,exclusionlimit,
			     parmax(parname)*0.960,exclusionlimit);
  oneLine->SetLineStyle(9);
  oneLine->SetLineColor(14);
  oneLine->Draw();
    
  c1->Draw();
  c1->Update();
  c1->Modified();
  c1->Update();

  c1->Print(Form("%s.pdf",plotprefix.Data()));
  c1->Print(Form("%s.eps",plotprefix.Data()));
  c1->Print(Form("%s.png",plotprefix.Data()));
}                                                         // draw1Dlimit

//======================================================================
// to plot only expected limit, just omit the observed limit file from
// the fileglob

void atgcplotLimit_lZ_dkz()
{
  atgcstyle();

  gStyle->SetPalette(1);

  vector<TString> fnames;
  //getFileNames(fileglob, fnames);

  TString   par1;
  TString   par2;
  par1 = TString("dkz");
  par2 = TString("dg1z");

  fnames.push_back("higgsCombine2Par_"+par1+"0_"+par2+"0_exp.MultiDimFit.mH120.root");
  fnames.push_back("higgsCombine2Par_"+par1+"0_"+par2+"0_obs.MultiDimFit.mH120.root");

  assert(fnames.size());



  
  float par1_bestfit	= 0;
  float par2_bestfit	= 0;
  TTree * tree	= (TTree*)TFile::Open(fnames[1])->Get("limit");
  tree->SetBranchAddress(par1, &par1_bestfit);
  tree->SetBranchAddress(par2, &par2_bestfit);
  tree->GetEntry(0); 



  if (!par1.Length() || !par2.Length() ) {
    cerr << "Unknown coupling parameters in name " << fnames[0] << endl;
    exit(-1);
  }

  TString method("Other");
  if (fnames[0].Contains("Asymptotic"))
    method = TString("asympCLs");
  else if (fnames[0].Contains("MultiDimFit"))
    method = TString("deltaNLL");

  // try this:
  //  method = TString("asympCLs");
  

  cout << "Plotting " << par2 << " versus " << par1 << ", method = " << method << endl;

  vector<string> keys;
  map<string,double> m_contourlevels;
  map<string,TGraph2D *> m_graphs;

  if (method.EqualTo("asympCLs")) {
    fillGraphsFromFilesAsymp(par1,par2,fnames,keys,m_graphs);
    for (size_t i=0; i<keys.size(); i++)
      m_contourlevels[keys[i]] = 1;
  }
  if (method.EqualTo("deltaNLL")) {
    fillGraphsFromFilesDeltaNLL(par1,par2,fnames,keys,m_graphs);
    m_contourlevels["exp68"] = 2.3;
    m_contourlevels["exp95"] = 5.99;
    m_contourlevels["exp99"] = 9.21;
    m_contourlevels["obs95"] = 5.99;
  }
  else if (fnames.size() == 1) {  
    fillGraphsFromTextTables          (fnames[0],keys,m_graphs);
    for (size_t i=0; i<keys.size(); i++)
      m_contourlevels[keys[i]] = 0.05;
#if 0
    m_graphs["-2s"]->Draw("COLZ TEXT");
#else
    TGraph2D *gr = m_graphs["-2s"];
    TH2D *h2 = new TH2D("h2d","h2d",
			parbins(par1),parmin(par1)-parinc(par1)/2,parmax(par1)+parinc(par1)/2,
			parbins(par2),parmin(par2)-parinc(par2)/2,parmax(par2)+parinc(par2)/2);
    h2->FillN(gr->GetN(),gr->GetX(),gr->GetY(),gr->GetZ());
    h2->GetXaxis()->SetTitle(par2latex(par1));
    h2->GetYaxis()->SetTitle(par2latex(par2));
    h2->GetYaxis()->SetTitleOffset(1.1);
    h2->Draw("COLZ TEXT");
#endif
  }  else {
    fillGraphsFromFiles        (par1,par2,fnames,keys,m_graphs);
    for (size_t i=0; i<keys.size(); i++)
      m_contourlevels[keys[i]] = 1;
  }

  // try this
  /*
    fillGraphsFromFiles        (par1,par2,fnames,keys,m_graphs);
    for (size_t i=0; i<keys.size(); i++)
      m_contourlevels[keys[i]] = 1;
  */
  // return;

  // for limit in limits:
  //   limits[limit]->Print()


#if 0
  // for a first look
  TCanvas *canv2 = new TCanvas("two","two",800,600);
  canv2->Divide(3,2);
  canv2->cd(1);  m_graphs["+2s"]->Draw("TRI"); gPad->SetLogz(1);
  canv2->cd(2);  m_graphs["+1s"]->Draw("TRI"); gPad->SetLogz(1);
  canv2->cd(3);  m_graphs["median"]->Draw("TRI"); gPad->SetLogz(1);
  canv2->cd(4);  m_graphs["-1s"]->Draw("TRI"); gPad->SetLogz(1);
  canv2->cd(5);  m_graphs["-2s"]->Draw("TRI"); gPad->SetLogz(1);
  canv2->cd(6);  m_graphs["obs"]->Draw("TRI"); gPad->SetLogz(1);
#else
  //m_graphs["obs"]->Draw("TRI");

  map<string,TList *> m_contours;

  collectContours(m_graphs,keys,m_contourlevels,m_contours);

  //  TLegend *legend = new TLegend(0.212,0.686,0.554,0.917,"","NDC");
  TLegend *legend = new TLegend(0.2,0.7,0.89,0.89,"","NDC");
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetHeader("CMS Preliminary");
  //legend->SetHeader("CMS");
  legend->SetTextFont(42);

  TString plotprefix=Form("%s_%s_2dlimit_%s",par1.Data(),par2.Data(),method.Data());

  if (method.EqualTo("deltaNLL"))
    draw2DLimitContours(m_contours,par1,par2,plotprefix,legend,par1_bestfit,par2_bestfit);
  else
    draw2DLimitBFstyle(m_contours,par1,par2,plotprefix,legend);

#if 0
  plotprefix=Form("%s_1dlimit_%s",par1.Data(),method.Data());
  draw1DLimit(m_graphs,par1,plotprefix,1000,0.15,exclusion_limit,true,legend);

  plotprefix=Form("%s_1dlimit_%s",par2.Data(),method.Data());
  draw1DLimit(m_graphs,par2,plotprefix,1000,0.15,exclusion_limit,false,legend);
#endif
#endif
}
