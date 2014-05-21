#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <algorithm>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/HistogramAnalyzer.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"

#include "CMGTools/HtoZZ2l2nu/interface/TMVAUtils.h"

using namespace std;
using namespace top;

struct bTagSorter{
  bool operator() (PhysicsObject_Jet a, PhysicsObject_Jet b)   {   return (a.btag1>b.btag1);  }
} sortByBtag;

struct ptSorter{
  bool operator() (PhysicsObject_Lepton a, PhysicsObject_Lepton b)   {   return (a.pt()>b.pt());  }
} sortByPt;

//
int main(int argc, char* argv[])
{
  SelectionMonitor controlHistos; //plot storage
  HistogramAnalyzer histoAnalyzer; 

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
    return 0;
  }
  
  //configure                                                                                                                                                                                                                 
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  TString evurl               = runProcess.getParameter<std::string>("input");
  TString outUrl              = runProcess.getParameter<std::string>("outdir");
  TString kindir              = runProcess.getParameter<std::string>("kindir");
  bool isMC                   = runProcess.getParameter<bool>("isMC");
  int mcTruthMode             = runProcess.getParameter<int>("mctruthmode");
  TString dirname             = runProcess.getParameter<std::string>("dirName");
  bool saveSummaryTree        = runProcess.getParameter<bool>("saveSummaryTree");
  bool useMVA                 = runProcess.getParameter<bool>("useMVA");
  if(mcTruthMode!=1)          useMVA=false;
  edm::ParameterSet tmvaInput = runProcess.getParameter<edm::ParameterSet>("tmvaInput");
  TString studyTag            = tmvaInput.getParameter<std::string>("studyTag");
  TString weightsDir          = tmvaInput.getParameter<std::string>("weightsDir");
  std::vector<std::string> methodList = tmvaInput.getParameter<std::vector<std::string> >("methodList");
  std::vector<std::string> varsList    = tmvaInput.getParameter<std::vector<std::string> >("varsList");
  bool trainMVA                        = tmvaInput.getParameter<bool>("doTrain");

  //book histos
  controlHistos.addHistogram( new TH1F ("njets", ";Jets;Events", 6, 0.,6.) );
  controlHistos.addHistogram( new TH1F ("btags", ";b-tag multiplicity;Events", 6, 0.,6.) );
  for(int ibin=1; ibin<=controlHistos.getHisto("njets","all")->GetXaxis()->GetNbins(); ibin++)
    {
      TString label(""); label += ibin-1; 
      if(controlHistos.getHisto("njets","all")->GetXaxis()->GetNbins()==ibin) label ="#geq" + label;
      controlHistos.getHisto("njets","all")->GetXaxis()->SetBinLabel(ibin,label + " jets");
      controlHistos.getHisto("btags","all")->GetXaxis()->SetBinLabel(ibin,label + " b-tags");
    }
  
  controlHistos.addHistogram( new TH1F ("leadjet", "; Leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("met", "; #slash{E}_{T} [GeV/c]; Events / (10 GeV/c)", 40, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("ht", "; #sum_{jets} [GeV/c]; Events / (20 GeV/c)",40, 0.,800.) );
  controlHistos.addHistogram( new TH1F ("st", "; #sum_{leptons,E_{T}^{miss}} p_{T} [GeV/c]; Events / (20 GeV/c)",40, 0.,800.) );
  controlHistos.addHistogram( new TH1F ("sumpt", "; #sum_{leptons} p_{T} [GeV/c]; Events / (20 GeV/c)",25, 0.,500.) );
  controlHistos.addHistogram( new TH1F ("htlep", "; #sum_{jets,leptons,E_{T}^{miss}} [GeV/c]; Events / (20 GeV/c)",70, 0.,1400.) );
  controlHistos.addHistogram( new TH1F ("mtop", "; m_{Top} [GeV/c^{2}]; Events / (15 GeV/c^{2})", 40, 0.,450.) );
  controlHistos.addHistogram( new TH1F ("ptttbar", "; p_{t#bar{t}} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("afb", "; #Delta #eta(t,#bar{t}) ); Events / (0.1)", 100, -5.,5.) );
  controlHistos.addHistogram( new TH1F ("mttbar", "; Mass(t,#bar{t}) [GeV/c^{2}]; Events / (50 GeV/c^{2})", 40, 0.,2000.) );
  controlHistos.addHistogram( new TH2F ("mtopvsdilmass", "; m_{Top} [GeV/c^{2}]; Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmlj", "; m_{Top} [GeV/c^{2}]; Mass(l,j) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmet", "; m_{Top} [GeV/c^{2}]; #slash{E}_{T} [GeV/c]; Events", 100, 0.,500.,50, 0.,250.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmttbar", "; m_{Top} [GeV/c^{2}]; Mass(t,#bar{t}) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,2000.) );
  controlHistos.addHistogram( new TH2F ("mtopvsafb", "; m_{Top} [GeV/c^{2}]; #Delta #eta(t,#bar{t}) ); Events", 100, 0.,500.,100, -5.,5.) );
  controlHistos.addHistogram( new TH2F ("mttbarvsafb", "; Mass(t,#bar{t}) [GeV/c^{2}];#Delta #eta(t,#bar{t}); Events", 100, 0.,2000.,100,-5.,5.) );
  controlHistos.addHistogram( new TH1F("assignmentdecision",";Good decisions",2,0.,2.) );
  
  //MVA analysis
  controlHistos.addHistogram( new TH1F("kIntegral",";Solutions; Lepton-jet assignments",100,0,5000) );
  controlHistos.addHistogram( new TH1F("kMPV",";(mpv-median)/median; Lepton-jet assignments",100,-1.5,1.5) );
  controlHistos.addHistogram( new TH1F("kMean",";(mean-median)/median; Lepton-jet assignments",100,-1.5,1.5) );
  controlHistos.addHistogram( new TH1F("kRMS",";rms/median; Lepton-jet assignments",100,0,1.) );
  controlHistos.addHistogram( new TH1F("kSkewness",";skewness/median; Lepton-jet assignments",100,-1.,1.) );
  controlHistos.addHistogram( new TH1F("kKurtosis",";kurtosis/median; Lepton-jet assignments",100,-1.,1.) );
  controlHistos.addHistogram( new TH1F("k10p",";(x_{10}-median)/median; Lepton-jet assignments",100,-0.5,0) );
  controlHistos.addHistogram( new TH1F("k25p",";(x_{25}-median)/median; Lepton-jet assignments",100,-1.0,0) );
  controlHistos.addHistogram( new TH1F("k75p",";(x_{75}-median)/median; Lepton-jet assignments",100,0,0.5) );
  controlHistos.addHistogram( new TH1F("k90p",";(x_{90}-median)/median; Lepton-jet assignments",100,0,1.0) );

  TString cats[]={"ee","mumu","emu"};//,"etau","mutau"};
  size_t ncats=sizeof(cats)/sizeof(TString);
  TString subcats[]={"","eq0btags","eq1btags","geq2btags","zcands","ss"};
  size_t nsubcats=sizeof(subcats)/sizeof(TString);
  for(size_t icat=0; icat<ncats; icat++)
    for(size_t jcat=0; jcat<nsubcats; jcat++)
      controlHistos.initMonitorForStep(cats[icat]+subcats[jcat]);
  

  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += gSystem->BaseName(evurl);
  TFile *file=TFile::Open(outUrl, "recreate");

  //
  // TMVA
  //  
  TMVA::Factory *tmvaFactory=0;
  TMVA::Reader *tmvaReader=0;
  const unsigned int nVariables = varsList.size()+1;
  std::vector<Double_t> tmvaVarsD( nVariables,0 );
  std::vector<Float_t> tmvaVarsF( nVariables,0 );
  if(useMVA)
    {
      TMVA::Tools::Instance();
      
      if(trainMVA)
	{
	  //create the factory object. 
	  tmvaFactory = new TMVA::Factory( studyTag.Data(), file,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
	 
	  //variables to study
	  TString varsString("");
	  for(std::vector<std::string>::iterator it = varsList.begin(); it != varsList.end(); it++) 
	    {
	      if(it != varsList.begin()) varsString += ":";
	      varsString += *it;
	      tmvaFactory->AddVariable( *it, *it, "", 'F');
	    }
	  tmvaFactory->AddSpectator( "eventCategory" );
	  
	  cout << "==> Start TMVAClassification with " << methodList.size() << " methods and " << nVariables-1 << " variables" << endl;
	}
      else
	{
	  //reader for the methods
	  tmvaReader = new TMVA::Reader( "!Color:!Silent" );
	  
	  //variables to use
	  for(size_t ivar=0; ivar<varsList.size(); ivar++)   tmvaReader->AddVariable( varsList[ivar], &tmvaVarsF[ivar] );
	  tmvaReader->AddSpectator("eventCategory", &tmvaVarsF[varsList.size()]);
	  
	  //read the methods already trained
	  for(size_t imet=0; imet<methodList.size(); imet++)
	    {
	      //open the file with the method description         
	      TString weightFile = weightsDir + "/" + studyTag + "_" + methodList[imet] + ".weights.xml";
	      gSystem->ExpandPathName(weightFile);
	      tmvaReader->BookMVA(methodList[imet], weightFile);
	      TH1 *h=tmva::getHistogramForDiscriminator( methodList[imet] );
	      controlHistos.addHistogram( h );
	      controlHistos.addHistogram( (TH1 *)h->Clone(h->GetName()+TString("correct")) );
	      controlHistos.addHistogram( (TH1 *)h->Clone(h->GetName()+TString("wrong")) );
	      controlHistos.addHistogram( new TH1F(methodList[imet]+TString("decision"),";Good decisions",2,0.,2.) );
	    }
	}
    }
  
  
  //fix entries flag
  ofstream *outf=0;
  //if(!isMC) outf=new ofstream("highmassevents.txt",ios::app);
  
  //process events file
  gSystem->ExpandPathName(evurl);
  TFile *evfile = TFile::Open(evurl);
  if(evfile==0) return -1;
  if(evfile->IsZombie()) return -1;
  EventSummaryHandler evSummaryHandler;
  if( !evSummaryHandler.attachToTree( (TTree *)evfile->Get(dirname) ) ) 
    {
      evfile->Close();
      return -1;
    }  
  TTree *evTree=evSummaryHandler.getTree();

  //init event spy
  EventSummaryHandler *spyEvents=0;
  TFile *spyFile=0;
  TDirectory *spyDir=0;  
  if(saveSummaryTree)
    {
      spyEvents = new EventSummaryHandler;
      spyFile = TFile::Open("EventSummaries.root","UPDATE");
      TString evtag=gSystem->BaseName(evurl);
      evtag.ReplaceAll(".root","");
      spyFile->rmdir(evtag);
      spyDir = spyFile->mkdir(evtag);
      TTree *outT = new TTree("data","Event summary");
      spyEvents->initTree(outT);
    }

  //process kin file
  TString kinUrl(evurl);
  kinUrl.ReplaceAll(".root","/"+kindir);
  gSystem->ExpandPathName(kinUrl);
  cout << "Kin results from " << kinUrl << " to be processed with summary from " << evurl << endl;
  KinResultsHandler kinHandler;
  kinHandler.init(kinUrl,false);
  TChain *t=kinHandler.getResultsChain();
 
  //loop over events
  std::map<TString,int> selEvents;
  for (int inum=0; inum < evTree->GetEntriesFast(); ++inum)
  {
    evTree->GetEvent(inum);
    EventSummary_t &ev = evSummaryHandler.getEvent();
    TString key(""); key+= ev.run; key+="-"; key += ev.lumi; key+="-"; key += ev.event;
    selEvents[key]=inum;
  }


  FILE * mtop_text_file;
  //mtop_text_file = fopen("csv.mtop.el.data.outs.correctconverter","a+");
  mtop_text_file = fopen("csv.mtop.el.data.outs.correctconverter","w");


  //loop over kin results
  int nresults(0),neventsused(0);
  int nsigtrain(0), nsigtest(0), nbkgtrain(0), nbkgtest(0);
  for (int inum=0; inum < t->GetEntries(); ++inum)
    {
      t->GetEvent(inum);
      
      //get original event
      Int_t irun,ievent,ilumi;
      kinHandler.getEventInfo(irun,ievent,ilumi);
      
      TString key("");  key+= irun; key+="-"; key += ilumi;  key+="-"; key += ievent;
      
      if(selEvents.find(key)==selEvents.end()) continue;

      nresults++;
      evTree->GetEntry( selEvents[key] );

      //get event summary
      EventSummary_t &ev = evSummaryHandler.getEvent();
      
      //fill histos
      float weight = ev.weight;

      if(isMC)
	{
	  if(mcTruthMode==1 && !ev.isSignal) continue;
	  if(mcTruthMode==2 && ev.isSignal) continue;
	}
    
      std::vector<TString> categs;
      categs.push_back("all");
      if(ev.cat==MUMU)  categs.push_back("mumu");
      if(ev.cat==EE)  categs.push_back("ee");
      if(ev.cat==EMU) categs.push_back("emu");
      if(ev.cat==ETAU) categs.push_back("etau");
      if(ev.cat==MUTAU) categs.push_back("mutau");
    
      PhysicsEvent_t phys = getPhysicsEventFrom(ev);
      sort(phys.jets.begin(),phys.jets.end(),sortByBtag);
      sort(phys.leptons.begin(),phys.leptons.end(),sortByPt);

      int nRecoBs( (fabs(phys.jets[0].flavid)==5) + (fabs(phys.jets[1].flavid)==5) );
      if(!ev.isSignal) nRecoBs=0;
      int iCorrectComb=0;
      if(nRecoBs>1)
	{
	  //the charge of the generator level matched particles must be opposite
	  int assignCode=(phys.leptons[0].genid*phys.jets[0].flavid);
	  if(assignCode<0) iCorrectComb=1;
	  else             iCorrectComb=2;
	 //  cout << iCorrectComb << " | " 
	  // 	       << phys.leptons[0].genid << " (" << phys.leptons[0].pt() << ") "
	  // 	       << phys.leptons[1].genid << " (" << phys.leptons[1].pt() << ") |"
	  // 	       << phys.jets[0].flavid << " (" << phys.jets[0].btag1 << ") "
	  // 	       << phys.jets[1].flavid << " (" << phys.jets[1].btag1 << ") |"
	  // 	       << endl;
 	}
      
      //btag counting
      int nbtags(0);
      for(size_t ijet=0; ijet<phys.jets.size(); ijet++) nbtags += (phys.jets[ijet].btag1>1.7);
     
      //get the combination preferred by KIN
      TH1F *h1=kinHandler.getHisto("mt",1), *h2=kinHandler.getHisto("mt",2);
      h1->Rebin(2); h2->Rebin(2);
    
      std::map<TH1*, std::map<TString,Double_t> > histoVars;
      histoVars[h1]=histoAnalyzer.analyzeHistogram(h1);
      histoVars[h2]=histoAnalyzer.analyzeHistogram(h2);
      for(std::map<TH1 *,std::map<TString,Double_t> >::iterator hvit=histoVars.begin(); hvit !=histoVars.end(); hvit++)
	{
	  std::map<TString,Double_t> &res = hvit->second;
	  if(res["kIntegral"]>0)
	    {
	      for(std::map<TString,Double_t>::iterator resIt = res.begin(); resIt != res.end(); resIt++)
		controlHistos.fillHisto(resIt->first,"all",resIt->second,weight);
	    } 
	}
  
      //
      // standard solution counting
      //
      Int_t icomb=(h1->Integral()< h2->Integral())+1;
      TH1F *mpref=kinHandler.getHisto("mt",icomb);
      double mtop = kinHandler.getMPVEstimate(mpref) [1];
      TH1F *mttbarpref=kinHandler.getHisto("mttbar",icomb);
      double mttbar = kinHandler.getMPVEstimate(mttbarpref)[1];
      TH1F *afbpref=kinHandler.getHisto("afb",icomb);
      double afb = kinHandler.getMPVEstimate(afbpref)[1];

      fprintf(mtop_text_file,"%f\n", mtop);
    }

  fclose(mtop_text_file);

}  
