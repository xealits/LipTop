#include <iostream>
#include <boost/shared_ptr.hpp>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;
using namespace top;

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  //check arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }
  
  const edm::ParameterSet &kinProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("kinProcess");
  TString url=kinProcess.getParameter<std::string>("input");
  TString output=kinProcess.getParameter<std::string>("output");
  int evStart=kinProcess.getParameter<int>("evStart");
  int evEnd=kinProcess.getParameter<int>("evEnd");
  TString dirname = kinProcess.getParameter<std::string>("dirName");
  TString scheme = kinProcess.getParameter<std::string>("kinScheme");
  int maxTries = kinProcess.getParameter<int>("maxTries");
  int maxJetMult = kinProcess.getParameter<int>("maxJetMult");
  
  float mw=kinProcess.getParameter<double>("mw");
  float mb=kinProcess.getParameter<double>("mb");

  TString etaFileName = kinProcess.getParameter<std::string>("etaResolFileName"); gSystem->ExpandPathName(etaFileName);
  JetResolution stdEtaResol(etaFileName.Data(),false);

  TString phiFileName = kinProcess.getParameter<std::string>("phiResolFileName"); gSystem->ExpandPathName(phiFileName);
  JetResolution stdPhiResol(phiFileName.Data(),false);

  TString ptFileName  = kinProcess.getParameter<std::string>("ptResolFileName");  gSystem->ExpandPathName(ptFileName);
  JetResolution stdPtResol(ptFileName.Data(),true); 
  
  //TString uncFile =  kinProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(uncFile);
  //JetCorrectionUncertainty jecUnc(uncFile.Data());

  TString jesUncSourcesUrl = kinProcess.getParameter<std::string>("jesUncFileName"); gSystem->ExpandPathName(jesUncSourcesUrl);

  // Instantiate uncertainty sources
  // FIX THE NAMES
  const int nsrc = 41;
  const char* srcnames[nsrc] =
    {"AbsoluteStat", "AbsoluteScale", "AbsoluteFlavMap", "AbsoluteMPFBias",
     "HighPtExtra", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "Time",
     "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB",
     "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeFSR",
     "RelativeStatEC2", "RelativeStatHF",
     "PileUpDataMC", "PileUpPtBB", "PileUpPtEC", "PileUpPtHF", "PileUpBias",
     "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalMC",
     "Total", "TotalNoFlavor",
     "FlavorZJet", "FlavorPhotonJet", "FlavorPureGluon", "FlavorPureQuark", "FlavorPureCharm", "FlavorPureBottom",
     "CorrelationGroupMPFInSitu", "CorrelationGroupIntercalibration", "CorrelationGroupbJES", "CorrelationGroupFlavor", "CorrelationGroupUncorrelated"};

     //{"Absolute", "HighPtExtra", "SinglePion", "Flavor", "Time",
     //"RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
     //"RelativeStatEC2", "RelativeStatHF", "RelativeFSR",
     //"PileUpDataMC", "PileUpOOT", "PileUpPt", "PileUpBias", "PileUpJetRate"};

  // Total JESunc for reference
  // CHECK THE NAME PileUpJetRate
  //JetCorrectionUncertainty *totalUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jesUncSourcesUrl.Data(), "SubTotalPileUp")));
  JetCorrectionUncertainty *totalUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jesUncSourcesUrl.Data(), "Total")));
  cout << jesUncSourcesUrl << endl;

  //open the file and get directory
  TFile *file = TFile::Open(url);
  if(file==0) return -1;
  if(file->IsZombie()) return -1;
  EventSummaryHandler evSummaryHandler;
  if( !evSummaryHandler.attachToTree( (TTree *)file->Get(dirname) ) ) 
    {
      file->Close();
      return -1;
    }
  
  //check run range
  if(evEnd<0 || evEnd>evSummaryHandler.getEntries() ) 
    evEnd=evSummaryHandler.getEntries();
  if(evStart > evEnd ) 
    {
      file->Close();
      return -1;
    }
  


  FILE * events_text_file;
  //events_text_file = fopen("csv.mtop.el.data.outs.correctconverter","a+");
  events_text_file = fopen("awk.spyntuple.events.particles.data","w");
  fprintf(events_text_file,"AWK style data: records of event data with fields storing data on each particle in the event (ID:px,py,pz,en).\nRecords (events) are separated with 1 blank linem fields (particles) are placed on separate lines.\nID codes: 0 - MET, 1 - jet, (+-)11/12/13 - (charge)electron/muon/tau.\n\n");

  for( int iev=evStart; iev<evEnd; iev++)
    {
      evSummaryHandler.getEntry(iev);
      EventSummary_t &ev = evSummaryHandler.getEvent();

      for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
         {
            fprintf(events_text_file,"%d:", ev.id[ipart]);
            fprintf(events_text_file,"%f,%f,%f,%f", ev.px[ipart], ev.py[ipart], ev.pz[ipart], ev.en[ipart]);
            fprintf(events_text_file,"\n");
         }
      fprintf(events_text_file,"\n");
    }

  fclose(events_text_file);
}
