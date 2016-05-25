#if !defined(__CINT__) || defined (__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TString.h"
#include "TSystem.h"
#include "TProof.h"
#include "TGrid.h"
#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TEnv.h"
#include "TObjArray.h"
#include "TObjString.h"

// STEER includes
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliCDBManager.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisDataContainer.h"

//parsing includes
#include <fstream>

#endif

using std::ifstream;

const char* const DELIMITER = "#";

Int_t ParseSettingFile(TString &runtype,TString &gridmode,bool &readAOD,bool &readMC,bool &embedding,bool &usePlugin,TString settingfilename){
	ifstream fin;
	fin.open(settingfilename.Data());
	if(!fin.good()) return -1;
	char buffer[256];

	fin.getline(buffer,256);
	runtype=TString(strtok(buffer, DELIMITER));
	if(runtype.Contains("gridterminate"))runtype="grid terminate";
	cout<<runtype<<endl;

	fin.getline(buffer,256);
	gridmode=TString(strtok(buffer, DELIMITER));
	cout<<gridmode<<endl;

	fin.getline(buffer,256);
	TString token;
	token=TString(strtok(buffer, DELIMITER));
	cout<<token<<endl;
	if(token.Contains("1")||token.Contains("true")||token.Contains("TRUE"))
		readAOD=true;
	else
		readAOD=false;
	cout<<readAOD<<endl;

	fin.getline(buffer,256);
	token=0x0;
	token=TString(strtok(buffer, DELIMITER));
	cout<<token<<endl;
	if((token.Contains("1")||token.Contains("true")||token.Contains("TRUE"))&&readAOD==false)
		readMC=true;
	else
		readMC=false;
	cout<<readMC<<endl;

	fin.getline(buffer,256);
	token=0x0;
	token=TString(strtok(buffer, DELIMITER));
	cout<<token<<endl;
	if(token.Contains("1")||token.Contains("true")||token.Contains("TRUE"))
		embedding=true;
	else
		embedding=false;
	cout<<embedding<<endl;

	fin.getline(buffer,256);
	token=0x0;
	token=TString(strtok(buffer, DELIMITER));
	cout<<token<<endl;
	if(token.Contains("1") || token.Contains("true") || token.Contains("TRUE"))
		usePlugin=true;
	else
		usePlugin=false;
	cout<<usePlugin<<endl;

	fin.close();
}

AliAnalysisGrid* CreateAlienHandler(const TString gridmode, const TString nRunFile);

void AliAnalysisTaskUpsilonLaunchGrid(const TString settingfilename, const TString nRunFile)
{
  TString runtype;                        // local, proof or grid //grid
  TString gridmode;                       // plugin mode ("full","test","offline","submit" or "terminate")
  const bool  readAOD;                    // 1 = AOD based analysis
  const bool  readMC;                     // 1 = Read MC
  const bool  embedding;                  // 1 = embedding
  const bool  usePlugin;                  // grid mode: use plugin

  ParseSettingFile(runtype,gridmode,readAOD,readMC,embedding,usePlugin,settingfilename);
    cout<<runtype<<endl;
    cout<<gridmode<<endl;
    cout<<readAOD<<endl;
    cout<<readMC<<endl;
    cout<<embedding<<endl;
    cout<<usePlugin<<endl;
    
  
  //load the required aliroot libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase"); //ok
  gSystem->Load("libESD"); //ok
  gSystem->Load("libAOD"); //ok
  gSystem->Load("libANALYSIS"); //ok
  gSystem->Load("libOADB.so"); //ok
  gSystem->Load("libANALYSISalice"); //ok
  gSystem->Load("libCORRFW.so"); //ok
  gSystem->Load("libPWGmuon.so"); //ok
  gSystem->Load("libPWGPPMUONlite.so"); //ok
  gSystem->Load("libPWGHFbase.so"); //ok
  
  gSystem->SetIncludePath("-I. -I$ALICE_PHYSICS/include ");
  //gROOT->LoadMacro("/home/gabriele/Tesi/Prova/AddTaskMTRchamberEfficiencyMOD.C");

  // check run type
  if(!runtype.Contains("local") && !runtype.Contains("proof")  && !runtype.Contains("grid") )
  {
    printf("Incorrect runtype! choose \"local\", \"prootf\", or \"grid\"\n");
    return;
  }
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("eff_task");
  
  TChain* analysisChain=0x0;

  // set run mode
  if(runtype.Contains("grid") && usePlugin)
  {
    AliAnalysisGrid *alienHandler = CreateAlienHandler(gridmode,nRunFile);
    if(!alienHandler)
      return;
    mgr->SetGridHandler(alienHandler);
    
    AliCDBManager::Instance()->SetDefaultStorage("alien://Folder=/alice/data/2011/OCDB");
  }
  else //local
  {
    if(readAOD)
    {
      analysisChain = new TChain("aodTree");
      if(readMC) //MC AOD
      {

      }
      
      else //real data AOD
      {

      }
    }
    else
    {
      analysisChain = new TChain("esdTree");
      if(readMC) //MC ESD
      {
      //         analysisChain->Add(Form("/home/marchiso/Desktop/AliESDs.root",nRun));
}
      
      else //real data ESD
      {
      	//analysisChain->Add("AliESDs_1.root");
      }  // local ESD RD
    }  // local ESD
  }  // local mode
  
  
//=============================================================================

  if(readAOD) //AOD
  {
    AliAODInputHandler *dataHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(dataHandler);
  }
  
  else //ESD
  {
    AliESDInputHandler *dataHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(dataHandler);
  }
  
	if ( readAOD==0 ) {
		gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
		AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(readMC==1);
		//physSelTask->GetPhysicsSelection()->SetPassName("passName");

		gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
		AliCentralitySelectionTask* centralityTask = AddTaskCentrality();
		if ( readMC==1 ) centralityTask->SetMCInput();
	}  

  // create the task
  gROOT->LoadMacro("./AliAnalysisTaskUpsilon.cxx++");
  AliAnalysisTaskUpsilon *task = new AliAnalysisTaskUpsilon("MyTask");

  if(readMC && !readAOD)
  {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
  }

  TString file=Form("%s:muonTrigger",mgr->GetCommonFileName()); 
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("list", TList::Class(),AliAnalysisManager::kOutputContainer, file.Data());

  mgr->ConnectInput(task,0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput);  

  if(!mgr->InitAnalysis())
    return;
  mgr->PrintStatus();
  

  mgr->StartAnalysis(runtype.Data(),analysisChain);
}




//=============================================================================
AliAnalysisGrid* CreateAlienHandler(const TString gridmode, const TString nRunFile)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  plugin->SetRunMode(gridmode.Data()); // AliEn plugin mode
  plugin->AddIncludePath("-I$ALICE_PHYSICS/include ");

  // JDL Merge
  if(gridmode.Contains("terminate")){
  	plugin->SetMergeViaJDL(kFALSE);
  	//plugin->SetOneStageMerging(kFALSE);  	
  }else{
  	plugin->SetMergeViaJDL(kTRUE);
  }
//plugin->SetMaxMergeStages(3);

  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20160203-1");

  /// Method 1: Create automatically XML collections using alien 'find' command.
  /*
  LHC11h 134 runs

  to do:
  170593, 170572, 170390, 170389, 170388, 170387, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170162, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170036, 170027, 169969, 169965, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169683, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169236, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168461, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168318, 168311, 168310, 168213, 168212, 168208, 168207, 168206, 168205, 168203, 168181, 168175, 168173, 168172, 168115, 168108, 168107, 168076, 168069, 168066, 167988, 167987, 167986, 167985, 167921, 167920, 167915, 167818, 
167814, 167813, 167808, 167807, 167806, 167713, 167706

  done:
  ---
  */
  
  plugin->SetGridDataDir("/alice/data/2015/LHC15o");
  plugin->SetDataPattern("pass2_muon/*/AliESDs.root"); //real data check reco pass and data base directory
  plugin->SetRunPrefix("000"); //real data
	ifstream fin;
	Int_t nRun=0;
	fin.open(nRunFile.Data());
	while(!fin.eof()){
		fin >> nRun;
		plugin->AddRunNumber(nRun);
		cout<<nRun<<endl;
	}
	fin.close();
  // plugin->AddRunNumber(nRun); //real data
  // plugin->AddRunNumber(170572); //real data
  // plugin->AddRunNumber(170593);
  //plugin->SetNtestFiles(10);
  plugin->SetDropToShell(kFALSE);

  /// Method 2: Declare existing data files (raw collections, xml collections, root file)
//   plugin->AddDataFile("Multiplicity_2290_a.xml");

  plugin->SetGridWorkingDir("MuonTrigger"); //mettere qui dentro . //path dell'output su alien (il mio)
//   plugin->CdWork(); //controlla se c'è la GridWorkingDir ed eventualmente la crea

  plugin->SetAnalysisSource("AliAnalysisTaskUpsilon.cxx");
//   plugin->SetAdditionalLibs("AliAnalysisTaskCS.h AliAnalysisTaskCS.cxx libTree.so libGeom.so libVMC.so libPhysics.so libSTEERBase.so libESD.so libAOD.so libANALYSIS.so libOADB.so libANALYSISalice.so libCORRFW.so libPWGHFbase.so libPWGmuon.so");
  plugin->SetAdditionalLibs("libPWGmuon.so AliAnalysisTaskUpsilon.h AliAnalysisTaskUpsilon.cxx");
  
  // Declare the output file names separated by blancs.
  //plugin->SetDefaultOutputs(kFALSE);
  plugin->SetOutputToRunNo();
  plugin->SetNumberOfReplicas(2);
  
  // Optionally define the files to be archived.
  // TString outputArchive = "log_archive.zip:stdout,stderr";
  // plugin->SetOutputArchive(outputArchive.Data());

  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("AnalysisMuonTrigger.C");
  
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);
  
  // Optionally set executable file (default: analysis.sh)
  //plugin->SetExecutable("CS.sh");
  
  // Optionally set time to live (default 30000 sec)
  //plugin->SetTTL(30000);
  
  // Optionally set input format (default xml-single)
  //plugin->SetInputFormat("xml-single");
  
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  //plugin->SetJDLName(Form("analysisCS_%d.jdl",nRun));
  
  // Optionally modify job price (default 1)
  //plugin->SetPrice(1);
  
  //plugin->SetUseSubmitPolicy();
  
  // Optionally modify split mode (default 'se')    
  //plugin->SetSplitMode("se");
  
  return plugin;
}