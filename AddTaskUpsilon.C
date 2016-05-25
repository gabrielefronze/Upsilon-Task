#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "AliAnalysisTaskUpsilon.h"
#endif

AliAnalysisTaskUpsilon* AddTaskUpsilon(Bool_t isMC = kFALSE, TString changeName = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskUpsilon", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddtaskUpsilon", "Upsilon task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":Upsilon" + changeName;
  else outputfile = "UpsilonAnalysis" + changeName + ".root";

  TString containerName = "UpsilonOut" + changeName;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  TString cutsName = "StdMuonPairCuts" + changeName;
  AliMuonPairCuts* muonPairCuts = new AliMuonPairCuts(cutsName.Data(), cutsName.Data());
  muonPairCuts->SetIsMC(isMC);


  // Create task
  TString taskName = "UpsilonTask" + changeName;
  AliMuonTrackCuts *muonCuts = new AliMuonTrackCuts("cuts","cuts");
  AliAnalysisTaskUpsilon *upsilonAnalysisTask = new AliAnalysisTaskUpsilon(taskName.Data(),muonCuts);
  if ( isMC ) upsilonAnalysisTask->SetTrigClassPatterns("ANY");
  mgr->AddTask(upsilonAnalysisTask);
   // Connect containers
   mgr->ConnectInput  (upsilonAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (upsilonAnalysisTask,  1, coutput1);

   return upsilonAnalysisTask;
}
