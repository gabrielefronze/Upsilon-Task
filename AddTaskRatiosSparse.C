#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonTrackCuts.h"
#include "AliAnalysisTaskRatiosSparse.h"
#endif

AliAnalysisTaskRatiosSparse* AddTaskRatiosSparse(Bool_t isMC = kFALSE, TString changeName = "")
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
  AliAnalysisTaskRatiosSparse *ratiosAnalysisTask = new AliAnalysisTaskRatiosSparse(taskName.Data(),muonCuts,isMC);
  if ( isMC ) ratiosAnalysisTask->SetTrigClassPatterns("ANY");
  mgr->AddTask(ratiosAnalysisTask);
   // Connect containers
   mgr->ConnectInput  (ratiosAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (ratiosAnalysisTask,  1, coutput1);

   return ratiosAnalysisTask;
}
