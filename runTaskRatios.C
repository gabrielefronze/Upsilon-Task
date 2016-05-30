int runTaskRatios (TString runMode="",TString analysisMode="",TString inputName="",TString inputOptions="",TString softVersions="",TString analysisOptions="",TString taskOptions=""){

	runMode ="full";
	analysisMode ="grid";
	inputName ="2015Ho_good_runs_175.txt";
	//inputName ="Find;BasePath=/alice/data/2015/LHC15o/000245766/muon_calo_pass1/AOD175/;FileName=AliAOD.root;";
	//inputName ="AliAOD.root"s
	inputOptions ="AOD";
	softVersions ="aliphysics=vAN-20160510-1";
	(analysisMode.Contains("local")) ? analysisOptions ="" : analysisOptions ="CENTR";
	//analysisOptions="CENTR";
	taskOptions ="";

  gROOT->LoadMacro(gSystem->ExpandPathName("$TASKDIR/runTaskUtilities.C"));

  SetupAnalysis(runMode,analysisMode,inputName,inputOptions,softVersions,analysisOptions, "libPWGmuon.so AliAnalysisTaskRatiosSparse.cxx AddTaskRatiosSparse.C",". $ALICE_ROOT/include $ALICE_PHYSICS/include","UpsilonTaskOutputRatiosSparse");

	AliAnalysisAlien* plugin = static_cast<AliAnalysisAlien*>(AliAnalysisManager::GetAnalysisManager()->GetGridHandler()); // Uncomment it if you want to configure the plugin...
	//plugin->SetGridWorkingDIr("/alice/data/2015/LHC15o");
	//plugin->SetMergeViaJDL(kTRUE);
	if(!analysisMode.Contains("local"))plugin->SetOverwriteMode(0);

  Bool_t isMC = IsMC(inputOptions);
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

  AliAnalysisTaskRatiosSparse* taskMy = AddTaskRatiosSparse(isMC, "first");
  taskMy->SelectCollisionCandidates(AliVEvent::kMUSPB);
  //taskMy->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kMUON);

// AliMuonEventCuts* eventCuts = BuildMuonEventCuts(map); // Pre-configured AliMuonEventCuts
// SetupMuonBasedTask(task,eventCuts,taskOptions,map); // Automatically setup "task" if it derives from AliVAnalysisMuon

  StartAnalysis(runMode,analysisMode,inputName,inputOptions);

  return 0;
}
