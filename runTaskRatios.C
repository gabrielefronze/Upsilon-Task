int runTaskRatios (TString runMode="",TString analysisMode="",TString inputName="",TString inputOptions="",TString softVersions="",TString analysisOptions="",TString taskOptions=""){

	runMode ="terminate";

	analysisMode ="local";
	//inputName ="Find;BasePath=/alice/cern.ch/user/a/alardeux/Sim/LHC15o/RealisticUps/Geant3_wVtx/MergedAOD/137runs/;FileName=AliAOD.Muons.root;";
	//inputName ="2015Ho_good_runs_175.txt";
	//inputName ="Find;BasePath=/alice/data/2015/LHC15o/000245766/muon_calo_pass1/AOD175/;FileName=AliAOD.root;";
	inputName ="/Users/Gabriele/AliAOD.Muons.root";
	inputOptions ="MC";
	softVersions ="aliphysics=vAN-20160510-1";
	(analysisMode.Contains("local")) ? analysisOptions ="" : analysisOptions ="CENTR";
	//analysisOptions="CENTR";
	taskOptions ="";

  gROOT->LoadMacro(gSystem->ExpandPathName("$TASKDIR/runTaskUtilities.C"));

  SetupAnalysis(runMode,analysisMode,inputName,inputOptions,softVersions,analysisOptions, "libPWGmuon.so AliAnalysisTaskRatiosSparse.cxx AddTaskRatiosSparse.C",". $ALICE_ROOT/include $ALICE_PHYSICS/include","UpsilonTaskOutputRatiosSparse");

	AliAnalysisAlien* plugin = static_cast<AliAnalysisAlien*>(AliAnalysisManager::GetAnalysisManager()->GetGridHandler()); // Uncomment it if you want to configure the plugin...
	//plugin->SetGridWorkingDir("UpsiOutRatiosSparseSim");
	//plugin->SetMergeViaJDL(kFALSE);
	//plugin->SetRunNumber("246994");
	//plugin->SetRunFromPath(246994);

	if(!analysisMode.Contains("local"))plugin->SetOverwriteMode(0);

  Bool_t isMC = IsMC(inputOptions);
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

	AliAnalysisTaskRatiosSparse* taskMy = AddTaskRatiosSparse(isMC, "first");

	if ( !isMC ){
	  taskMy->SelectCollisionCandidates(AliVEvent::kMUSPB);
	  //taskMy->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kMUON);
	}

// AliMuonEventCuts* eventCuts = BuildMuonEventCuts(map); // Pre-configured AliMuonEventCuts
// SetupMuonBasedTask(task,eventCuts,taskOptions,map); // Automatically setup "task" if it derives from AliVAnalysisMuon

  StartAnalysis(runMode,analysisMode,inputName,inputOptions);

  return 0;
}
