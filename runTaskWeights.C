int runTaskWeights (TString runMode="",TString analysisMode="",TString inputName="",TString inputOptions="",TString softVersions="",TString analysisOptions="",TString taskOptions=""){

	runMode ="terminate";

	analysisMode ="local";
	//inputName ="Find;BasePath=/alice/cern.ch/user/a/alardeux/Sim/LHC15o/RealisticUps/Geant3_wVtx/MergedAOD/137runs/;FileName=AliAOD.Muons.root;";
	//inputName ="2015Ho_good_runs_175.txt";
	//inputName ="Find;BasePath=/alice/data/2015/LHC15o/000245766/muon_calo_pass1/AOD175/;FileName=AliAOD.root;";
	inputName ="/Users/Gabriele/AliAOD.Muons.root";
	inputOptions ="MC";
	softVersions ="aliphysics=vAN-20160510-1";
	(analysisMode.Contains("local")) ? analysisOptions ="NOPHYSSEL" : analysisOptions ="CENTR NOPHYSSEL";
	//analysisOptions="CENTR";
	taskOptions ="";

  gROOT->LoadMacro(gSystem->ExpandPathName("$TASKDIR/runTaskUtilities.C"));

  SetupAnalysis(runMode,analysisMode,inputName,inputOptions,softVersions,analysisOptions, "libPWGmuon.so AliAnalysisTaskWeightedSpectrum.cxx AddTaskRatiosWeights.C",". $ALICE_ROOT/include $ALICE_PHYSICS/include","UpsilonTaskOutputWeightedSpectrum");

	AliAnalysisAlien* plugin = static_cast<AliAnalysisAlien*>(AliAnalysisManager::GetAnalysisManager()->GetGridHandler()); // Uncomment it if you want to configure the plugin...
	//plugin->SetGridWorkingDir("analysis/LHC15o/UpsilonTaskOutputRatiosSparse");
	//plugin->SetMergeViaJDL(kFALSE);
	//plugin->SetRunNumber("246994");
	//plugin->SetRunFromPath(246994);
//plugin->SetOverwriteMode(1);
	if(!analysisMode.Contains("local"))plugin->SetOverwriteMode(1);

  Bool_t isMC = IsMC(inputOptions);
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

	TString inputFileName="/Users/Gabriele/cernbox/AlirootXcode/Upsilon/Task/full_analysis_outputs/FullRealData/AptOverLptHistograms.root";
	TString mode="F";

	AliAnalysisTaskWeightedSpectrum* taskMy = AddTaskRatiosWeights(isMC,inputFileName,mode,"first");

	if ( !isMC ){
	  taskMy->SelectCollisionCandidates(AliVEvent::kMUSPB);
	  //taskMy->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kMUON);
	}

// AliMuonEventCuts* eventCuts = BuildMuonEventCuts(map); // Pre-configured AliMuonEventCuts
// SetupMuonBasedTask(task,eventCuts,taskOptions,map); // Automatically setup "task" if it derives from AliVAnalysisMuon

  StartAnalysis(runMode,analysisMode,inputName,inputOptions);

  return 0;
}
