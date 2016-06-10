int runTaskTree (TString runMode="",TString analysisMode="",TString inputName="",TString inputOptions="",TString softVersions="",TString analysisOptions="",TString taskOptions=""){

	runMode ="terminate";
	analysisMode ="local";
	//inputName="/Users/Gabriele/cernbox/AlirootXcode/Upsilon/Task/AliAOD.Muons.root";
	inputName ="/Users/Gabriele/cernbox/AlirootXcode/Upsilon/Task/MC_upsilon_modOCDB.txt";
	//inputName ="Find;BasePath=/alice/cern.ch/user/g/gfronze/sim/LHC15o/upsilon/244982/;FileName=AliAOD.Muons.root;";
	//inputName ="Find;BasePath=/alice/data/2015/LHC15o/000246087/muon_calo_pass1/AOD/;FileName=AliAOD.Muons.root;";
	//inputName ="Find;BasePath=/alice/cern.ch/user/g/gfronze/sim/LHC15o/upsilon_chEffSyst/;FileName=AliAOD.Muons.root;";
	//inputName ="AliAOD.root"s
	inputOptions ="MC AOD";
	softVersions ="aliphysics=vAN-20160203-1";
	analysisOptions ="NOPHYSSEL";
	taskOptions ="";

  gROOT->LoadMacro(gSystem->ExpandPathName("$TASKDIR/runTaskUtilities.C"));

  //SetupAnalysis(runMode,analysisMode,inputName,inputOptions,softVersions,analysisOptions, "libPWGmuon.so AliAnalysisTaskUpsilonTreeAEff.cxx AddTaskUpsilonTreeAEff.C",". $ALICE_ROOT/include $ALICE_PHYSICS/include","UpsilonTaskOutputFullNoPDCAProof");
	SetupAnalysis(runMode,analysisMode,inputName,inputOptions,softVersions,analysisOptions, "libPWGmuon.so AliAnalysisTaskUpsilonTreeAEff.cxx AddTaskUpsilonTreeAEff.C",". $ALICE_ROOT/include $ALICE_PHYSICS/include","MCModOCDBAEff");

	if ( !analysisMode.Contains("local") ){
		AliAnalysisAlien* plugin = static_cast<AliAnalysisAlien*>(AliAnalysisManager::GetAnalysisManager()->GetGridHandler()); // Uncomment it if you want to configure the plugin...
		//plugin->SetGridWorkingDIr("/alice/data/2015/LHC15o");
		//plugin->SetMergeViaJDL(kTRUE);
		plugin->SetOverwriteMode(0);
	}

  Bool_t isMC = IsMC(inputOptions);
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

  AliAnalysisTaskUpsilonTreeAEff* taskMy = AddTaskUpsilonTreeAEff(isMC, "first");
  //taskMy->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB | AliVEvent::kINT7);
  //taskMy->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kMUON);

// AliMuonEventCuts* eventCuts = BuildMuonEventCuts(map); // Pre-configured AliMuonEventCuts
// SetupMuonBasedTask(task,eventCuts,taskOptions,map); // Automatically setup "task" if it derives from AliVAnalysisMuon

  StartAnalysis(runMode,analysisMode,inputName,inputOptions);

  return 0;
}
