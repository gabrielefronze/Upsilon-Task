#include <AliMuonAccEffSubmitter.h>

void LaunchSimulation(TString templateDir, TString generatorFile, TString outDir, TString runList, TString runMode){
  AliMuonAccEffSubmitter sub("GenParamCustom");

  sub.ShouldOverwriteFiles(true);
  sub.UseOCDBSnapshots(kFALSE);
  sub.SetCompactMode(11);

  sub.SetTemplateDir(templateDir.Data());
  sub.SetGenerator(generatorFile.Data());
  sub.SetRemoteDir(outDir.Data());
  sub.SetRunList(runList.Data());

  sub.SetAliPhysicsVersion("VO_ALICE@AliPhysics::v5-06-19-01");
  sub.SetVar("VAR_MUONMCMODE","3");
  sub.SetVar("VAR_REC_ALIGNDATA","alien://folder=/alice/simulation/2008/v4-15-Release/Residuals"); // IN REALTA' il TEMPLATE va commentato
  sub.SetVar("VAR_SIM_ALIGNDATA","alien://folder=/alice/simulation/2008/v4-15-Release/Full");

  sub.SetOCDBPath(const char *ocdbPath);

  Double_t fracOfEvents = 0.05;
  sub.MakeNofEventsPropToTriggerCount("CMUL7-B-NOPF-MUFAST",fracOfEvents);

  Int_t maxEventsPerChunk = ( runMode == "LOCALTEST" ) ? 10 : 1000;
  sub.SetMaxEventsPerChunk(maxEventsPerChunk);
  if ( runMode == "LOCALTEST" ) sub.MakeNofEventsFixed(sub.MaxEventsPerChunk());

  sub.Print();
  sub.Run(runMode.Data());
}
