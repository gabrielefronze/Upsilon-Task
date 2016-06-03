#ifndef AliAnalysisTaskWeightedSpectrum_H
#define AliAnalysisTaskWeightedSpectrum_H
#include "AliAnalysisTaskSE.h"
//#include "AliMuonTrackCuts.h"

class AliAODEvent;
class AliESDEvent;
class TTree;
class AliAODDimuon;
// class AliAODTrack;
// class AliESDMuonTrack;
class AliVParticle;
class AliCFContainer;
class AliCFManager;
class AliMuonTrackCuts;
class AliMultSelection;
class TObjArray;

class AliAnalysisTaskWeightedSpectrum : public AliAnalysisTaskSE
{
  public:
    virtual void NotifyRun();

    AliAnalysisTaskWeightedSpectrum();
    AliAnalysisTaskWeightedSpectrum(AliMuonTrackCuts *cuts, Bool_t isMC, TString inputFileName, TString mode);
    virtual ~AliAnalysisTaskWeightedSpectrum();

    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *);
    //Bool_t   MatchTriggerDigitsNChambers (AliVParticle *track, Int_t nchambers) const;                    //  Muon track matches trigger digits in nchambers chambers
    //Bool_t   MatchTriggerDigitsNBending  (AliVParticle *track, Int_t nbending,Int_t nnotbending) const;   //  Muon track matches trigger digits in nbending bending sides and nnotbending notbending sides

  private:
    TList *fOutput; //!<TList output object
    AliMuonTrackCuts *fCuts;
    ULong64_t fNEvents;
    Bool_t fMode;
    TAxis *fRapidityAxis;
    TObjArray *fInputResponseFunctionsMC;
    TObjArray *fInputResponseFunctionsData;

    enum{
        kTransverse,
        kRapidity,
        kMass,
    };

    Bool_t fIsMC;

  ClassDef(AliAnalysisTaskWeightedSpectrum,2);
};

#endif
