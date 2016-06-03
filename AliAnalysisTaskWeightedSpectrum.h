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

    AliAnalysisTaskWeightedSpectrum(Bool_t isMC, TString inputFileName, TString mode);
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
    TString fMode;
    TAxis *fRapidityAxis;
    TObjArray *fInputResponseFunctions;

    enum{
        kMomentum,
        kTransverse,
        kRapidity,
        kCentrality,
        kPt1,
        kPt2,
        kMass,
        kWeight
    };

    Float_t *fTreeData;
    Bool_t fIsMC;

  ClassDef(AliAnalysisTaskWeightedSpectrum,2);
};

#endif
