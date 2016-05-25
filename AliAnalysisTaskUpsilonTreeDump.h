#ifndef AliAnalysisTaskUpsilonTreeDump_H
#define AliAnalysisTaskUpsilonTreeDump_H
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

class AliAnalysisTaskUpsilonTreeDump : public AliAnalysisTaskSE
{
  public:
    virtual void NotifyRun();
    
    AliAnalysisTaskUpsilonTreeDump();
    AliAnalysisTaskUpsilonTreeDump(const char *name, AliMuonTrackCuts *cuts);
    virtual ~AliAnalysisTaskUpsilonTreeDump();

    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *);
    //Bool_t   MatchTriggerDigitsNChambers (AliVParticle *track, Int_t nchambers) const;                    //  Muon track matches trigger digits in nchambers chambers
    //Bool_t   MatchTriggerDigitsNBending  (AliVParticle *track, Int_t nbending,Int_t nnotbending) const;   //  Muon track matches trigger digits in nbending bending sides and nnotbending notbending sides

  private:
    AliAODEvent* fAODEvent;
    AliESDEvent* fESDEvent;
    TList *fOutput; //!<TList output object
    AliMuonTrackCuts *fCuts;
    ULong64_t fNEvents;

    enum{
        kMomentum,
        kTransverse,
        kRapidity,
        kCentrality,
        kCharge,
        kSparseDimension
    };

    Float_t *fTreeData;

    enum{
        kPhysicsSelection,
        kTriggerSelection,
        kHasTrackerData,
        kLowPtCut,
        kEtaCut,
        kRAbsCut,
        kStandardMuonCuts,
        kDimension
    };

    Bool_t *fPassedCutsMuon;
  
  ClassDef(AliAnalysisTaskUpsilonTreeDump,2);
};

#endif
