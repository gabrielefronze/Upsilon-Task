#ifndef AliAnalysisTaskUpsilonTree2_H
#define AliAnalysisTaskUpsilonTree2_H
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

class AliAnalysisTaskUpsilonTree2 : public AliAnalysisTaskSE
{
  public:
    virtual void NotifyRun();

    AliAnalysisTaskUpsilonTree2();
    AliAnalysisTaskUpsilonTree2(const char *name, AliMuonTrackCuts *cuts, Bool_t isMC, Bool_t fillTree);
    virtual ~AliAnalysisTaskUpsilonTree2();

    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *);
    //Bool_t   MatchTriggerDigitsNChambers (AliVParticle *track, Int_t nchambers) const;                    //  Muon track matches trigger digits in nchambers chambers
    //Bool_t   MatchTriggerDigitsNBending  (AliVParticle *track, Int_t nbending,Int_t nnotbending) const;   //  Muon track matches trigger digits in nbending bending sides and nnotbending notbending sides

  private:
    TList *fOutput; //!<TList output object
    AliMuonTrackCuts *fCuts;
    ULong64_t fNEvents;

    enum{
        kMomentum,
        kTransverse,
        kRapidity,
        kCentrality,
        kLowMomentumMuonpt,
        kHighMomentumMuonpt,
        kMass,
        kSparseDimension
    };

    Bool_t fFillTree;
    Bool_t fIsMC;

  ClassDef(AliAnalysisTaskUpsilonTree2,4);
};

#endif
