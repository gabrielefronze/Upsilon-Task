#ifndef AliAnalysisTaskRatiosSparse_H
#define AliAnalysisTaskRatiosSparse_H
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

class AliAnalysisTaskRatiosSparse : public AliAnalysisTaskSE
{
  public:
    virtual void NotifyRun();

    AliAnalysisTaskRatiosSparse();
    AliAnalysisTaskRatiosSparse(const char *name, AliMuonTrackCuts *cuts);
    virtual ~AliAnalysisTaskRatiosSparse();

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
      kTriggerFlag,
      kPt,
      kRapidity,
      kCentrality,
      kLocalBoard,
      kSparseDimension
    };

  ClassDef(AliAnalysisTaskRatiosSparse,1);
};

#endif
