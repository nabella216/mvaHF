#ifndef __PHOD_DTREE__
#define __PHOD_DTREE__

#include <TTree.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2F.h>

#include <vector>
#include <map>
#include <string>
#include <type_traits>

#include "xjjcuti.h"

#define MAX_XB       20000


enum class VariedCutType {
  Dtrk1Pt,
  Dtrk2Pt,
  DsvpvSig,
  Dchi2cl,
  Dalpha,
  Ddtheta
};

namespace hfupc
{

  class dtree
  {
  public:
    dtree(TTree* nt, bool ishi);
    dtree(TFile* outf, std::string name, bool ishi);
    TTree* nt() { return nt_; }

    // read
    template<typename T> T val(std::string br, int j);
    float* &operator[](std::string br) {//PROTECT AGAINST MISSING BRANCH BY DUMPING OUTPUT ON ACCESS
      if(!bvs_[br]) { std::cout<<__PRETTY_FUNCTION__<<" error: bad branch ["<<br<<"]."<<std::endl;}
      return bvf_[br];
    }
    bool status(std::string br) { return bvs_[br]; }
    int GetEntries() { return nt_->GetEntries(); }
    void GetEntry(int i) { nt_->GetEntry(i); }

    // tools
    int Dsize() { return Dsize_; }
    bool presel(int j);
    bool tightsel(int j);
    void varySelection(std::vector<VariedCutType>);
    void fillDcut();

    // fill new file tree
    void ClearDsize() { if(newtree_) { Dsize_ = 0; } }
    void Fillall(dtree* nt, int j);
    void Fillone(std::string br, float val) { if(newtree_) { bvf_[br][Dsize_] = val; } }
    void Dsizepp() { if(newtree_) Dsize_++; }
    void Fill() { if(newtree_ ){ dr_->cd(); nt_->Fill(); } }
  private:
    //Members set in constructor
    TTree* nt_;
    bool ishi_;

    TDirectory* dr_;
    bool newtree_;
    void setbranchaddress();
    bool cutProtectionMsg(std::string funcStr);
    void branch();
    //List of non-array variables, all datatypes
    std::vector<std::string> tbfio_ =
      {
	"Dsize"
      };
    //TTree variables, non-array, all datatypes
    int Dsize_;

    //List of variables, arrays tied to Dsize -> floats, ints, bools
    std::vector<std::string> tbvf_ =
      {
	"Dmass",
	"Dpt",
	"Deta",
	"Dphi",
	"Dy",
	"Dchi2cl",
	"Ddtheta",
	"Dalpha",
	"DsvpvDistance",
	"DsvpvDisErr",
	"DlxyBS",
	"DlxyBSErr",
	"DdthetaBS",
	"DdthetaBScorr",
	"DsvpvDistance_2D",
	"DsvpvDisErr_2D",
	"Dtrk1Pt",
	"Dtrk2Pt",
	"Dtrk1PtErr",
	"Dtrk2PtErr",
	"Dtrk1Eta",
	"Dtrk2Eta",
	"Dtrk1Dz1",
	"Dtrk2Dz1",
	"Dtrk1DzError1",
	"Dtrk2DzError1",
	"Dtrk1Dxy1",
	"Dtrk2Dxy1",
	"Dtrk1DxyError1",
	"Dtrk2DxyError1",
	"Dtrk1PixelHit",
	"Dtrk2PixelHit",
	"Dtrk1StripHit",
	"Dtrk2StripHit",
	"Dtrk1nStripLayer",
	"Dtrk2nStripLayer",
	"Dtrk1nPixelLayer",
	"Dtrk2nPixelLayer",
	"Dtrk1Chi2ndf",
	"Dtrk2Chi2ndf",
	"BDT",
	"Dgenpt",
	"Dgeneta",
	"Dgenphi",
	"Dgeny"
      };
    std::vector<std::string> tbvi_ =
      {
	"Dgen",
	"DgencollisionId",
      };
    std::vector<std::string> tbvo_ =
      {
	"Dtrk1highPurity",
	"Dtrk2highPurity"
      };

    std::map<std::string, float*> bvf_;
    std::map<std::string, int*> bvi_;
    std::map<std::string, bool*> bvo_;

    std::map<std::string, bool> bvs_; //

    //Vectors defined for function protections
    //If you define a function that is a complex combination of variables, add this internal protection
    //Default bool to true, if single var not found switch to false
    //For bool presel(int j)
    bool preselVarsFound_ = true;
    std::vector<std::string> preselVarHI_ =
      {
	"Dtrk1Eta",
	"Dtrk1highPurity",
	"Dtrk1Pt",
	"Dtrk1PtErr",
	"Dtrk1PixelHit",
	"Dtrk1StripHit",
	"Dtrk1Chi2ndf",
	"Dtrk1nStripLayer",
	"Dtrk1nPixelLayer",
	"Dtrk2Eta",
	"Dtrk2highPurity",
	"Dtrk2Pt",
	"Dtrk2PtErr",
	"Dtrk2PixelHit",
	"Dtrk2StripHit",
	"Dtrk2Chi2ndf",
	"Dtrk2nStripLayer",
	"Dtrk2nPixelLayer"
      };
    std::vector<std::string> preselVarPP_ =
      {
	"Dtrk1Eta",
	"Dtrk1highPurity",
	"Dtrk1Pt",
	"Dtrk1PtErr",
	"Dtrk2Eta",
	"Dtrk2highPurity",
	"Dtrk2Pt",
	"Dtrk2PtErr",
      };
    //For bool tightsel(int j)
    bool tightselVarsFound_ = true;
    std::vector<std::string> tightselVarHI_ =
      {
      };
    std::vector<std::string> tightselVarPP_ =
      {
	"DsvpvDistance",
	"DsvpvDisErr",
	"Dalpha",
	"Ddtheta",
	"Dchi2cl"
      };

    // Cut value based on pT and y
    static const int nPt = 3;
    static const int nY = 6;
    float ptBins[nPt + 1] = {2.0, 4.0, 8.0, 999.0};
    float yBins[nY + 1] = {-2.4, -2.0, -1.0, 0.0, 1.0, 2.0, 2.4};
    TH2F DsvpvSigCut, DalphaCut, DdthetaCut, Dchi2clCut;
    TH2F Dtrk1PtCut, Dtrk2PtCut;
    using floatVector2D = std::vector<std::vector<float> >;

    // Nominal cut values
    floatVector2D DsvpvSigCutValue = {
      {3.5, 3.5, 3.5, 3.5, 3.5, 3.5},
      {3.5, 3.5, 3.5, 3.5, 3.5, 3.5},
      {3.5, 3.5, 3.5, 3.5, 3.5, 3.5}
    };
    floatVector2D Dchi2clCutValue = {
      {0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1, 0.1, 0.1, 0.1}
    };
    floatVector2D DalphaCutValue = {
      {0.2, 0.2, 0.4, 0.4, 0.2, 0.2},
      {0.25, 0.25, 0.35, 0.35, 0.25, 0.25},
      {0.25, 0.25, 0.4, 0.4, 0.25, 0.25}
    };
    floatVector2D DdthetaCutValue = {
      {0.3, 0.3, 0.5, 0.5, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3}
    };
    floatVector2D Dtrk1PtCutValue = {
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
    };
    floatVector2D Dtrk2PtCutValue = {
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
    };

    // Dfinder or looser cut values
    floatVector2D VariedDsvpvSigCutValue = {
      {2.5, 2.5, 2.5, 2.5, 2.5, 2.5},
      {2.5, 2.5, 2.5, 2.5, 2.5, 2.5},
      {2.5, 2.5, 2.5, 2.5, 2.5, 2.5}
    };
    floatVector2D VariedDchi2clCutValue = {
      {0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
      {0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
      {0.05, 0.05, 0.05, 0.05, 0.05, 0.05}
    };
    floatVector2D VariedDalphaCutValue = {
      {5.0, 5.0, 5.0, 5.0, 5.0, 5.0},
      {5.0, 5.0, 5.0, 5.0, 5.0, 5.0},
      {5.0, 5.0, 5.0, 5.0, 5.0, 5.0}
    };
    floatVector2D VariedDdthetaCutValue = {
      {5.0, 5.0, 5.0, 5.0, 5.0, 5.0},
      {5.0, 5.0, 5.0, 5.0, 5.0, 5.0},
      {5.0, 5.0, 5.0, 5.0, 5.0, 5.0}
    };
    floatVector2D VariedDtrk1PtCutValue = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };
    floatVector2D VariedDtrk2PtCutValue = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };

    std::vector<floatVector2D*> cutValues = {
      &Dtrk1PtCutValue,
      &Dtrk2PtCutValue,
      &DsvpvSigCutValue,
      &Dchi2clCutValue,
      &DalphaCutValue,
      &DdthetaCutValue
    };
    std::vector<floatVector2D> variedCutValues = {
      VariedDtrk1PtCutValue,
      VariedDtrk2PtCutValue,
      VariedDsvpvSigCutValue,
      VariedDchi2clCutValue,
      VariedDalphaCutValue,
      VariedDdthetaCutValue
    };

  };
}

hfupc::dtree::dtree(TTree* nt, bool ishi) :
  nt_(nt), ishi_(ishi),
  Dtrk1PtCut("Dtrk1PtCut", "Dtrk1PtCut", nPt, ptBins, nY, yBins),
  Dtrk2PtCut("Dtrk2PtCut", "Dtrk2PtCut", nPt, ptBins, nY, yBins),
  DsvpvSigCut("DsvpvSigCut", "DsvpvSigCut", nPt, ptBins, nY, yBins),
  DalphaCut("DalphaCut", "DalphaCut", nPt, ptBins, nY, yBins),
  DdthetaCut("DdthetaCut", "DdthetaCut", nPt, ptBins, nY, yBins),
  Dchi2clCut("Dchi2clCut", "Dchi2clCut", nPt, ptBins, nY, yBins)
{
  std::cout<<"\e[32;1m -- "<<__PRETTY_FUNCTION__<<"\e[0m"<<std::endl;
  fillDcut();

  for(auto& i : tbvf_) { bvf_[i] = new float[MAX_XB]; }
  for(auto& i : tbvi_) { bvi_[i] = new int[MAX_XB]; }
  for(auto& i : tbvo_) { bvo_[i] = new bool[MAX_XB]; }

  // newtree_ = !nt_->FindBranch("Dsize");
  newtree_ = false;

  nt_->SetBranchStatus("*", 0);

  //Check all defined branches exist, and set status
  std::vector<std::vector<std::string> > tbSuperVect = {tbfio_, tbvf_, tbvi_, tbvo_};
  for(unsigned int i = 0; i < tbSuperVect.size(); ++i){
    for(auto& b : tbSuperVect[i]){ bvs_[b] = false;
      if(nt_->FindBranch(b.c_str())) { nt_->SetBranchStatus(b.c_str(), 1); bvs_[b] = true; }
    }
  }

  setbranchaddress();

  //Check branches for presence of selections for complex cuts
  //Do it with two maps of pointers, each with strings as keys, since addresses are unique
  std::map<std::string, std::vector<std::string>*> functionStrToVarVect;
  std::map<std::string, bool*> functionStrToBool;
  if(ishi_){//init hi varvect maps
    functionStrToVarVect["presel"] = &preselVarHI_;
    functionStrToVarVect["tightsel"] = &tightselVarHI_;
  }
  else{//init pp varvect maps
    functionStrToVarVect["presel"] = &preselVarPP_;
    functionStrToVarVect["tightsel"] = &tightselVarPP_;
  }

  //Init bool maps (no hi dependency)
  functionStrToBool["presel"] = &preselVarsFound_;
  functionStrToBool["tightsel"] = &tightselVarsFound_;

  //Iterate over the map
  for(auto const & val : functionStrToVarVect){
    std::vector<std::string>* vect_p = val.second;
    for(unsigned int vI = 0; vI < vect_p->size(); ++vI){
      //std::cout << " " << vect_p->at(vI) << std::endl;//For Debug

      if(!bvs_[vect_p->at(vI)]){
	std::cout << __PRETTY_FUNCTION__ << ": VAR '" << vect_p->at(vI) << "' IS NOT FOUND IN GIVEN TTREE. " << std::endl;
	std::cout << " Function '" << val.first << "' will not work. setting false" << std::endl;
	*(functionStrToBool[val.first]) = false;
      }
    }
  }
}

hfupc::dtree::dtree(TFile* outf, std::string name, bool ishi) : ishi_(ishi)
{
  std::cout<<"\e[32;1m -- "<<__PRETTY_FUNCTION__<<"\e[0m"<<std::endl;

  for(auto& i : tbvf_) { bvf_[i] = new float[MAX_XB]; bvs_[i] = false; }
  for(auto& i : tbvi_) { bvi_[i] = new int[MAX_XB]; bvs_[i] = false; }
  for(auto& i : tbvo_) { bvo_[i] = new bool[MAX_XB]; bvs_[i] = false; }

  newtree_ = true;

  std::vector<std::string> p = xjjc::str_divide(name, "/");
  if(p.size() > 1)
    {
      if(!outf->cd(p[0].c_str())) dr_ = outf->mkdir(p[0].c_str());
      else { dr_ = outf->GetDirectory(p[0].c_str()); }
    }
  else { dr_ = outf; }
  dr_->cd();
  nt_ = new TTree(p.back().c_str(), "");
  branch();
  outf->cd();
}

void hfupc::dtree::branch()
{
  nt_->Branch("Dsize", &Dsize_, "Dsize/I");
  for(auto& b : tbvf_) nt_->Branch(b.c_str(), bvf_[b], Form("%s[Dsize]/F", b.c_str()));
  for(auto& b : tbvi_) nt_->Branch(b.c_str(), bvi_[b], Form("%s[Dsize]/I", b.c_str()));
  for(auto& b : tbvo_) nt_->Branch(b.c_str(), bvo_[b], Form("%s[Dsize]/O", b.c_str()));
}

void hfupc::dtree::setbranchaddress()
{
  if(bvs_["Dsize"]) nt_->SetBranchAddress("Dsize", &Dsize_);
  for(auto& b : tbvf_) { if(bvs_[b]) { nt_->SetBranchAddress(b.c_str(), bvf_[b]); } }
  for(auto& b : tbvi_) { if(bvs_[b]) { nt_->SetBranchAddress(b.c_str(), bvi_[b]); } }
  for(auto& b : tbvo_) { if(bvs_[b]) { nt_->SetBranchAddress(b.c_str(), bvo_[b]); } }
}

bool hfupc::dtree::cutProtectionMsg(std::string funcStr)
{
  std::cout << funcStr << ": ERROR! Not all required variables in TTree! PLEASE FIX. return false" << std::endl;
  return false;
}

template<typename T> T hfupc::dtree::val(std::string br, int j)
{
  if(!bvs_[br]) { std::cout<<__PRETTY_FUNCTION__<<" error: bad branch ["<<br<<"]."<<std::endl;}
  if(std::is_same<T, float>::value) { return bvf_[br][j]; }
  if(std::is_same<T, int>::value) { return bvi_[br][j]; }
  if(std::is_same<T, bool>::value) { return bvo_[br][j]; }
  return (T)0;
}

void hfupc::dtree::Fillall(dtree* nt, int j)
{
  if(!newtree_) return;
  for(auto& b : tbvf_) { if(nt->status(b)) { bvf_[b][Dsize_] = nt->val<float>(b, j); } }
  for(auto& b : tbvi_) { if(nt->status(b)) { bvi_[b][Dsize_] = nt->val<int>(b, j); } }
  for(auto& b : tbvo_) { if(nt->status(b)) { bvo_[b][Dsize_] = nt->val<bool>(b, j); } }
}

bool hfupc::dtree::presel(int j)
{
  //IF YOU CHANGE VARIABLES USED BY THIS FUNCTION,
  //PLEASE EDIT preselVarHI_, preselVarPP_ APPROPRIATELY!
  //Protections for cut
  if(!preselVarsFound_) return cutProtectionMsg(__PRETTY_FUNCTION__);
  float Dpt = bvf_["Dpt"][j];
  float Dy = bvf_["Dy"][j];
  int bin = DalphaCut.FindBin(Dpt, Dy);

  bool trkcut = fabs(bvf_["Dtrk1Eta"][j]) < 2.4 && fabs(bvf_["Dtrk2Eta"][j]) < 2.4 &&
                                                                               bvo_["Dtrk1highPurity"][j] && bvo_["Dtrk2highPurity"][j];
  if(ishi_ && trkcut &&
     bvf_["Dtrk1Pt"][j] > 1.0 && bvf_["Dtrk2Pt"][j] > 1.0 &&
     fabs(bvf_["Dtrk1PtErr"][j]/bvf_["Dtrk1Pt"][j]) < 0.1 && fabs(bvf_["Dtrk2PtErr"][j]/bvf_["Dtrk2Pt"][j]) < 0.1 &&
     (bvf_["Dtrk1PixelHit"][j]+bvf_["Dtrk1StripHit"][j]) >= 11 && (bvf_["Dtrk2PixelHit"][j]+bvf_["Dtrk2StripHit"][j]) >= 11 &&
     (bvf_["Dtrk1Chi2ndf"][j]/(bvf_["Dtrk1nStripLayer"][j]+bvf_["Dtrk1nPixelLayer"][j])) < 0.18 && (bvf_["Dtrk2Chi2ndf"][j]/(bvf_["Dtrk2nStripLayer"][j]+bvf_["Dtrk2nPixelLayer"][j])) < 0.18
     ) return true;
  if(!ishi_ && trkcut &&
     bvf_["Dtrk1Pt"][j] > Dtrk1PtCut[bin] &&
     bvf_["Dtrk2Pt"][j] > Dtrk2PtCut[bin] &&
     fabs(bvf_["Dtrk1PtErr"][j]/bvf_["Dtrk1Pt"][j]) < 0.1 && fabs(bvf_["Dtrk2PtErr"][j]/bvf_["Dtrk2Pt"][j]) < 0.1
     ) return true;
  return false;
}

void hfupc::dtree::fillDcut() {
  for (int iPt = 0; iPt < nPt; ++iPt) {
    for (int iY = 0; iY < nY; ++iY) {
      DsvpvSigCut.SetBinContent(iPt+1, iY+1, DsvpvSigCutValue[iPt][iY]);
      DalphaCut.SetBinContent(iPt+1, iY+1, DalphaCutValue[iPt][iY]);
      DdthetaCut.SetBinContent(iPt+1, iY+1, DdthetaCutValue[iPt][iY]);
      Dchi2clCut.SetBinContent(iPt+1, iY+1, Dchi2clCutValue[iPt][iY]);
      Dtrk1PtCut.SetBinContent(iPt+1, iY+1, Dtrk1PtCutValue[iPt][iY]);
      Dtrk2PtCut.SetBinContent(iPt+1, iY+1, Dtrk2PtCutValue[iPt][iY]);
    }
  }
}

bool hfupc::dtree::tightsel(int j)
{
  //IF YOU CHANGE VARIABLES USED BY THIS FUNCTION,
  //PLEASE EDIT tightselVarHI_, tightselVarPP_ APPROPRIATELY!
  //Protections for cut
  if(!tightselVarsFound_) return cutProtectionMsg(__PRETTY_FUNCTION__);
  float Dpt = bvf_["Dpt"][j];
  float Dy = bvf_["Dy"][j];
  int bin = DalphaCut.FindBin(Dpt, Dy);
  // 2023.12.14: DdethetaBScorr is missing in current production
  // Replace w/ cut not using DdethetaBScorr + opening issue in repo for a better long term fix
  //  bool cut = bvf_["DlxyBS"][j]/bvf_["DlxyBSErr"][j] > 3.5 && bvf_["DdthetaBScorr"][j] < 0.2;
  bool cut = bvf_["DsvpvDistance"][j]/bvf_["DsvpvDisErr"][j] > DsvpvSigCut[bin] &&
    bvf_["Dalpha"][j] < DalphaCut[bin] &&
    bvf_["Ddtheta"][j] < DdthetaCut[bin] &&
                         bvf_["Dchi2cl"][j] > Dchi2clCut[bin]
    ;
  if(ishi_) cut = true;

  return cut;
}

/** Use varied cut values for specified variables
 */
void hfupc::dtree::varySelection(std::vector<VariedCutType> variedCut){
  for (auto& c : variedCut) {
    int cuttype = static_cast<int>(c);
    auto cutvec = cutValues[cuttype];
    std::copy(std::begin(variedCutValues[cuttype]),
              std::end(variedCutValues[cuttype]), cutvec->begin());
    fillDcut();
  }
}

#endif
