#include "header.h"

class SpectraClass{
  public:
    SpectraClass();
    SpectraClass(bool mc, TString Period = "13TeV", bool light = true, int skipping = 0);
    ~SpectraClass();

    unsigned int GetNClusInEvt(unsigned int i);
    void CreateHistos();
    void Process(unsigned int nclus = 0);
    void NextEvt();
    void CalculatePi0s(unsigned int kCorrection);
    Cluster *GetCluster(unsigned int i = 0);
    std::array<float, 3> GetMomVect(unsigned int i = 0);
    float GetE(unsigned int i)                          {T->GetEvent(i); return applyNL(fenergy, Cluster_NumCells);};
    int GetClusterPDG(unsigned int i)                          {T->GetEvent(i); if(isMC == 0){ return 0;} else { return ClusterPDG;};};
    unsigned int GetClusNCells(unsigned int i)          {T->GetEvent(i); return   Cluster_NumCells;};
    inline float EtaToTheta(float eta = 0);
    void FillNCellVsEGammas();
    void WriteToFile();
    void FillHitmap();
    bool IsIsolated(unsigned int i);
    bool IsBorderCell(unsigned int clus, unsigned int NCellDist);
    bool IsTrackMatched(unsigned int i);
    bool IsAcceptedByNCellEffi(unsigned int i, unsigned int kCorrection = 0);
    bool IsExoticClus(unsigned int i, float eThresh = 0.07);
    void FillElectronHist();
    void FillAllCLusterHist();
    float applyNL(float e, int numcells);
    bool IsClusAcceptedByThreshold(unsigned int i, float agg, float thr);

  private:

    Event Evt;
    int fdebug = 0;
    bool fDoLight = true;

    AliEMCALGeometry * geom = nullptr;

    TRandom randomEng;

    double fNCellEfficiencyParams[10] = {0};


    static constexpr int fnEMCalGapsPhi = 12;
    static constexpr int fnEMCalGapsEta = 3;
    std::array<float, fnEMCalGapsPhi> fEMCalGapsPhi = {1.4, 1.75, 2.1, 2.45, 2.8, 3.15, 3.27, 4.55, 4.9, 5.25, 5.6, 5.72};
    std::array<float, fnEMCalGapsEta> fEMCalGapsEta = {-0.65, 0, 0.65};

  Float_t FunctionNL_8TeV_MC(float energy){
    Double_t fNonLinearityParams[5];
    fNonLinearityParams[0] = 1.014;
    fNonLinearityParams[1] =-0.03329;
    fNonLinearityParams[2] =-0.3853;
    fNonLinearityParams[3] = 0.5423;
    fNonLinearityParams[4] =-0.4335;
     energy *= (fNonLinearityParams[0]*exp(-fNonLinearityParams[1]/energy))+
                ((fNonLinearityParams[2]/(fNonLinearityParams[3]*2.*TMath::Pi())*
                  exp(-(energy-fNonLinearityParams[4])*(energy-fNonLinearityParams[4])/(2.*fNonLinearityParams[3]*fNonLinearityParams[3]))));
    return energy;
  }

  Float_t FunctionNL_8TeV_data(float energy){
    Double_t fNonLinearityParams[7];
    fNonLinearityParams[0] =  0.99078;
    fNonLinearityParams[1] =  0.161499;
    fNonLinearityParams[2] =  0.655166;
    fNonLinearityParams[3] =  0.134101;
    fNonLinearityParams[4] =  163.282;
    fNonLinearityParams[5] =  23.6904;
    fNonLinearityParams[6] =  0.978;
      energy *= fNonLinearityParams[6]/(fNonLinearityParams[0]*(1./(1.+fNonLinearityParams[1]*exp(-energy/fNonLinearityParams[2]))*1./(1.+fNonLinearityParams[3]*exp((energy-fNonLinearityParams[4])/fNonLinearityParams[5]))));

    return energy;
  }

  Float_t FunctionNL_OfficialTB_100MeV_MC_V2(Float_t e){// 1.5% shift instead of 5% to make the scale correct
    Double_t funcParams[5] = {1.09357, 0.0192266, 0.291993, 370.927, 694.656};
    return ( 1.00 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
  }
  // Float_t FunctionNL_OfficialTB_100MeV_MC_V2(Float_t e){
  //   Double_t funcParams[5] = {1.09357, 0.0192266, 0.291993, 370.927, 694.656};
  //   return ( 1.00 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
  // }
  //________________________________________________________________________
  Float_t FunctionNL_kSDM(Float_t e, Float_t p0 = 1., Float_t p1 = 1., Float_t p2 = 1., Float_t p3 = 1.){
    return ( p0 + p3 * exp( p1 + ( p2 * e ) ) );
  }
  //________________________________________________________________________
  Float_t FunctionNL_OfficialTB_100MeV_Data_V2(Float_t e){ // 1.5% shift instead of 5% to make the scale correct
    Double_t funcParams[5] = {1.91897, 0.0264988, 0.965663, -187.501, 2762.51};
    // return ( 1.015 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
    return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
  }
  // Float_t FunctionNL_OfficialTB_100MeV_Data_V2(Float_t e){
  //   Double_t funcParams[5] = {1.91897, 0.0264988, 0.965663, -187.501, 2762.51};
  //   return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
  // }

  inline bool fileExists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

  //--------------------------------------
  //--- Switch between energies
  //--------------------------------------
    int period = 0; // 0 = 13TeV, 1 = 8TeV
    enum enumPeriod { kpp13TeV, kpp8TeV };

    //--------------------------------------
    //--- Tree and File
    //--------------------------------------
    TFile *dataFile = nullptr;
    TChain *T = nullptr;

    //--------------------------------------
    //--- Cluster properties from tree
    //--------------------------------------
    float fenergy;
    float eta;
    float phi;
    float *CellEnergy = new float[100];
    Int_t Cluster_NumCells;
    int Cluster_CellsID[50];
    bool isEMC = false;

    bool fDoSkipping = false;
    int SkippingFrac = 1;

    int fnumNCelVar = 1;


    // MC stuff
    int ClusterPDG;
    int MCLabel = 0;


    //--------------------------------------
    //--- Surrounding tracks from tree
    //--------------------------------------
    int Surrounding_Cells_N;
    float Surrounding_Cells_R[1000];
    float Surrounding_Cells_E[1000];
    int Surrounding_Cells_ID[1000];
    float Surrounding_Cells_relEta[1000];
    float Surrounding_Cells_relPhi[1000];
    // surrounding stuff tracks
    int Surrounding_Tracks_N;
    float Surrounding_Tracks_R[1000];
    float Surrounding_Tracks_P[1000];
    // float Surrounding_Tracks_ID[1000];
    float Surrounding_Tracks_relEta[1000];
    float Surrounding_Tracks_relPhi[1000];
    float Surrounding_Tracks_nSigdEdx[1000];
    bool  Surrounding_Tracks_IsV0[1000];

    //--------------------------------------
    //--- Keep track of the current event and the number of custers in that event
    //--------------------------------------
    unsigned int currEvt = 0;
    unsigned int nclus = 0;
    // /event stuff
    float vertexZ;


    //--------------------------------------
    //--- TH2 NCell Vs E histograms
    //--------------------------------------
    TH2F* hNCellVsEGammas = nullptr;
    TH2F* hNCellVsEGammasNL = nullptr;
    TH2F* hNCellVsEGammasNLTrue = nullptr;
    TH2F* hNCellVsEGammasNLTrueElec = nullptr;

    TH2F* hNCellVsEGammasLeft = nullptr;
    TH2F* hNCellVsEGammasNLLeft = nullptr;
    TH2F* hNCellVsEGammasNLTrueLeft = nullptr;
    TH2F* hNCellVsEGammasNLTrueElecLeft = nullptr;

    TH2F* hNCellVsEGammasWide = nullptr;
    TH2F* hNCellVsEGammasNLWide = nullptr;
    TH2F* hNCellVsEGammasNLTrueWide = nullptr;
    TH2F* hNCellVsEGammasNLTrueElecWide = nullptr;

    TH2F* hNCellVsEGammasSideBand = nullptr;
    TH2F* hNCellVsEGammasNLSideBand = nullptr;
    TH2F* hNCellVsEGammasNLTrueSideBand = nullptr;
    TH2F* hNCellVsEGammasNLTrueElecSideBand = nullptr;

    TH2F* hNCellVsEGammasOnlyHighClus = nullptr;
    TH2F* hNCellVsEGammasNLOnlyHighClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueOnlyHighClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueElecOnlyHighClus = nullptr;

    TH2F* hNCellVsEGammasOnlyLowClus = nullptr;
    TH2F* hNCellVsEGammasNLOnlyLowClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueOnlyLowClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueElecOnlyLowClus = nullptr;

    TH2F* hNCellVsEGammasSideBandOnlyHighClus = nullptr;
    TH2F* hNCellVsEGammasNLSideBandOnlyHighClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueSideBandOnlyHighClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueElecSideBandOnlyHighClus = nullptr;

    //--------------------------------------
    //--- TH2 NCell Vs E histograms for clusterizer studies
    //--------------------------------------
    TH2F* hNCellVsENL = nullptr;
    TH2F* hNCellVsETMNL = nullptr;
    TH2F* hNCellVsENoExotic50 = nullptr;
    TH2F* hNCellVsENoExotic75 = nullptr;
    TH2F* hNCellVsETMNLS500A100 = nullptr;
    TH2F* hNCellVsETMNLS500A105 = nullptr;
    TH2F* hNCellVsETMNLS500A110 = nullptr;

    //--------------------------------------
    //--- Inv. Mass histograms
    //--------------------------------------
    std::vector<TH2F*> hInvMassVsPt = {nullptr};
    std::vector<TH2F*> hInvMassVsPtBack = {nullptr};
    // TH2F* hInvMassVsPtBack = nullptr;
    // TH2F* hInvMassVsPtGG = nullptr;
    // TH2F* hInvMassVsPtGC = nullptr;
    // TH2F* hInvMassVsPtCC = nullptr;
    // TH2F* hInvMassVsHighGammaPt = nullptr;
    // TH2F* hInvMassVsHighGammaPtBack = nullptr;
    // TH2F* hInvMassVsBothGammaPtBack = nullptr;
    // TH2F* hInvMassVsLowGammaPtBack = nullptr;
    // TH2F* hInvMassVsPtHigh = nullptr;
    // TH2F* hInvMassVsPtHighGamma = nullptr;
    // TH2F* hInvMassVsPtHighElec = nullptr;
    // TH2F* hInvMassVsPtBoth = nullptr;
    // TH2F* hInvMassVsPtBothGamma = nullptr;
    // TH2F* hInvMassVsPtBothElec = nullptr;
    // TH2F* hInvMassVsPtLow = nullptr;
    // TH2F* hInvMassVsPtLowGamma = nullptr;
    // TH2F* hInvMassVsPtLowElec = nullptr;
    // TH2F* hInvMassVsPtHigh1cell = nullptr;
    // TH2F* hInvMassVsPtHigh2cell = nullptr;
    // TH2F* hInvMassVsPtHigh3cell = nullptr;
    //
    // TH2F* hInvMassVsPtDoubleCount = nullptr;

    //--------------------------------------
    //--- NCell vs. E for electron clusters
    //--------------------------------------
    TH2F* hNCellVsEelec = nullptr;
    TH2F* hNCellVsEelecNL = nullptr;
    TH2F* hNCellVsEelecNLTrue = nullptr;

    //--------------------------------------
    //--- Random stuff
    //--------------------------------------
    TH2F* hHitmap = nullptr;
    TH1F* hClusE = nullptr;
    TH1F* hClusENL = nullptr;

    //--------------------------------------
    //--- vectors for storing the photon ids selected by the mass cut
    //--------------------------------------
    std::vector<unsigned int> vecGoodGammas;
    std::vector<unsigned int> vecGoodGammasLeft;
    std::vector<unsigned int> vecGoodGammasOnlyHigh;
    std::vector<unsigned int> vecGoodGammasOnlyLow;
    std::vector<unsigned int> vecGoodGammasSideBandOnlyHigh;
    std::vector<unsigned int> vecGoodGammasWide;
    std::vector<unsigned int> vecGoodGammasSideBand;

    //--------------------------------------
    //--- For mass window estimation
    //--------------------------------------
    TGraph * grMass = nullptr;
    TF1 *fMassPos = nullptr;

    //--------------------------------------
    //--- Binning for histograms
    //--------------------------------------
    const int nBinsE = 23;
    Double_t arrEbins[24] = {0.0, 0.7, 0.8, 0.9, 1.0,   1.2, 1.4, 1.6, 1.8, 2.0,
                             2.5, 3.0, 3.5, 4.0, 5.0,   6.0, 7.0, 8.0, 10., 12.0,
                             14.0, 16.0, 20.0, 25.0};
    // const int nBinsE = 18;
    // Double_t arrEbins[19] = {0.0, 0.7, 0.8, 0.9, 1.0,   1.2, 1.4, 1.6, 1.8, 2.0,
    //                          2.5, 3.0, 3.5, 4.0, 5.0,   6.0, 7.0, 8.0, 10.};


    bool isMC;

};

//----------------------------
SpectraClass::SpectraClass(): randomEng(){

}
//----------------------------
SpectraClass::~SpectraClass(){

}
//----------------------------
SpectraClass::SpectraClass(bool mc, TString Period, bool light, int skipping): randomEng(){
  isMC = mc;
  fDoLight = light;
  SkippingFrac = skipping;
  if(skipping > 0) fDoSkipping = true;
  if(Period.Contains("13TeVNomB")) period = 0;
  else if(Period.Contains("13TeVLowB")) period = 1;
  else if(Period.Contains("8TeV")) period = 2;

  fnumNCelVar = 2;
  if(mc) fnumNCelVar = 10;

  hInvMassVsPt.clear();
  hInvMassVsPtBack.clear();

  //--------------------------------------
  //--- Put together input files
  //--------------------------------------

  //------------------------ 13TeV nomB (18m)------------------------------------//
  if(period == 0){

    TString NameTree = "ClusterQA_00010113_4117900050020000000";
    // if(fDoLight) NameTree = "2ClusterTree";
    T = new TChain(NameTree);
    // std::vector<int> runs = {285471};
    std::vector<int> runs = {290323, 290324, 290327, 290350, 290374, 290375, 290376, 290399, 290401, 290411, 290412, 290425, 290426, 290427, 290428, 290456, 290458, 290459, 290469, 290499, 290500, 290538, 290539, 290540, 290544, 290549, 290550, 290553, 290588, 290590, 290612, 290613, 290614, 290615, 290627, 290632, 290645, 290658, 290660, 290665, 290687, 290689, 290692, 290696, 290699, 290721, 290742, 290764, 290766, 290769, 290774, 290787, 290790, 290841, 290843, 290846, 290848, 290860, 290862, 290886, 290887, 290892, 290894, 290895, 290932, 290935, 290941, 290943, 290944, 290948, 290974, 290975, 290976, 290979, 290980, 291002, 291003, 291004, 291005, 291035, 291037, 291041, 291065, 291066, 291069, 291093, 291100, 291101, 291110, 291111, 291116, 291143, 291188, 291209, 291240, 291257, 291262, 291263, 291265, 291266, 291282, 291283, 291284, 291285, 291286, 291360, 291361, 291362, 291373, 291375, 291377, 291397, 291399, 291402, 291416, 291417, 291419, 291420, 291424, 291446, 291447, 291453, 291456, 291457, 291481, 291482, 291484, 291485, 291590, 291614, 291615, 291618, 291622, 291624, 291626, 291657, 291661, 291665, 291690, 291692, 291694, 291697, 291698, 291702, 291706, 291729, 291755, 291756, 291760, 291768, 291769, 291795, 291796, 291803, 291942, 291943, 291944, 291945, 291946, 291948, 291953, 291976, 291977, 291982, 292012, 292040, 292060, 292061, 292062, 292067, 292075, 292080, 292081, 292106, 292107, 292108, 292109, 292114, 292115, 292140, 292160, 292161, 292162, 292163, 292164, 292166, 292167, 292168, 292192, 292218, 292240, 292242, 292265, 292270, 292273, 292274, 292298, 292397, 292398, 292405, 292406, 292428, 292429, 292430, 292432, 292434, 292456, 292457, 292460, 292461, 292495, 292496, 292497, 292500, 292521, 292523, 292524, 292526, 292553, 292554, 292557, 292559, 292560, 292563, 292584, 292586, 292693, 292695, 292696, 292698, 292701, 292704, 292737, 292739, 292747, 292748, 292750, 292803, 292804, 292809, 292810, 292811, 292831, 292832, 292836, 292839, 292832, 292836, 292839};
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(runs), std::end(runs), rng);
    if(isMC){
      for(const int& i : runs){
        if(!fileExists(Form("/media/joshua/Elements/Data/2021_02_23_18m_Tree/MC/%i/AnalysisResults.root",i))) continue;
        T->Add(Form("/media/joshua/Elements/Data/2021_02_23_18m_Tree/MC/%i/AnalysisResults.root",i));
      }
    } else {
      for(const int& i : runs){
        if(!fileExists(Form("/media/joshua/Elements/Data/2021_02_23_18m_Tree/data/%i/AnalysisResults.root",i))) continue;
        T->Add(Form("/media/joshua/Elements/Data/2021_02_23_18m_Tree/data/%i/AnalysisResults.root", i));
      }
    }
  }
  //------------------------ 13TeV lowB (18c) ------------------------------------//
  if(period == 1){

    TString NameTree = "ClusterQA_00010113_4117900050020000000";
    // if(fDoLight) NameTree = "2ClusterTree";
    T = new TChain(NameTree);
    // std::vector<int> runs = {285471};
    std::vector<int> runs = {285471, 285481, 285486, 285496, 285497, 285545, 285550, 285557, 285575, 285576, 285577, 285578, 285599, 285601, 285602, 285603, 285639, 285640, 285641, 285642, 285643, 285662, 285663, 285664, 285666, 285697, 285698, 285722, 285753, 285754, 285755, 285777, 285778, 285781, 285804, 285810, 285811, 285812, 285830, 285851, 285869, 285893, 285917, 285946, 285957, 285958};
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(runs), std::end(runs), rng);
    if(isMC){
      for(const int& i : runs){
        if(!fDoLight) T->Add(Form("/media/joshua/external_drive/Data/2021_01_18_pp13TeVlowBTree/LHC18/LHC18h1_fast/%i/AnalysisResults.root",i));
        else T->Add(Form("/media/joshua/external_drive/Data/2021_01_18_pp13TeVlowBTree/LHC18/LHC18h1_fast/%i/SkimmedTree.root",i));
      }
    } else {
      for(const int& i : runs){
        if(!fDoLight) T->Add(Form("/media/joshua/external_drive/Data/2021_01_18_pp13TeVlowBTree/LHC18/LHC18c_child1/%i/AnalysisResults.root", i));
        else T->Add(Form("/media/joshua/external_drive/Data/2021_01_18_pp13TeVlowBTree/LHC18/LHC18c_child1/%i/SkimmedTree.root", i));
      }
    }
  }

  //------------------------ 8TeV ------------------------------------//
  else if(period == 2){
    TString NameTree = "ClusterQA_00010113_4117900050020000000";
    // if(fDoLight) NameTree = "2ClusterTree";
    T = new TChain(NameTree);
    std::vector<int> runs = {183916, 183932, 183933, 183934, 183935, 183936, 183937, 183938, 183942, 183946, 184126, 184127, 184131, 184132, 184134, 184135, 184137, 184138, 184140, 184144, 184145, 184147, 184183, 184188, 184208, 184209, 184210, 184215, 184216, 184371, 184383, 184389, 184673, 184678, 184682, 184687, 184784, 184786, 185029, 185031, 185116, 185126, 185127, 185132, 185133, 185134, 185157, 185160, 185164, 185189, 185196, 185198, 185203, 185206, 185208, 185217, 185221, 185282, 185284, 185289, 185291, 185292, 185293, 185296, 185299, 185300, 185302, 185303, 185349, 185350, 185351, 185356, 185359, 185360, 185361, 185362, 185363, 185368, 185371, 185375, 185378, 185461, 185465, 185474, 185582, 185583, 185588, 185589, 185659, 185687, 185738, 185764, 185765, 185768, 185775, 185776, 185778, 185784, 186007, 186009, 186011, 186163, 186164, 186165, 186167, 186205, 186208, 186319, 186320};
    if(!isMC) runs = { 183913, 183916, 183932, 183933, 183934, 183935, 183936, 183937, 183938, 183942, 184126, 184127, 184131, 184132, 184134, 184135, 184137, 184138, 184140, 184144, 184145, 184183, 184188, 184208, 184209, 184210, 184216, 184371, 184383, 184389, 184673, 184678, 184682, 184687, 184784, 184786, 185029, 185031, 185116, 185126, 185127, 185132, 185133, 185134, 185157, 185160, 185164, 185189, 185196, 185198, 185203, 185206, 185208, 185217, 185221, 185282, 185284, 185291, 185292, 185293, 185296, 185299, 185300, 185302, 185303, 185349, 185351, 185356, 185359, 185360, 185361, 185362, 185363, 185368, 185371, 185375, 185378, 185461, 185465, 185474, 185582, 185583, 185588, 185589, 185659, 185687, 185738, 185764, 185765, 185768, 185775, 185776, 185778, 185784, 186007, 186009, 186011, 186163, 186164, 186165, 186167, 186205, 186208, 186319, 186320, 186668, 186688, 186689, 186690, 186692, 186694, 186811, 186814, 186938, 186939, 186966, 186969, 186990, 186992, 186994, 187143, 187145, 187146, 187147, 187148, 187149, 187150, 187151, 187152, 187202, 187203, 187339, 187340, 187341, 187487, 187488, 187489, 187510, 187623, 187624, 187627, 187656, 187698, 187739, 187749, 187783, 187785, 187791, 187796, 188093, 188101, 189306, 189310, 189315, 189316, 189350, 189351, 189352, 189353, 189400, 189407, 189409, 189410, 189411, 189577, 189578, 189602, 189603, 189605, 189610, 189611, 189612, 189616, 189621, 189623, 189647, 189648, 189650, 189654, 189656, 189658, 189659, 189696, 189697, 189698, 190150, 190209, 190210, 190212, 190213, 190214, 190215, 190216, 190240, 190303, 190305, 190307, 190337, 190338, 190340, 190341, 190342, 190344, 190386, 190388, 190389, 190390, 190392, 190393, 190416, 190417, 190418, 190419, 190421, 190422, 190424, 190425, 190895, 190898, 190903, 190904, 190968, 190970, 190974, 190975, 190979, 190981, 190983, 190984, 191129, 191227, 191229, 191230, 191231, 191245, 191247, 191248, 191450, 191451, 192072, 192073, 192075, 192128, 192136, 192140, 192141, 192172, 192174, 192177, 192194, 192197, 192200, 192201, 192202, 192205, 192246, 192344, 192347, 192348, 192349, 192415, 192417, 192453, 192461, 192468, 192471, 192492, 192499, 192505, 192510, 192535, 192542, 192548, 192551, 192729, 192731, 192732};
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(runs), std::end(runs), rng);
    if(isMC){
      for(const int& i : runs){
        if(!fDoLight) T->Add(Form("/media/joshua/Elements/Data/2021_01_27_8TeVTree/MC/%i/AnalysisResults.root",i));
        else T->Add(Form("/media/joshua/Elements/Data/2021_01_27_8TeVTree/MC/%i/SkimmedTree.root",i));
      }
    } else {
      for(const int& i : runs){
        if(!fDoLight) T->Add(Form("/media/joshua/Elements/Data/2021_01_27_8TeVTree/data2/%i/AnalysisResults.root", i));
        else T->Add(Form("/media/joshua/Elements/Data/2021_01_27_8TeVTree/data2/%i/SkimmedTree.root", i));
      }
    }
  }

  //--------------------------------------
  //--- Setting up the tree
  //--------------------------------------

  // ClusterQA_00000113_4117900050020000000
  // T->Print();

  T->SetBranchAddress("Cluster_E", &fenergy);
  cerr<<__LINE__<<endl;
  T->SetBranchAddress("Cluster_IsEMCAL", &isEMC);
  T->SetBranchAddress("Cluster_Eta", &eta);
  T->SetBranchAddress("Cluster_Phi", &phi);
  T->SetBranchAddress("Cluster_Cells_E", CellEnergy);
  T->SetBranchAddress("Cluster_NumCells", &Cluster_NumCells);
  T->SetBranchAddress("Cluster_Cells_ID", Cluster_CellsID);
  T->SetBranchAddress("Event_Vertex_X", &vertexZ);

  if(isMC)T->SetBranchAddress("Cluster_MC_FirstLabel", &MCLabel);



  T->SetBranchAddress("Surrounding_NCells", &Surrounding_Cells_N);
  T->SetBranchAddress("Surrounding_Cells_R", Surrounding_Cells_R);
  T->SetBranchAddress("Surrounding_Cells_E", Surrounding_Cells_E);
  T->SetBranchAddress("Surrounding_Cells_ID", Surrounding_Cells_ID);
  T->SetBranchAddress("Surrounding_Cells_RelativeEta", Surrounding_Cells_relEta);
  T->SetBranchAddress("Surrounding_Cells_RelativePhi", Surrounding_Cells_relPhi);

  T->SetBranchAddress("Surrounding_NTracks", &Surrounding_Tracks_N);
  T->SetBranchAddress("Surrounding_Tracks_R", Surrounding_Tracks_R);
  T->SetBranchAddress("Surrounding_Tracks_P", Surrounding_Tracks_P);
  T->SetBranchAddress("Surrounding_Tracks_RelativeEta", Surrounding_Tracks_relEta);
  T->SetBranchAddress("Surrounding_Tracks_RelativePhi", Surrounding_Tracks_relPhi);
  T->SetBranchAddress("Surrounding_Tracks_nSigdEdxE", Surrounding_Tracks_nSigdEdx);
  T->SetBranchAddress("Surrounding_Tracks_V0Flag", Surrounding_Tracks_IsV0);
  if(isMC){
    T->SetBranchAddress("Mother_MC_Label", &ClusterPDG);
  }
}

//--------------------------------------
//--- Setting up the histos
//--------------------------------------
void SpectraClass::CreateHistos(){

  hHitmap = new TH2F("hHitmap", "", 160, -0.8, 0.8, 126*2, 0, 6.2);

}

Cluster *SpectraClass::GetCluster(unsigned int i){
  T->GetEvent(i);
  // for(int k = 0; k < Surrounding_Tracks_N; k++){
  //   cout<<"Surrounding_Cells_ID["<<k<<"] = "<<Surrounding_Cells_ID[k]<<endl;
  // }
  Cluster *tmp;
  if(!isMC){
      tmp = new Cluster(applyNL(fenergy, Cluster_NumCells), eta, phi, Cluster_NumCells , CellEnergy, Cluster_CellsID,
      Surrounding_Cells_N, Surrounding_Cells_E, Surrounding_Cells_R, Surrounding_Cells_ID,
      Surrounding_Tracks_N, Surrounding_Tracks_P, Surrounding_Tracks_R, Surrounding_Tracks_nSigdEdx, Surrounding_Tracks_IsV0,
      false, 0, 0);
  } else {
      tmp = new Cluster(applyNL(fenergy, Cluster_NumCells), eta, phi, Cluster_NumCells , CellEnergy, Cluster_CellsID,
      Surrounding_Cells_N, Surrounding_Cells_E, Surrounding_Cells_R, Surrounding_Cells_ID,
      Surrounding_Tracks_N, Surrounding_Tracks_P, Surrounding_Tracks_R, Surrounding_Tracks_nSigdEdx, Surrounding_Tracks_IsV0,
      isMC, ClusterPDG, MCLabel);
  }
  return tmp;
}

//--------------------------------------
//--- get the number of clusters in current event, use he vertex position as an indicator
//--------------------------------------
unsigned int SpectraClass::GetNClusInEvt(unsigned int i){
  T->GetEvent(i);
  float savedVtx = vertexZ;
  unsigned int nClus = 0;
  while(savedVtx == vertexZ){
    nClus++;
    T->GetEvent(i+nClus);
  }
  T->GetEvent(i);
  return nClus;
}

//--------------------------------------
//--- Switch to the next event
//--------------------------------------
void SpectraClass::NextEvt(){
  nclus = GetNClusInEvt(currEvt);
  if(fdebug) cout<<"cluster in evt: "<<nclus<<endl;
  T->GetEvent(currEvt + nclus);
  currEvt += nclus;
  nclus = GetNClusInEvt(currEvt);
  if(fdebug == 2)cout<<__LINE__<<endl;
  Evt.Reset();
  if(fdebug == 2)cout<<__LINE__<<endl;
  for(unsigned int i = 0; i < nclus; ++i){
    Evt.InsertCluster(GetCluster(currEvt + i), i);
  }
}

//--------------------------------------
//--- Main function
//--------------------------------------
void SpectraClass::Process(unsigned int maxNumClus){
  cerr<<"max events: "<<T->GetEntries()<<endl;

  CreateHistos();
  // unsigned int counter = 0;
  if(fdebug == 2) cerr<<__LINE__<<endl;
  while(currEvt < maxNumClus && currEvt < T->GetEntries()*0.99){
    // counter++;
    if(currEvt%10000 == 0) cerr<<currEvt<<" events processed..."<<endl;
    if(fdebug == 2) cerr<<__LINE__<<endl;
    NextEvt();
    if(fdebug)cout<<"energy: "<<fenergy<<endl;
    if(fdebug == 2) cerr<<__LINE__<<endl;

    for(int i = 0; i < fnumNCelVar; ++i){
      CalculatePi0s(i);
    }

  }

}

//--------------------------------------
//--- Eta to theta
//--------------------------------------
inline float SpectraClass::EtaToTheta(float Eta){
 return (float) 2.*atan(exp(-Eta));
}

//--------------------------------------
//--- Get Momntum vector from cluster
//--------------------------------------
std::array<float, 3> SpectraClass::GetMomVect(unsigned int i){
  // T->GetEvent(i);
  if(fdebug == 2)cout<<__LINE__<<endl;
  float theta = EtaToTheta(Evt.clus[i]->eta);
  float px = Evt.clus[i]->e * sin(theta) * cos(Evt.clus[i]->phi);
  float py = Evt.clus[i]->e * sin(theta) * sin(Evt.clus[i]->phi);
  float pz = Evt.clus[i]->e * cos(theta);
  std::array<float, 3> P = {px, py, pz};
  return P;
}

//--------------------------------------
//--- Select photons if the lie in the pi0 mass window together with another cluster
//--------------------------------------
void SpectraClass::CalculatePi0s(unsigned int kCorrection){
  if(fdebug)cout<<"kCorrection: "<<kCorrection<<"  hInvMassVsPt.size(): "<<hInvMassVsPt.size()<<endl;
  if(kCorrection >= hInvMassVsPt.size()){
    hInvMassVsPt.push_back( new TH2F(Form("hInvMassVsPt_%i", kCorrection), Form("hInvMassVsPt_%i", kCorrection), 500, 0, 1., 200, 0, 20));
    hInvMassVsPtBack.push_back( new TH2F(Form("hInvMassVsPtBack_%i", kCorrection), Form("hInvMassVsPtBack_%i", kCorrection), 500, 0, 1., 200, 0, 20));
  }

  std::vector<unsigned int> fDoubleCount;
  // vecGoodGammasSideBandOnlyHigh.resize(0);
  if(nclus < 2) return;
   if(fdebug) cerr<<"startEvt:"<<endl;
  for(unsigned int ig1 = 0; ig1 < nclus; ++ig1){
    // get momentum and energy
    if(fdebug) cerr<<__LINE__<<endl;
    // if(!IsIsolated(ig1)) continue;
    if(fdebug) cerr<<__LINE__<<endl;
    if(IsTrackMatched(ig1)) continue;

    if(!IsAcceptedByNCellEffi(ig1, kCorrection)) continue;
    if(fdebug) cerr<<__LINE__<<endl;
    // if(IsBorderCell(ig1, 3)) continue;
    if(fdebug) cerr<<__LINE__<<endl;
    int PDGig1 = Evt.clus[ig1]->ClusterPDG;
    // int PDGig1 = GetClusterPDG(ig1);
    std::array<float, 3> Pg1 = GetMomVect(ig1);
    float Eg1 = Evt.clus[ig1]->e;

    // write them into TLorentzvector
    TLorentzVector LVg1;
    LVg1.SetX(Pg1[0]);//, Pg1[1], Pg1[2], Eg1);
    LVg1.SetY(Pg1[1]);//, Pg1[1], Pg1[2], Eg1);
    LVg1.SetZ(Pg1[2]);//, Pg1[1], Pg1[2], Eg1);
    LVg1.SetE(Eg1);//, Pg1[1], Pg1[2], Eg1);
    // cerr<<"Mass Photon: "<<LVg1.M()<<endl;
    if(fdebug) cerr<<__LINE__<<endl;
    for(unsigned int ig2 = ig1 + 1; ig2 < nclus; ++ig2){
      // if(!IsIsolated(ig2)) continue;
      if(IsTrackMatched(ig2)) continue;

      if(!IsAcceptedByNCellEffi(ig2, kCorrection)) continue;
      // if(IsBorderCell(ig2, 3)) continue;

      int PDGig2 = Evt.clus[ig2]->ClusterPDG;
      // int PDGig2 = GetClusterPDG(ig2);
      std::array<float, 3> Pg2 = GetMomVect(ig2);
      float Eg2 = Evt.clus[ig2]->e;
      TLorentzVector LVg2;
      LVg2.SetX(Pg2[0]);//, Pg1[1], Pg1[2], Eg1);
      LVg2.SetY(Pg2[1]);//, Pg1[1], Pg1[2], Eg1);
      LVg2.SetZ(Pg2[2]);//, Pg1[1], Pg1[2], Eg1);
      LVg2.SetE(Evt.clus[ig2]->e);
      // combine to mother particle
      TLorentzVector Pi0 = (LVg1 + LVg2);
      // cerr<<Pi0.M()<<endl;

      hInvMassVsPt.at(kCorrection)->Fill(Pi0.M(), Pi0.Pt());

      //-----------------------------------------------
      // rotation background
      //-----------------------------------------------
      if(fdebug) cerr<<__LINE__<<endl;
      // rotate around mother axis
      TVector3 pi0Vec = Pi0.Vect();
      LVg1.Rotate(TMath::Pi()/2., pi0Vec);
      LVg2.Rotate(TMath::Pi()/2., pi0Vec);
      for(unsigned int ig3 = 0; ig3 < nclus; ++ig3){
        if(ig3 == ig2 || ig3 == ig1) continue;
        // if(!IsIsolated(ig3)) continue;
        if(IsTrackMatched(ig3)) continue;

        if(!IsAcceptedByNCellEffi(ig3, kCorrection)) continue;
        // int PDGig3 = GetClusterPDG(currEvt+ig3);
        std::array<float, 3> Pg3 = GetMomVect(ig3);
        // float Eg3 = GetE(ig3);


        // write them into TLorentzvector
        TLorentzVector LVg3;
        LVg3.SetX(Pg3[0]);//, Pg1[1], Pg1[2], Eg1);
        LVg3.SetY(Pg3[1]);//, Pg1[1], Pg1[2], Eg1);
        LVg3.SetZ(Pg3[2]);//, Pg1[1], Pg1[2], Eg1);
        LVg3.SetE(Evt.clus[ig3]->e);

        // combine to mother particle
        TLorentzVector BackPi01 = (LVg1 + LVg3);
        TLorentzVector BackPi02 = (LVg2 + LVg3);
        // cerr<<Pi0.M()<<endl;

        hInvMassVsPtBack.at(kCorrection)->Fill(BackPi01.M(), BackPi01.Pt());
        hInvMassVsPtBack.at(kCorrection)->Fill(BackPi02.M(), BackPi02.Pt());

      }
    }
  }
}

//--------------------------------------
//--- Hitmap
//--------------------------------------
void SpectraClass::FillHitmap(){
  for(unsigned int i = 0; i < nclus; ++i){
    if(fdebug == 2) cerr<<__LINE__<<endl;
    // T->GetEvent(currEvt + i);
    if(IsBorderCell(i, 3)) continue;
    float theta = EtaToTheta(Evt.clus[i]->eta);
    if(fdebug == 2) cerr<<__LINE__<<endl;
    // cerr<<"theta: "<<theta - TMath::Pi()*0.5<<endl;
    hHitmap->Fill(theta - TMath::Pi()*0.5, Evt.clus[i]->phi);
    if(fdebug == 2) cerr<<__LINE__<<endl;
    hClusE->Fill(Evt.clus[i]->e);
    if(fdebug == 2) cerr<<__LINE__<<endl;
    hClusENL->Fill(Evt.clus[i]->e);
    if(fdebug == 2) cerr<<__LINE__<<endl;
  }
}


//--------------------------------------
// return true if 1 cell cluster is accepted, return true for all 2+ cell clusters
//--------------------------------------
bool SpectraClass::IsAcceptedByNCellEffi(unsigned int clus, unsigned int kCorrection){
  if(Evt.clus[clus]->ncells >= 2) return true;

  int clusPDG = Evt.clus[clus]->ClusterPDG;
  float clusE = Evt.clus[clus]->e;

  float randNr = randomEng.Rndm();

  switch(kCorrection){
    case 0: //kNoCorrection:
    {
      return false;
    }
    case 1: //kAcceptAll:
    {
      return true;
    }
    case 2: //kNCeAllClusters:
    {
      // based on all clusters in data and MC
      // data clusters influenced by exotics above 2 GeV
      fNCellEfficiencyParams[0] = 2.71596e-01;
      fNCellEfficiencyParams[1] = 1.80393;
      fNCellEfficiencyParams[2] = 6.50026e-01;
      // printf("kNCeAllClusters\n");
      Float_t val = fNCellEfficiencyParams[0]*exp(
                    -0.5*((clusE-fNCellEfficiencyParams[1])/fNCellEfficiencyParams[2])*
                    ((clusE-fNCellEfficiencyParams[1])/fNCellEfficiencyParams[2]));
      if(randNr < val) return kTRUE;
      else return kFALSE;
      break;
    }

    case 3: //kNCeTestBeam:
    {
      // based on test beam measurements
      // should behave like pure photon clusters
      fNCellEfficiencyParams[0] = -4.23138e-02;
      fNCellEfficiencyParams[1] = 1.95466e-01;
      Float_t val = fNCellEfficiencyParams[0] + fNCellEfficiencyParams[1]*clusE;
      // printf("kNCeTestBeam\n");
      if(randNr < val) return kTRUE;
      else return kFALSE;
      break;
    }

    case 4: //kNCeGammaAndElec applied to all clusters. Based on PCMEDC results: file:///home/joshua/PCG_Software/EMCal_NCellEffi/PlotsCalibTrain/Corr_PCMEMC_wFit.pdf
    {
      // NonLin for obtained correction: TB with scale + FT
      // based on clusters which are part of a cluster pair with a mass of: [M(Pi0) - 0.05;M(Pi0) + 0.02]
      // mostly photon clusters and electron(conversion) clusters (purity about 95%)
      // exotics should be cancled by that
      fNCellEfficiencyParams[0] = -0.0377925;
      fNCellEfficiencyParams[1] = 0.160758;
      fNCellEfficiencyParams[2] = -0.00357992;
      Float_t val = fNCellEfficiencyParams[0]*clusE*clusE +
                    fNCellEfficiencyParams[1]*clusE +
                    fNCellEfficiencyParams[2];

      // printf("kNCeGammaAndElec\n");
      if(randNr < val) return kTRUE;
      else return kFALSE;
      break;
    }

    case 5: //kNCeGammaAndElec applied to all clusters. Based on PCMEDC results: file:///home/joshua/PCG_Software/EMCal_NCellEffi/PlotsCalibTrain/Corr_PCMEMC_wFit.pdf
    {
      // NonLin for obtained correction: TB with scale
      // based on clusters which are part of a cluster pair with a mass of: [M(Pi0) - 0.05;M(Pi0) + 0.02]
      // mostly photon clusters and electron(conversion) clusters (purity about 95%)
      // exotics should be cancled by that
      fNCellEfficiencyParams[0] = -0.0387877;
      fNCellEfficiencyParams[1] = 0.104607;
      fNCellEfficiencyParams[2] = 0.0793534;
      Float_t val = fNCellEfficiencyParams[0]*clusE*clusE +
                    fNCellEfficiencyParams[1]*clusE +
                    fNCellEfficiencyParams[2];

      // printf("kNCeGammaAndElec\n");
      if(randNr < val) return kTRUE;
      else return kFALSE;
      break;
    }

    case 6:
      //kNCeGammaAndElec applied to only gamma clusters. Based on PCMEDC results: file:///home/joshua/PCG_Software/EMCal_NCellEffi/PlotsCalibTrain/Corr_PCMEMC_wFit.pdf
      {
        if(fabs(clusPDG) != 22) return false;
        // NonLin for obtained correction: TB with scale + FT
        // based on clusters which are part of a cluster pair with a mass of: [M(Pi0) - 0.05;M(Pi0) + 0.02]
        // mostly photon clusters and electron(conversion) clusters (purity about 95%)
        // exotics should be cancled by that
        fNCellEfficiencyParams[0] = -0.0377925;
        fNCellEfficiencyParams[1] = 0.160758;
        fNCellEfficiencyParams[2] = -0.00357992;
        Float_t val = fNCellEfficiencyParams[0]*clusE*clusE +
                      fNCellEfficiencyParams[1]*clusE +
                      fNCellEfficiencyParams[2];

        // printf("kNCeGammaAndElec\n");
        if(randNr < val) return kTRUE;
        else return kFALSE;
        break;
      }

    case 7: //kNCeGammaAndElec applied to all clusters. Based on EDC results: file:///home/joshua/PCG_Software/EMCal_NCellEffi/PlotsCalibTrain/Corr_PCMEMC_wFit.pdf
    {
      // NonLin for obtained correction: TB with scale + FT
      // based on clusters which are part of a cluster pair with a mass of: [M(Pi0) - 0.05;M(Pi0) + 0.02]
      // mostly photon clusters and electron(conversion) clusters (purity about 95%)
      // exotics should be cancled by that
      fNCellEfficiencyParams[0] = 0.0110515;
      fNCellEfficiencyParams[1] = -0.035678;
      fNCellEfficiencyParams[2] = 0.077142;
      Float_t val = fNCellEfficiencyParams[0]*clusE*clusE +
                    fNCellEfficiencyParams[1]*clusE +
                    fNCellEfficiencyParams[2];

      // printf("kNCeGammaAndElec\n");
      if(randNr < val) return kTRUE;
      else return kFALSE;
      break;
    }
    case 8: //kNCeGammaAndElec applied to all clusters. Based on EDC results: file:///home/joshua/PCG_Software/EMCal_NCellEffi/PlotsCalibTrain/Corr_PCMEMC_wFit.pdf
    {
      // NonLin for obtained correction: TB with scale + FT
      // based on clusters which are part of a cluster pair with a mass of: [M(Pi0) - 0.05;M(Pi0) + 0.02]
      // mostly photon clusters and electron(conversion) clusters (purity about 95%)
      // exotics should be cancled by that
      fNCellEfficiencyParams[0] = 0.0124651;
      fNCellEfficiencyParams[1] = -0.0961889;
      fNCellEfficiencyParams[2] = 0.16844;
      Float_t val = fNCellEfficiencyParams[0]*clusE*clusE +
                    fNCellEfficiencyParams[1]*clusE +
                    fNCellEfficiencyParams[2];

      // printf("kNCeGammaAndElec\n");
      if(randNr < val) return kTRUE;
      else return kFALSE;
      break;
    }

    case 9: //kNCeGammaAndElec applied to only gamma clusters. Based on EDC results: file:///home/joshua/PCG_Software/EMCal_NCellEffi/PlotsCalibTrain/Corr_PCMEMC_wFit.pdf
    {
      if(fabs(clusPDG) != 22) return false;
      // NonLin for obtained correction: TB with scale + FT
      // based on clusters which are part of a cluster pair with a mass of: [M(Pi0) - 0.05;M(Pi0) + 0.02]
      // mostly photon clusters and electron(conversion) clusters (purity about 95%)
      // exotics should be cancled by that
      fNCellEfficiencyParams[0] = 0.0110515;
      fNCellEfficiencyParams[1] = -0.035678;
      fNCellEfficiencyParams[2] = 0.077142;
      Float_t val = fNCellEfficiencyParams[0]*clusE*clusE +
                    fNCellEfficiencyParams[1]*clusE +
                    fNCellEfficiencyParams[2];

      // printf("kNCeGammaAndElec\n");
      if(randNr < val) return kTRUE;
      else return kFALSE;
      break;
    }
  }
  return false;
}
//--------------------------------------
// return true if cell is isolated, false if its not isolated
//--------------------------------------
bool SpectraClass::IsIsolated(unsigned int clus){
  // T->GetEvent(clus);
  if(fdebug) cout<<"Evt.clus[clus]->NsurrCells: "<<Evt.clus[clus]->NsurrCells<<endl;
  for(int i = 0; i <Evt.clus[clus]->NsurrCells; ++i){
    bool isSameCell = false;
    for(int ii = 0; ii < Evt.clus[clus]->ncells; ++ii){
      if(Evt.clus[clus]->CellsID[ii] == Evt.clus[clus]->surrCellsID[i]) {
        isSameCell = true;
        break;
      }
    }
    if(isSameCell) continue;
    if(fdebug) cout<<"Evt.clus[clus]->surrCellsE[i]: "<<Evt.clus[clus]->surrCellsE[i]<<endl;
    if(Evt.clus[clus]->surrCellsE[i] > 0.1){
      if(fdebug) cout<<"Evt.clus[clus]->surrCellsR[i]: "<<Evt.clus[clus]->surrCellsR[i]<<endl;
      if(Evt.clus[clus]->surrCellsR[i] < 0.01){
        return false;
      }
    }
  }
  return true;
}

//--------------------------------------
// return true if cell is less than NCells away from SM border
//--------------------------------------
bool SpectraClass::IsBorderCell(unsigned int clus, unsigned int NCellDist){
  // T->GetEvent(clus);
  float etaphiDist = NCellDist * 0.013;
  for(int ie = 0; ie < fnEMCalGapsEta; ++ie){
    if(fEMCalGapsEta.at(ie) - etaphiDist < Evt.clus[clus]->eta && fEMCalGapsEta.at(ie) + etaphiDist > Evt.clus[clus]->eta) return true;
  }
  for(int ip = 0; ip < fnEMCalGapsPhi; ++ip){
    if(fEMCalGapsPhi.at(ip) - etaphiDist*1.5 < Evt.clus[clus]->phi && fEMCalGapsPhi.at(ip) + etaphiDist*1.5 > Evt.clus[clus]->phi) return true;
  }
  return false;
}

//-----------------------------------------------------
// return true if cluster is matched with track
//-----------------------------------------------------
bool SpectraClass::IsTrackMatched(unsigned int clus){
  // T->GetEvent(clus);
  for(int i = 0; i <Evt.clus[clus]->NsurrTracks; ++i){
    // if(Surrounding_Tracks_P[i] > 0.1){
      if(Evt.clus[clus]->surrTracksR[i] < 0.05){
        return true;
      // }
    }
  }
  return false;
}

//-----------------------------------------------------
// Fill pure electron hists (V0 tracks)
//-----------------------------------------------------
void SpectraClass::FillElectronHist(){

  for(unsigned int ic = 0; ic < Evt.NClus; ++ic){
    // T->GetEvent(currEvt + ic);
    for(int i = 0; i <Evt.clus[ic]->NsurrTracks; ++i){
      if(Evt.clus[ic]->surrTracksR[i] < 0.02){
        if(Evt.clus[ic]->surrTracksIsV0[i] == true){
          // if(fabs(Surrounding_Tracks_nSigdEdx[i]) < 2){
            // hNCellVsEelec->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
            hNCellVsEelecNL->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
            if(isMC){
              if(fabs(Evt.clus[ic]->ClusterPDG) == 11){
                hNCellVsEelecNLTrue->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
              }
            }
          // }
        }
      }
    }
  }
}

//-----------------------------------------------------
// Check if the cluster has a cell with at least eThresh next to it
// If not its considered an exotic cluster
//-----------------------------------------------------
bool SpectraClass::IsExoticClus(unsigned int ic, float eThresh){

  if(Evt.clus[ic]->e < 3) return false;

  if(!geom) geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");

  if(!geom)
  {
    cerr<<"DEBUG: NO geometry found..."<<endl;
    return false;
  }

  // no exotic if its above 1 cell
  float eCross = 0;
  int largestCellIndex = -1;
  float largestE = 0;
  for(int i = 0; i < Evt.clus[ic]->ncells; ++i){
    if(Evt.clus[ic]->CellEnergy[i] > largestE) {
      largestE = Evt.clus[ic]->CellEnergy[i];
      largestCellIndex = Evt.clus[ic]->CellsID[i];
    }
  }
  int absID = largestCellIndex;

  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
  geom->GetCellIndex(absID,imod,iTower,iIphi,iIeta);
  geom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);

  // Get close cells index, energy and time, not in corners

  Int_t absID1 = -1;
  Int_t absID2 = -1;

  if ( iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = geom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if ( iphi > 0 )                                absID2 = geom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);

  // cout<<"absID1: "<<absID1<<endl;
  for(int i = 0; i < Evt.clus[ic]->NsurrCells; ++i){
    // cout<<"Evt.clus[ic]->surrCellsID[i]: "<<Evt.clus[ic]->surrCellsID[i]<<endl;
    if(Evt.clus[ic]->surrCellsID[i] == absID1){
      eCross += Evt.clus[ic]->surrCellsE[i];
      // if(Evt.clus[ic]->surrCellsE[i] < eThresh){
      //   return true;
      // }
    }
    if(Evt.clus[ic]->surrCellsID[i] == absID2){
      eCross += Evt.clus[ic]->surrCellsE[i];
      // if(Evt.clus[ic]->surrCellsE[i] < eThresh){
      //   return true;
      // }
    }
  }


  // In case of cell in eta = 0 border, depending on SM shift the cross cell index

  Int_t absID3 = -1;
  Int_t absID4 = -1;

  if ( ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2) )
  {
    absID3 = geom-> GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = geom-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
  }
  else if ( ieta == 0 && imod%2 )
  {
    absID3 = geom-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = geom-> GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
  }
  else
  {
    if ( ieta < AliEMCALGeoParams::fgkEMCALCols-1 )
      absID3 = geom-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if ( ieta > 0 )
      absID4 = geom-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
  }


  for(int i = 0; i < Evt.clus[ic]->NsurrCells; ++i){
    if(Evt.clus[ic]->surrCellsID[i] == absID3){
      eCross += Evt.clus[ic]->surrCellsE[i];
      // if(Evt.clus[ic]->surrCellsE[i] < eThresh){
      //   return true;
      // }
    }
    if(Evt.clus[ic]->surrCellsID[i] == absID4){
      eCross += Evt.clus[ic]->surrCellsE[i];
      // if(Evt.clus[ic]->surrCellsE[i] < eThresh){
      //   return true;
      // }
    }
  }
  if(Evt.clus[ic]->ncells == 1){
    // if(eCross < 0.1)cout<<"eCross: "<<eCross<<"  eThresh: "<<eThresh<<"  Evt.clus[ic]->e: "<< Evt.clus[ic]->e<<endl;
    if(eCross < eThresh) return true;
    else return false;
  }
  else if(eCross < 0.1* Evt.clus[ic]->e){
    return true;
  }

  // all 4 neighbours were checked, cell is not next to cluster
  return false;
}


//-----------------------------------------------------
// FIll all clusters
//-----------------------------------------------------
void SpectraClass::FillAllCLusterHist(){

  for(unsigned int ic = 0; ic < Evt.NClus; ++ic){
    // T->GetEvent(currEvt + ic);
    float Eclus = Evt.clus[ic]->e;
    if(fdebug) cout<<__LINE__<<endl;
    hNCellVsENL->Fill(Eclus, Evt.clus[ic]->ncells);
    if(fdebug) cout<<__LINE__<<endl;
    if(!IsTrackMatched(ic)){
      if(fdebug) cout<<__LINE__<<endl;
      hNCellVsETMNL->Fill(Eclus, Evt.clus[ic]->ncells);
      if(fdebug) cout<<__LINE__<<endl;
      if(IsClusAcceptedByThreshold(ic, 0.5, 0.1)) hNCellVsETMNLS500A100->Fill(Eclus, Evt.clus[ic]->ncells);
      if(fdebug) cout<<__LINE__<<endl;
      if(IsClusAcceptedByThreshold(ic, 0.5, 0.105)) hNCellVsETMNLS500A105->Fill(Eclus, Evt.clus[ic]->ncells);
      if(IsClusAcceptedByThreshold(ic, 0.5, 0.11)) hNCellVsETMNLS500A110->Fill(Eclus, Evt.clus[ic]->ncells);

      if(!IsExoticClus(ic, 0.05)){
        hNCellVsENoExotic50->Fill(Eclus, Evt.clus[ic]->ncells);
      }
      if(!IsExoticClus(ic, 0.075)){
        hNCellVsENoExotic75->Fill(Eclus, Evt.clus[ic]->ncells);
      }

    }
  }
}

//-----------------------------------------------------
// apply NonLinearity
//-----------------------------------------------------
float SpectraClass::applyNL(float e, int numcells){

  //--------------- pp 13TeV -------------------------
  // if (period == 0){

    if(isMC){
      e /= FunctionNL_OfficialTB_100MeV_MC_V2(e);
      // fine tuning based on gaussian fits on PCMEMC in pPb5TeV
      // e /= FunctionNL_kSDM(e, 0.987912, -2.94105, -0.273207) ;  // nom. B
      if(period == 1){
        e /= FunctionNL_kSDM(e, 1.02451, -3.49297, -0.420027); // lowB
        e /= FunctionNL_kSDM(e, 0.986634, -4.12191, -0.321714);
        // if(period == 1) e *= 1.035;
      } else if(period == 0){
        e /= FunctionNL_kSDM(e, 0.987912, -2.94105, -0.273207);
      } else if(period == 2){
        e /= FunctionNL_kSDM(e, 0.987912, -2.94105, -0.273207);
        e /= 0.98;
      }
      if(numcells == 1){ // additional fine tuning for 1 cell clusters
        e /= FunctionNL_kSDM(e,0, -0.002069903, -0.00669839);
        e /= 0.995;
      }
    } else {
      e /= FunctionNL_OfficialTB_100MeV_Data_V2(e);
      if(period == 2){
        e /= 16.3/16.;
      }
    }

    return e;


  //--------------- pp 8TeV -------------------------
  // } else if (period == 1){
  //   if(isMC){
  //     e = FunctionNL_8TeV_MC(e);
  //     e *= 1.04;
  //     if(numcells == 1){ // additional fine tuning for 1 cell clusters
  //       e /= FunctionNL_kSDM(e,0, -0.002069903, -0.00669839);
  //       e /= 0.995;
  //     }
  //   } else {
  //     e = FunctionNL_8TeV_data(e);
  //   }
  //   return e;
  //
  // }
}

//--------------------------------------
//--- check if cluster passes aggregation and seed cuts
//--------------------------------------
bool SpectraClass::IsClusAcceptedByThreshold(unsigned int i, float agg, float thr){
  // T->GetEvent(i);
  bool hasAgg = false;
  if(fdebug) cout<<__LINE__<<endl;
  for(int ic = 0; ic < Evt.clus[i]->ncells; ++ic){
    if(fdebug) cout<<__LINE__<<endl;
    if(Evt.clus[i]->CellEnergy[ic] < thr) return false;
    if(fdebug) cout<<__LINE__<<endl;
    if(Evt.clus[i]->CellEnergy[ic] >= agg) hasAgg = true;
  }
  if(fdebug) cout<<__LINE__<<endl;
  return hasAgg;
}


//--------------------------------------
//--- Write everything to file
//--------------------------------------
void SpectraClass::WriteToFile(){
  gSystem->Exec("mkdir rootFiles_InvMass");

  TString filename = "rootFiles_InvMass/InvMassCorrected_13TeV_nom_04_28_MConly.root";
  if(period ==1) filename = "rootFiles_InvMass/InvMassCorrected_13TeV_low_03_25_TBNL.root";
  else if(period ==2) filename = "rootFiles_InvMass/InvMassCorrected_8TeV_03_08_ShaperTBNL.root";
  std::ifstream ftest(filename);

  TFile fout(filename, (ftest.good() == 0) ? "Recreate" : "Update");
  fout.cd();
  cout<<hInvMassVsPt.size()<<endl;
  for(unsigned int i = 0; i < hInvMassVsPt.size(); ++i){
    hInvMassVsPt.at(i)->Write(Form("hInvMassVsPt_%i_%s", i, (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
    hInvMassVsPtBack.at(i)->Write(Form("hInvMassVsPtBack_%i_%s", i,  (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  }
  fout.Close();


}



//--------------------------------------
//--- Start the thing :)
//--------------------------------------
void MakeInvMassSpectraCorected(bool isMC = 1, TString period = "13TeVNomB", bool light = false, int skip = 0){
  TStopwatch watch;
  watch.Start();
  unsigned int nEvt = 1.2e7;
  if(isMC){
    SpectraClass myAna(true, period, light, skip);
    myAna.Process(nEvt*2.5);
    myAna.WriteToFile();
  } else {
    SpectraClass myAna2(false, period, light, skip);
    myAna2.Process(nEvt);
    myAna2.WriteToFile();
  }
  watch.Stop();
  watch.Print();

}
