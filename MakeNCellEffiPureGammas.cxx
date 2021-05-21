#include "header.h"


class DataTree{
  public:
    DataTree();
    DataTree(bool mc, TString Period = "13TeV", bool light = true, int skipping = 0);
    ~DataTree();

    unsigned int GetNClusInEvt(unsigned int i);
    void CreateHistos();
    void Process(unsigned int nclus = 0);
    void NextEvt();
    void GetGammasViaPi0();
    Cluster *GetCluster(unsigned int i = 0);
    std::array<float, 3> GetMomVect(unsigned int i = 0);
    float GetE(unsigned int i)                          {T->GetEvent(i); return applyNL(energy, Cluster_NumCells);};
    int GetClusterPDG(unsigned int i)                          {T->GetEvent(i); if(isMC == 0){ return 0;} else { return ClusterPDG;};};
    unsigned int GetClusNCells(unsigned int i)          {T->GetEvent(i); return   Cluster_NumCells;};
    void SetDoNoTRD(bool tmp)                           {fDoNoTRD = tmp;};
    void SetRejectBehindTRDSupport(bool tmp)                           {fDoRejectTRDSupport = tmp;};
    inline float EtaToTheta(float eta = 0);
    void FillNCellVsEGammas();
    void WriteToFile();
    void FillHitmap();
    bool IsIsolated(unsigned int i);
    bool IsBorderCell(unsigned int clus, unsigned int NCellDist);
    void RejectBorderCells(bool tmp)                  {fRejectBorderCells = tmp;};
    bool IsTrackMatched(unsigned int i);
    bool IsExoticClus(unsigned int i, float eThresh = 0.07);
    void FillElectronHist();
    void FillAllCLusterHist();
    float applyNL(float e, int numcells);
    bool IsClusAcceptedByThreshold(unsigned int i, float agg, float thr);
    void FillMassHists(TString name);
    bool SMwithNoTRD(unsigned int clus);
    bool IsBehindTRDSupport(unsigned int clus);

  private:

    Event Evt;
    int fdebug = 0;
    bool fDoLight = true;
    bool fDoNoTRD = false;
    bool fDoRejectTRDSupport = false;
    bool fRejectBorderCells = false;

    AliEMCALGeometry * geom = nullptr;

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

  Float_t FunctionNL_OfficialTB_100MeV_MC_V2(Float_t e){
    Double_t funcParams[5] = {1.09357, 0.0192266, 0.291993, 370.927, 694.656};
    return ( /*0.965*/  (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
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
    return ( 1.015 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
    // return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
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
    float energy;
    float eta;
    float phi;
    float *CellEnergy = new float[100];
    Int_t Cluster_NumCells;
    int Cluster_CellsID[50];
    bool isEMC = false;

    bool fDoSkipping = false;
    int SkippingFrac = 1;


    // MC stuff
    int ClusterPDG;
    int MCLabel = 0;
    float MCTrueEnergy = 0;


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

    TH2F* hNCellVsEGammasNLTrue_TrueE = nullptr;
    TH2F* hNCellVsEGammasNLTrue_RecE = nullptr;

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

    TH2F* hNCellVsEGammasSideBandOnlyLowClus = nullptr;
    TH2F* hNCellVsEGammasNLSideBandOnlyLowClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueSideBandOnlyLowClus = nullptr;
    TH2F* hNCellVsEGammasNLTrueElecSideBandOnlyLowClus = nullptr;

    //--------------------------------------
    //--- TH2 NCell Vs E histograms for clusterizer studies
    //--------------------------------------
    TH2F* hNCellVsENL = nullptr;
    TH2F* hNCellVsETMNL = nullptr;
    TH2F* hNCellVsENoExotic50 = nullptr;
    TH2F* hNCellVsENoExotic75 = nullptr;
    TH2F* hNCellVsETMNLS500A100 = nullptr;
    TH2F* hNCellVsETMNLS500A105 = nullptr;
    TH2F* hNCellVsETMNLS525A105 = nullptr;
    TH2F* hNCellVsETMNLS550A110 = nullptr;
    TH2F* hNCellVsETMNLS510A102 = nullptr;
    TH2F* hNCellVsETMNLS500A110 = nullptr;

    //--------------------------------------
    //--- True particles (gamma, elec, hadrons)
    //--------------------------------------

    TH2F* hNCellVsETrueGamma = nullptr;
    TH2F* hNCellVsETrueElec = nullptr;
    TH2F* hNCellVsETrueHadr = nullptr;

    //--------------------------------------
    //--- True vs. rec. energy
    //--------------------------------------
    //
    TH2F* hRecVsTrueE = nullptr;
    TH2F* hRecDivTrueE = nullptr;
    TH2F* hRecDivTrueEOneCell = nullptr;
    TH2F* hRecDivTrueETwoCell = nullptr;
    TH2F* hRecDivTrueEThreeCell = nullptr;

    //--------------------------------------
    //--- Inv. Mass histograms
    //--------------------------------------
    TH2F* hInvMassVsPt = nullptr;
    TH2F* hInvMassVsPtBack = nullptr;
    TH2F* hInvMassVsPtGG = nullptr;
    TH2F* hInvMassVsPtGC = nullptr;
    TH2F* hInvMassVsPtCC = nullptr;
    TH2F* hInvMassVsPtHadrons = nullptr;
    TH2F* hInvMassVsHighGammaPt = nullptr;
    TH2F* hInvMassVsHighGammaPtBack = nullptr;
    TH2F* hInvMassVsBothGammaPtBack = nullptr;
    TH2F* hInvMassVsLowGammaPtBack = nullptr;
    TH2F* hInvMassVsPtHigh = nullptr;
    TH2F* hInvMassVsPtHighGamma = nullptr;
    TH2F* hInvMassVsPtHighElec = nullptr;
    TH2F* hInvMassVsPtBoth = nullptr;
    TH2F* hInvMassVsPtBothGamma = nullptr;
    TH2F* hInvMassVsPtBothElec = nullptr;
    TH2F* hInvMassVsPtPOnly1cell = nullptr;
    TH2F* hInvMassVsPtLow = nullptr;
    TH2F* hInvMassVsPtLowGamma = nullptr;
    TH2F* hInvMassVsPtLowElec = nullptr;
    TH2F* hInvMassVsPtHigh1cell = nullptr;
    TH2F* hInvMassVsPtHigh2cell = nullptr;
    TH2F* hInvMassVsPtHigh3cell = nullptr;

    TH2F* hInvMassVsPtDoubleCount = nullptr;

    TH2F* CorrelationLowHigh = nullptr;

    TH2F* hInvMassCheck = nullptr;
    //--------------------------------------
    //--- NCell vs. E for electron clusters
    //--------------------------------------
    TH2F* hNCellVsEelec = nullptr;
    TH2F* hNCellVsEelecNL = nullptr;
    TH2F* hNCellVsEelecNLTrue = nullptr;

    //--------------------------------------
    //--- E over P
    //--------------------------------------

    TH2F* hEOverPElec = nullptr;
    TH2F* hEOverPElecTrueElec = nullptr;

    //--------------------------------------
    //--- Random stuff
    //--------------------------------------
    TH2F* hHitmap = nullptr;
    TH2F* hHitmapGammas = nullptr;
    TH2F* hHitmapElec = nullptr;
    TH2F* hHitmapHadrons = nullptr;
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
    std::vector<unsigned int> vecGoodGammasSideBandOnlyLow;
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
DataTree::DataTree(){

}
//----------------------------
DataTree::~DataTree(){

}
//----------------------------
DataTree::DataTree(bool mc, TString Period, bool light, int skipping){
  isMC = mc;
  fDoLight = light;
  SkippingFrac = skipping;
  if(skipping > 0) fDoSkipping = true;
  if(Period.Contains("13TeVNomB")) period = 0;
  else if(Period.Contains("13TeVLowB")) period = 1;
  else if(Period.Contains("8TeV")) period = 2;

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

  T->SetBranchAddress("Cluster_E", &energy);
  cerr<<__LINE__<<endl;
  T->SetBranchAddress("Cluster_IsEMCAL", &isEMC);
  T->SetBranchAddress("Cluster_Eta", &eta);
  T->SetBranchAddress("Cluster_Phi", &phi);
  T->SetBranchAddress("Cluster_Cells_E", CellEnergy);
  T->SetBranchAddress("Cluster_NumCells", &Cluster_NumCells);
  T->SetBranchAddress("Cluster_Cells_ID", Cluster_CellsID);
  T->SetBranchAddress("Event_Vertex_X", &vertexZ);

  if(isMC)T->SetBranchAddress("Cluster_MC_FirstLabel", &MCLabel);
  if(isMC)T->SetBranchAddress("Cluster_MC_TrueEFirstLabel", &MCTrueEnergy);



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
    // T->SetBranchAddress("Mother_MC_Label", &ClusterPDG);
    T->SetBranchAddress("Cluster_MC_FirstLabel", &ClusterPDG);
  }
}

//--------------------------------------
//--- Setting up the histos
//--------------------------------------
void DataTree::CreateHistos(){

  cerr<<__LINE__<<endl;
  hNCellVsEGammas = new TH2F("EvsNCell", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNL = new TH2F("EvsNCellNL", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrue = new TH2F("hNCellVsEGammasNLTrue", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElec = new TH2F("hNCellVsEGammasNLTrueElec", "", nBinsE, arrEbins , 20, -0.5, 19.5);

  hNCellVsEGammasNLTrue_TrueE = new TH2F("hNCellVsEGammasNLTrue_TrueE", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrue_RecE = new TH2F("hNCellVsEGammasNLTrue_RecE", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  cerr<<__LINE__<<endl;
  hNCellVsEGammasLeft = new TH2F("EvsNCellLeft", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLLeft = new TH2F("EvsNCellNLLeft", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueLeft = new TH2F("hNCellVsEGammasNLTrueLeft", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElecLeft = new TH2F("hNCellVsEGammasNLTrueElecLeft", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  cerr<<__LINE__<<endl;
  hNCellVsEGammasWide = new TH2F("EvsNCellWide", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLWide = new TH2F("EvsNCellNLWide", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueWide = new TH2F("hNCellVsEGammasNLTrueWide", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElecWide = new TH2F("hNCellVsEGammasNLTrueElecWide", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  cerr<<__LINE__<<endl;
  hNCellVsEGammasSideBand = new TH2F("EvsNCellSideBand", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLSideBand = new TH2F("EvsNCellNLSideBand", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueSideBand = new TH2F("hNCellVsEGammasNLTrueSideBand", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElecSideBand = new TH2F("hNCellVsEGammasNLTrueElecSideBand", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  cerr<<__LINE__<<endl;
  hNCellVsEGammasOnlyHighClus = new TH2F("hNCellVsEGammasOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLOnlyHighClus = new TH2F("hNCellVsEGammasNLOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueOnlyHighClus = new TH2F("hNCellVsEGammasNLTrueOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElecOnlyHighClus = new TH2F("hNCellVsEGammasNLTrueElecOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  cerr<<__LINE__<<endl;
  hNCellVsEGammasOnlyLowClus = new TH2F("hNCellVsEGammasOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLOnlyLowClus = new TH2F("hNCellVsEGammasNLOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueOnlyLowClus = new TH2F("hNCellVsEGammasNLTrueOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElecOnlyLowClus = new TH2F("hNCellVsEGammasNLTrueElecOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  cerr<<__LINE__<<endl;
  hNCellVsEGammasSideBandOnlyHighClus = new TH2F("hNCellVsEGammasSideBandOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLSideBandOnlyHighClus = new TH2F("hNCellVsEGammasNLSideBandOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueSideBandOnlyHighClus = new TH2F("hNCellVsEGammasNLTrueSideBandOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElecSideBandOnlyHighClus = new TH2F("hNCellVsEGammasNLTrueElecSideBandOnlyHighClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  cerr<<__LINE__<<endl;
  hNCellVsEGammasSideBandOnlyLowClus = new TH2F("hNCellVsEGammasSideBandOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLSideBandOnlyLowClus = new TH2F("hNCellVsEGammasNLSideBandOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueSideBandOnlyLowClus = new TH2F("hNCellVsEGammasNLTrueSideBandOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEGammasNLTrueElecSideBandOnlyLowClus = new TH2F("hNCellVsEGammasNLTrueElecSideBandOnlyLowClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);

  hNCellVsENL = new TH2F("hNCellVsENL_AllClus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETMNL = new TH2F("hNCellVsETMNL_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsENoExotic50 = new TH2F("hNCellVsENoExotic50_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsENoExotic75 = new TH2F("hNCellVsENoExotic75_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETMNLS500A100 = new TH2F("hNCellVsETMNLS500A100_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETMNLS500A105 = new TH2F("hNCellVsETMNLS500A105_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETMNLS525A105 = new TH2F("hNCellVsETMNLS525A105_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETMNLS550A110 = new TH2F("hNCellVsETMNLS550A110_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETMNLS510A102 = new TH2F("hNCellVsETMNLS510A102_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETMNLS500A110 = new TH2F("hNCellVsETMNLS500A110_AllCLus", "", nBinsE, arrEbins , 20, -0.5, 19.5);


  hNCellVsETrueGamma = new TH2F("hNCellVsETrueGamma", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETrueElec = new TH2F("hNCellVsETrueElec", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsETrueHadr = new TH2F("hNCellVsETrueHadr", "", nBinsE, arrEbins , 20, -0.5, 19.5);

  hRecVsTrueE = new TH2F("hRecVsTrueE", "", 100, 0, 10, 100, 0, 10);
  hRecDivTrueE = new TH2F("hRecDivTrueE", "", 100, 0, 10, 200, -2, 2);
  hRecDivTrueEOneCell = new TH2F("hRecDivTrueEOneCell", "", 100, 0, 10, 200, -2, 2);
  hRecDivTrueETwoCell = new TH2F("hRecDivTrueETwoCell", "", 100, 0, 10, 200, -2, 2);
  hRecDivTrueEThreeCell = new TH2F("hRecDivTrueEThreeCell", "", 100, 0, 10, 200, -2, 2);

  hInvMassVsPt = new TH2F("hInvMassVsPt", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtGG = new TH2F("hInvMassVsPtGG", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtGC = new TH2F("hInvMassVsPtGC", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtCC = new TH2F("hInvMassVsPtCC", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtHadrons = new TH2F("hInvMassVsPtHadrons", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtBack = new TH2F("hInvMassVsPtBack", "", 500, 0, 1., 200, 0, 20);

  hInvMassVsPtHigh = new TH2F("hInvMassVsPtHigh", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtHighGamma = new TH2F("hInvMassVsPtHighGamma", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtHighElec = new TH2F("hInvMassVsPtHighElec", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtHigh1cell = new TH2F("hInvMassVsPtHigh1cell", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtHigh2cell = new TH2F("hInvMassVsPtHigh2cell", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtHigh3cell = new TH2F("hInvMassVsPtHigh3cell", "", 500, 0, 1., 200, 0, 20);

  hInvMassVsPtBoth = new TH2F("hInvMassVsPtBoth", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtBothGamma = new TH2F("hInvMassVsPtBothGamma", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtBothElec = new TH2F("hInvMassVsPtBothElec", "", 500, 0, 1., 200, 0, 20);

  hInvMassVsPtPOnly1cell = new TH2F("hInvMassVsPtPOnly1cell", "", 500, 0, 1., 200, 0, 20);

  hInvMassVsPtDoubleCount = new TH2F("hInvMassVsPtDoubleCount", "", 500, 0, 1., 200, 0, 20);

  hInvMassCheck = new TH2F("hInvMassCheck", "", 500, 0, 1., 200, 0, 20);

  CorrelationLowHigh = new TH2F("CorrelationLowHigh", "", 200, 0, 20, 200, 0, 20);

  hInvMassVsPtLowGamma = new TH2F("hInvMassVsPtLowGamma", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtLowElec = new TH2F("hInvMassVsPtLowElec", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsPtLow = new TH2F("hInvMassVsPtLow", "", 500, 0, 1., 200, 0, 20);

  hInvMassVsHighGammaPt = new TH2F("hInvMassVsHighGammaPt", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsHighGammaPtBack = new TH2F("hInvMassVsHighGammaPtBack", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsBothGammaPtBack = new TH2F("hInvMassVsBothGammaPtBack", "", 500, 0, 1., 200, 0, 20);
  hInvMassVsLowGammaPtBack = new TH2F("hInvMassVsLowGammaPtBack", "", 500, 0, 1., 200, 0, 20);
  hHitmap = new TH2F("hHitmap", "", 160, -0.8, 0.8, 126*2, 0, 6.2);
  hHitmapGammas = new TH2F("hHitmapGammas", "", 160, -0.8, 0.8, 126*2, 0, 6.2);
  hHitmapElec = new TH2F("hHitmapElec", "", 160, -0.8, 0.8, 126*2, 0, 6.2);
  hHitmapHadrons = new TH2F("hHitmapHadrons", "", 160, -0.8, 0.8, 126*2, 0, 6.2);
  hClusE = new TH1F("hClusE", "", nBinsE, arrEbins );
  hClusENL = new TH1F("hClusENL", "", nBinsE, arrEbins );

  hNCellVsEelec = new TH2F("hNCellVsEelec", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEelecNL = new TH2F("hNCellVsEelecNL", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hNCellVsEelecNLTrue = new TH2F("hNCellVsEelecNLTrue", "", nBinsE, arrEbins , 20, -0.5, 19.5);
  hEOverPElec = new TH2F("hEOverPElec", "" , 100, 0, 10, 100, 0.5, 2);
  hEOverPElecTrueElec = new TH2F("hEOverPElecTrueElec", "" , 100, 0, 10, 100, 0.5, 2);

}

Cluster *DataTree::GetCluster(unsigned int i){
  T->GetEvent(i);
  // for(int k = 0; k < Surrounding_Tracks_N; k++){
  //   cout<<"Surrounding_Cells_ID["<<k<<"] = "<<Surrounding_Cells_ID[k]<<endl;
  // }
  Cluster *tmp;
  if(!isMC){
      tmp = new Cluster(applyNL(energy, Cluster_NumCells), eta, phi, Cluster_NumCells , CellEnergy, Cluster_CellsID,
      Surrounding_Cells_N, Surrounding_Cells_E, Surrounding_Cells_R, Surrounding_Cells_ID,
      Surrounding_Tracks_N, Surrounding_Tracks_P, Surrounding_Tracks_R, Surrounding_Tracks_nSigdEdx, Surrounding_Tracks_IsV0,
      false, 0, 0);
  } else {
      tmp = new Cluster(applyNL(energy, Cluster_NumCells), eta, phi, Cluster_NumCells , CellEnergy, Cluster_CellsID,
      Surrounding_Cells_N, Surrounding_Cells_E, Surrounding_Cells_R, Surrounding_Cells_ID,
      Surrounding_Tracks_N, Surrounding_Tracks_P, Surrounding_Tracks_R, Surrounding_Tracks_nSigdEdx, Surrounding_Tracks_IsV0,
      isMC, ClusterPDG, MCLabel, MCTrueEnergy);
  }
  return tmp;
}

//--------------------------------------
//--- get the number of clusters in current event, use he vertex position as an indicator
//--------------------------------------
unsigned int DataTree::GetNClusInEvt(unsigned int i){
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
void DataTree::NextEvt(){
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
void DataTree::Process(unsigned int maxNumClus){
  cerr<<"max events: "<<T->GetEntries()<<endl;

  CreateHistos();
  FillMassHists("fMassPos.root");
  // unsigned int counter = 0;
  if(fdebug == 2) cerr<<__LINE__<<endl;
  while(currEvt < maxNumClus && currEvt < T->GetEntries()*0.99){
    // counter++;
    if(currEvt%10000 == 0) cerr<<currEvt<<" events processed..."<<endl;
    // if(fdebug == 2) cerr<<__LINE__<<endl;
    NextEvt();
    if(Evt.NClus < 2) continue;
    // if(fdebug)cout<<"energy: "<<energy<<endl;
    // if(fdebug == 2) cerr<<"In Process: "<<__LINE__<<endl;
    GetGammasViaPi0();
    // if(fdebug == 2) cerr<<"In Process: "<<__LINE__<<endl;
    FillNCellVsEGammas();
    // if(fdebug == 2) cerr<<"In Process: "<<__LINE__<<endl;
    FillHitmap();
    // if(fdebug == 2) cerr<<"In Process: "<<__LINE__<<endl;
    FillElectronHist();
    // if(fdebug == 2) cerr<<"In Process: "<<__LINE__<<endl;
    FillAllCLusterHist();
    // if(fdebug == 2) cerr<<"In Process: "<<__LINE__<<endl;
  }

}

//--------------------------------------
//--- Eta to theta
//--------------------------------------
inline float DataTree::EtaToTheta(float Eta){
 return (float) 2.*atan(exp(-Eta));
}

//--------------------------------------
//--- Get Momntum vector from cluster
//--------------------------------------
std::array<float, 3> DataTree::GetMomVect(unsigned int i){
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
void DataTree::GetGammasViaPi0(){
  vecGoodGammas.clear();
  // vecGoodGammas.resize(0);
  vecGoodGammasLeft.clear();
  // vecGoodGammasLeft.resize(0);
  vecGoodGammasSideBand.clear();
  // vecGoodGammasSideBand.resize(0);
  vecGoodGammasWide.clear();
  // vecGoodGammasWide.resize(0);
  vecGoodGammasOnlyHigh.clear();
  vecGoodGammasOnlyLow.clear();
  // vecGoodGammasOnlyHigh.resize(0);
  vecGoodGammasSideBandOnlyHigh.clear();
  vecGoodGammasSideBandOnlyLow.clear();

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
    if(!SMwithNoTRD(ig1)) continue;
    if(IsBehindTRDSupport(ig1)) continue;
    if(fdebug) cerr<<__LINE__<<endl;
    if(IsBorderCell(ig1, 3)) continue;
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
      if(!SMwithNoTRD(ig2)) continue;
      if(IsBehindTRDSupport(ig2)) continue;
      if(IsBorderCell(ig2, 3)) continue;

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

      // check the opening angle
      if(LVg1.Angle(LVg2.Vect()) < 0.0202) continue;

      // cerr<<Pi0.M()<<endl;

      hInvMassVsPt->Fill(Pi0.M(), Pi0.Pt());
      // both gammas
      // get high and low gamma pT
      float gammapTHigh = (LVg2.Pt() > LVg1.Pt() ? LVg2.Pt() : LVg1.Pt());
      float gammapTLow = (LVg2.Pt() > LVg1.Pt() ? LVg1.Pt() : LVg2.Pt());
      if(PDGig1 == PDGig2 && PDGig2 == 22) {
        hInvMassVsPtGG->Fill(Pi0.M(), Pi0.Pt());
      }
      else if( (fabs(PDGig1) == 11 && PDGig2 == 22) || (PDGig1 == 22 && fabs(PDGig2) == 11)) {
        hInvMassVsPtGC->Fill(Pi0.M(), Pi0.Pt());
      }
      else if((PDGig1 == PDGig2) && (fabs(PDGig2) == 11)) {
        hInvMassVsPtCC->Fill(Pi0.M(), Pi0.Pt());
      }
      else{
        hInvMassVsPtHadrons->Fill(Pi0.M(), Pi0.Pt());
      }

      hInvMassVsHighGammaPt->Fill(Pi0.M(), gammapTHigh);
      hInvMassVsPtHigh->Fill(Pi0.M(), gammapTHigh);
      hInvMassVsPtBoth->Fill(Pi0.M(), gammapTHigh);
      hInvMassVsPtLow->Fill(Pi0.M(), gammapTLow);
      hInvMassVsPtBoth->Fill(Pi0.M(), gammapTLow);

      if(Pi0.M() > 0.09 && Pi0.M() < 0.17){
        CorrelationLowHigh->Fill(gammapTHigh, gammapTLow);
      }




      if(!fDoLight){

        unsigned int iclusHigh = (LVg2.Pt() > LVg1.Pt()) ? ig2 : ig1;
        unsigned int ClusNCells = Evt.clus[iclusHigh]->ncells;
        if(ClusNCells == 1) {
          hInvMassVsPtHigh1cell->Fill(Pi0.M(), gammapTHigh);
          hInvMassVsPtPOnly1cell->Fill(Pi0.M(), Pi0.Pt());
        }
        if(ClusNCells == 2) hInvMassVsPtHigh2cell->Fill(Pi0.M(), gammapTHigh);
        if(ClusNCells > 2) hInvMassVsPtHigh3cell->Fill(Pi0.M(), gammapTHigh);

        if(isMC){
          // cerr<<"PDGig1: "<<PDGig1<<"   PDGig2"<<PDGig2<<endl;
          int iPDGHigh = (LVg2.Pt() > LVg1.Pt()) ? PDGig2 : PDGig1;
          int iPDGLow = (LVg2.Pt() < LVg1.Pt()) ? PDGig2 : PDGig1;
          if(iPDGHigh == 22) {
            hInvMassVsPtHighGamma->Fill(Pi0.M(), gammapTHigh);
            hInvMassVsPtBothGamma->Fill(Pi0.M(), gammapTHigh);
          }
          else if(fabs(iPDGHigh) == 11) {
            hInvMassVsPtHighElec->Fill(Pi0.M(), gammapTHigh);
            hInvMassVsPtBothElec->Fill(Pi0.M(), gammapTHigh);
          }
          if(iPDGLow == 22) {
            hInvMassVsPtLowGamma->Fill(Pi0.M(), gammapTLow);
            hInvMassVsPtBothGamma->Fill(Pi0.M(), gammapTLow);
          }
          else if(fabs(iPDGLow) == 11) {
            hInvMassVsPtLowElec->Fill(Pi0.M(), gammapTLow);
            hInvMassVsPtBothElec->Fill(Pi0.M(), gammapTLow);
          }

        }
      }
      if(fdebug) cerr<<__LINE__<<endl;
      // check if pi0
      float massPos = 0.12;
      // if(grMass) massPos = grMass->Eval(Pi0.Pt());
      if(fMassPos) massPos = fMassPos->Eval(Pi0.Pt());
      // cerr<<Pi0.Pt()<<"\t"<<massPos<<endl;
      //-----------------------------------------------
      // Right side of peak
      //-----------------------------------------------
      if(!fDoLight){
        if(Pi0.M() > massPos && Pi0.M() < massPos + 0.02){
        // if(Pi0.M() > massPos && Pi0.M() < massPos + 0.02){
          vecGoodGammas.push_back(ig1);
          vecGoodGammas.push_back(ig2);
        }
        //-----------------------------------------------
        // left side of peak
        //-----------------------------------------------
        if(Pi0.M() > massPos - 0.05 && Pi0.M() < massPos){
        // if(Pi0.M() > 0.1 && Pi0.M() < 0.6){
          vecGoodGammasLeft.push_back(ig1);
          vecGoodGammasLeft.push_back(ig2);
        }
      }
      //-----------------------------------------------
      // both side of peak
      //-----------------------------------------------
      if(Pi0.M() > 0.09 && Pi0.M() < 0.17){
      // if(Pi0.M() > massPos - 0.05 && Pi0.M() < massPos + 0.02){
        vecGoodGammasWide.push_back(ig1);
        vecGoodGammasWide.push_back(ig2);

        // only fill higher cluster (suppress conversions?)
        if(Eg2 > Eg1) {
          vecGoodGammasOnlyHigh.push_back(ig2);
          vecGoodGammasOnlyLow.push_back(ig1);
        }
        else{
          vecGoodGammasOnlyHigh.push_back(ig1);
          vecGoodGammasOnlyLow.push_back(ig2);
        }

        // Double counting
        if(std::find(fDoubleCount.begin(), fDoubleCount.end(), ig1) != fDoubleCount.end()){
          hInvMassVsPtDoubleCount->Fill(Pi0.M(), LVg1.Pt());
        }
        if(std::find(fDoubleCount.begin(), fDoubleCount.end(), ig2) != fDoubleCount.end()){
          hInvMassVsPtDoubleCount->Fill(Pi0.M(), LVg2.Pt());
        }
        fDoubleCount.push_back(ig1);
        fDoubleCount.push_back(ig2);

      }

      //-----------------------------------------------
      // right side sideband
      //-----------------------------------------------
      if(Pi0.M() > 0.2 && Pi0.M() < 0.3){
      // if(Pi0.M() > massPos + 0.05 && Pi0.M() < massPos + 0.2){
        vecGoodGammasSideBand.push_back(ig1);
        vecGoodGammasSideBand.push_back(ig2);
        if(Eg2 > Eg1) {
          vecGoodGammasSideBandOnlyHigh.push_back(ig2);
          vecGoodGammasSideBandOnlyLow.push_back(ig1);
        }
        else {
          vecGoodGammasSideBandOnlyHigh.push_back(ig1);
          vecGoodGammasSideBandOnlyLow.push_back(ig2);
        }
      }

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

        bool isAcc1 = (LVg1.Angle(LVg3.Vect()) > 0.0202) ? true : false;
        bool isAcc2 = (LVg2.Angle(LVg3.Vect()) > 0.0202) ? true : false;
        // cerr<<Pi0.M()<<endl;

        if(isAcc1) hInvMassVsPtBack->Fill(BackPi01.M(), BackPi01.Pt());
        if(isAcc2) hInvMassVsPtBack->Fill(BackPi02.M(), BackPi02.Pt());

        if(isAcc1) {
          hInvMassVsHighGammaPtBack->Fill(BackPi01.M(), (LVg1.Pt() > LVg3.Pt() ? LVg1.Pt() : LVg3.Pt()));
          hInvMassVsLowGammaPtBack->Fill(BackPi01.M(), (LVg1.Pt() > LVg3.Pt() ? LVg3.Pt() : LVg1.Pt()));
          hInvMassVsBothGammaPtBack->Fill(BackPi01.M(), (LVg1.Pt() > LVg3.Pt() ? LVg3.Pt() : LVg1.Pt()));
          hInvMassVsBothGammaPtBack->Fill(BackPi01.M(), (LVg1.Pt() > LVg3.Pt() ? LVg1.Pt() : LVg3.Pt()));
        }
        if(isAcc2) {
          hInvMassVsHighGammaPtBack->Fill(BackPi02.M(), (LVg2.Pt() > LVg3.Pt() ? LVg2.Pt() : LVg3.Pt()));
          hInvMassVsLowGammaPtBack->Fill(BackPi02.M(), (LVg2.Pt() > LVg3.Pt() ? LVg3.Pt() : LVg2.Pt() ));
          hInvMassVsBothGammaPtBack->Fill(BackPi02.M(), (LVg2.Pt() > LVg3.Pt() ? LVg2.Pt() : LVg3.Pt()));
          hInvMassVsBothGammaPtBack->Fill(BackPi02.M(), (LVg2.Pt() > LVg3.Pt() ? LVg3.Pt() : LVg2.Pt() ));
        }

      }
    }
  }
}

//--------------------------------------
//--- Fill the histograms from the vectors for each event
//--------------------------------------
void DataTree::FillNCellVsEGammas(){
  if(!fDoLight){
    sort( vecGoodGammas.begin(), vecGoodGammas.end() );
    vecGoodGammas.erase( unique( vecGoodGammas.begin(), vecGoodGammas.end() ), vecGoodGammas.end() );
    for(auto &i : vecGoodGammas){
    // for(unsigned int i = 0; i < vecGoodGammas.size(); ++i){
      // T->GetEvent(vecGoodGammas.at(i));
      hNCellVsEGammas->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      hNCellVsEGammasNL->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(isMC){
        if(Evt.clus[i]->ClusterPDG == 22) hNCellVsEGammasNLTrue->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
        if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElec->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      }
    }
    sort( vecGoodGammasLeft.begin(), vecGoodGammasLeft.end() );
    vecGoodGammasLeft.erase( unique( vecGoodGammasLeft.begin(), vecGoodGammasLeft.end() ), vecGoodGammasLeft.end() );
    // for(unsigned int i = 0; i < vecGoodGammasLeft.size(); ++i){
    for(auto &i : vecGoodGammasLeft){
      // T->GetEvent(vecGoodGammasLeft.at(i));
      hNCellVsEGammasLeft->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      hNCellVsEGammasNLLeft->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(isMC){
        if(Evt.clus[i]->ClusterPDG == 22) hNCellVsEGammasNLTrueLeft->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
        if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElecLeft->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      }
    }
  }
  sort( vecGoodGammasSideBand.begin(), vecGoodGammasSideBand.end() );
  vecGoodGammasSideBand.erase( unique( vecGoodGammasSideBand.begin(), vecGoodGammasSideBand.end() ), vecGoodGammasSideBand.end() );
  // for(unsigned int i = 0; i < vecGoodGammasSideBand.size(); ++i){
  for(auto &i : vecGoodGammasSideBand){
    // T->GetEvent(vecGoodGammasSideBand.at(i));
    hNCellVsEGammasSideBand->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    hNCellVsEGammasNLSideBand->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    if(isMC){
      if(Evt.clus[i]->ClusterPDG == 22) hNCellVsEGammasNLTrueSideBand->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElecSideBand->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    }
  }
  sort( vecGoodGammasWide.begin(), vecGoodGammasWide.end() );
  vecGoodGammasWide.erase( unique( vecGoodGammasWide.begin(), vecGoodGammasWide.end() ), vecGoodGammasWide.end() );
  // for(unsigned int i = 0; i < vecGoodGammasWide.size(); ++i){
    for(auto &i : vecGoodGammasWide){
    // T->GetEvent(vecGoodGammasWide.at(i));
    hNCellVsEGammasWide->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    hNCellVsEGammasNLWide->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    if(isMC){
      if(Evt.clus[i]->ClusterPDG == 22) hNCellVsEGammasNLTrueWide->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElecWide->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    }
  }


    //
    // for(int i = 0; i < vecGoodGammasWide.size(); ++i){
    // for(int j = i+1; j < vecGoodGammasWide.size(); ++j){
    //   int ig1 = vecGoodGammasWide.at(i);
    //   int ig2 = vecGoodGammasWide.at(j);
    //   // write them into TLorentzvector
    //   std::array<float, 3> Pg = GetMomVect(ig1);
    //   TLorentzVector LVg1;
    //   LVg1.SetX(Pg[0]);//, Pg1[1], Pg1[2], Eg1);
    //   LVg1.SetY(Pg[1]);//, Pg1[1], Pg1[2], Eg1);
    //   LVg1.SetZ(Pg[2]);//, Pg1[1], Pg1[2], Eg1);
    //   LVg1.SetE(Evt.clus[ig1]->e);
    //   std::array<float, 3> Pg2 = GetMomVect(ig2);
    //   TLorentzVector LVg2;
    //   LVg2.SetX(Pg2[0]);//, Pg1[1], Pg1[2], Eg1);
    //   LVg2.SetY(Pg2[1]);//, Pg1[1], Pg1[2], Eg1);
    //   LVg2.SetZ(Pg2[2]);//, Pg1[1], Pg1[2], Eg1);
    //   LVg2.SetE(Evt.clus[ig2]->e);
    //   TLorentzVector LPi0 = (LVg1 + LVg2);
    //
    //   hInvMassCheck->Fill(LPi0.M(), LPi0.Pt());
    //
    // }
    // }
  sort( vecGoodGammasOnlyHigh.begin(), vecGoodGammasOnlyHigh.end() );
  vecGoodGammasOnlyHigh.erase( unique( vecGoodGammasOnlyHigh.begin(), vecGoodGammasOnlyHigh.end() ), vecGoodGammasOnlyHigh.end() );
  // for(unsigned int i = 0; i < vecGoodGammasOnlyHigh.size(); ++i){
    for(auto &i : vecGoodGammasOnlyHigh){
    // T->GetEvent(vecGoodGammasOnlyHigh.at(i));
    hNCellVsEGammasOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    hNCellVsEGammasNLOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    if(isMC){
      if(Evt.clus[i]->ClusterPDG == 22) hNCellVsEGammasNLTrueOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElecOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    }
  }
  sort( vecGoodGammasOnlyLow.begin(), vecGoodGammasOnlyLow.end() );
  vecGoodGammasOnlyLow.erase( unique( vecGoodGammasOnlyLow.begin(), vecGoodGammasOnlyLow.end() ), vecGoodGammasOnlyLow.end() );
  // for(unsigned int i = 0; i < vecGoodGammasOnlyLow.size(); ++i){
    for(auto &i : vecGoodGammasOnlyLow){
    // T->GetEvent(vecGoodGammasOnlyHigh.at(i));
    hNCellVsEGammasOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    hNCellVsEGammasNLOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    if(isMC){
      if(Evt.clus[i]->ClusterPDG == 22) hNCellVsEGammasNLTrueOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElecOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    }
  }
  sort( vecGoodGammasSideBandOnlyHigh.begin(), vecGoodGammasSideBandOnlyHigh.end() );
  vecGoodGammasSideBandOnlyHigh.erase( unique( vecGoodGammasSideBandOnlyHigh.begin(), vecGoodGammasSideBandOnlyHigh.end() ), vecGoodGammasSideBandOnlyHigh.end() );
  // for(unsigned int i = 0; i < vecGoodGammasSideBandOnlyHigh.size(); ++i){
    for(auto &i : vecGoodGammasSideBandOnlyHigh){
    // T->GetEvent(vecGoodGammasSideBandOnlyHigh.at(i));
    hNCellVsEGammasSideBandOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    hNCellVsEGammasNLSideBandOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    if(isMC){
      if(fabs(Evt.clus[i]->ClusterPDG == 22)) hNCellVsEGammasNLTrueSideBandOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElecSideBandOnlyHighClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    }
  }
  sort( vecGoodGammasSideBandOnlyLow.begin(), vecGoodGammasSideBandOnlyLow.end() );
  vecGoodGammasSideBandOnlyLow.erase( unique( vecGoodGammasSideBandOnlyLow.begin(), vecGoodGammasSideBandOnlyLow.end() ), vecGoodGammasSideBandOnlyLow.end() );
  // for(unsigned int i = 0; i < vecGoodGammasSideBandOnlyLow.size(); ++i){
    for(auto &i : vecGoodGammasSideBandOnlyLow){
    // T->GetEvent(vecGoodGammasSideBandOnlyHigh.at(i));
    hNCellVsEGammasSideBandOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    hNCellVsEGammasNLSideBandOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    if(isMC){
      if(fabs(Evt.clus[i]->ClusterPDG == 22)) hNCellVsEGammasNLTrueSideBandOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
      if(fabs(Evt.clus[i]->ClusterPDG) == 11) hNCellVsEGammasNLTrueElecSideBandOnlyLowClus->Fill(Evt.clus[i]->e, Evt.clus[i]->ncells);
    }
  }

}

//--------------------------------------
//--- Hitmap
//--------------------------------------
void DataTree::FillHitmap(){
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
    if(isMC){
      if(Evt.clus[i]->ClusterPDG == 22) hHitmapGammas->Fill(theta - TMath::Pi()*0.5, Evt.clus[i]->phi);
      else if(fabs(Evt.clus[i]->ClusterPDG) == 11)hHitmapElec->Fill(theta - TMath::Pi()*0.5, Evt.clus[i]->phi);
      else hHitmapHadrons->Fill(theta - TMath::Pi()*0.5, Evt.clus[i]->phi);

    }
  }
}

//--------------------------------------
// return true if cell is isolated, false if its not isolated
//--------------------------------------
bool DataTree::IsIsolated(unsigned int clus){
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
      if(Evt.clus[clus]->surrCellsR[i] < 0.02){
        return false;
      }
    }
  }
  return true;
}

//--------------------------------------
// return true if cell is less than NCells away from SM border
//--------------------------------------
bool DataTree::IsBorderCell(unsigned int clus, unsigned int NCellDist){
  // T->GetEvent(clus);
  if(!fRejectBorderCells) return false;
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
bool DataTree::IsTrackMatched(unsigned int clus){
  // T->GetEvent(clus);
  for(int i = 0; i <Evt.clus[clus]->NsurrTracks; ++i){
    // if(Surrounding_Tracks_P[i] > 0.1){
      if(Evt.clus[clus]->surrTracksR[i] < 0.08){
        return true;
      // }
    }
  }
  return false;
}

//-----------------------------------------------------
// Fill pure electron hists (V0 tracks)
//-----------------------------------------------------
void DataTree::FillElectronHist(){

  for(unsigned int ic = 0; ic < Evt.NClus; ++ic){
    // T->GetEvent(currEvt + ic);

    unsigned int ntrack = 0;
    for(int i = 0; i <Evt.clus[ic]->NsurrTracks; ++i){
      if(Evt.clus[ic]->surrTracksR[i] < 12./430){
        ntrack++;
      }
    }
    for(int i = 0; i <Evt.clus[ic]->NsurrTracks; ++i){
      if(Evt.clus[ic]->surrTracksR[i] < 6./430){
        // has to be in electron band and cluster e and track p have to match within 10%
        if(fabs(Evt.clus[ic]->surrTracksnSig[i]) < 3 /*&& fabs((Evt.clus[ic]->surrTracksP[i]/Evt.clus[ic]->e) - 1) < 0.1*/  ){
          if(Evt.clus[ic]->surrTracksIsV0[i] == true){
            // if(fabs(Surrounding_Tracks_nSigdEdx[i]) < 2){
              // hNCellVsEelec->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
              hNCellVsEelecNL->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
            // }
            if(isMC){
              // cout<<"ntrack: "<<ntrack<<endl;
              // cout<<"Evt.clus[ic]->ClusterPDG: "<<Evt.clus[ic]->ClusterPDG<<endl;
              // cout<<"Evt.clus[ic]->surrTracksP[i]/Evt.clus[ic]->e = "<<Evt.clus[ic]->surrTracksP[i]/Evt.clus[ic]->e<<endl;
              // cout<<"Evt.clus[ic]->surrTracksnSig[i] = "<<Evt.clus[ic]->surrTracksnSig[i]<<endl;
              if(fabs(Evt.clus[ic]->ClusterPDG) == 11){
                hNCellVsEelecNLTrue->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
              }
            }
            hEOverPElec->Fill(Evt.clus[ic]->e, Evt.clus[ic]->e / Evt.clus[ic]->surrTracksP[i]);
            if(isMC){
              if(fabs(Evt.clus[ic]->ClusterPDG) == 11){
                hEOverPElecTrueElec->Fill(Evt.clus[ic]->e, Evt.clus[ic]->e / Evt.clus[ic]->surrTracksP[i]);
              }
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------
// Check if the cluster has a cell with at least eThresh next to it
// If not its considered an exotic cluster
//-----------------------------------------------------
bool DataTree::IsExoticClus(unsigned int ic, float eThresh){

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


bool DataTree::SMwithNoTRD(unsigned int clus){
  // only check if analysis with no TRD is needed
  if(!fDoNoTRD) return true;
  // only look for 8 TeV
  if(period != 2) return true;

  // SM 0, 1, 2 and 3 have no TRD in front for pp 8TeV
  if(Evt.clus[clus]->SM < 4) return true;
  else return false;

}

bool DataTree::IsBehindTRDSupport(unsigned int clus){
  // only check if analysis with no TRD is needed
  if(!fDoRejectTRDSupport) return false;

  float clusEta = Evt.clus[clus]->eta;
  if(fabs(clusEta) < 0.54 && fabs(clusEta) > 0.49) return true;
  else if(fabs(clusEta) < 0.19 && fabs(clusEta) > 0.14) return true;
  else return false;

}

//-----------------------------------------------------
// FIll all clusters
//-----------------------------------------------------
void DataTree::FillAllCLusterHist(){

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
      if(IsClusAcceptedByThreshold(ic, 0.525, 0.105)) hNCellVsETMNLS525A105->Fill(Eclus, Evt.clus[ic]->ncells);
      if(IsClusAcceptedByThreshold(ic, 0.51, 0.102)) hNCellVsETMNLS510A102->Fill(Eclus, Evt.clus[ic]->ncells);
      if(IsClusAcceptedByThreshold(ic, 0.55, 0.11)) hNCellVsETMNLS550A110->Fill(Eclus, Evt.clus[ic]->ncells);

      if(!IsExoticClus(ic, 0.05)){
        hNCellVsENoExotic50->Fill(Eclus, Evt.clus[ic]->ncells);
      }
      if(!IsExoticClus(ic, 0.075)){
        hNCellVsENoExotic75->Fill(Eclus, Evt.clus[ic]->ncells);
      }

      // FIll rec vs true energy for true photons
      if(isMC){
        if( Evt.clus[ic]->ClusterPDG == 22){
          hRecVsTrueE->Fill(Evt.clus[ic]->TrueE, Evt.clus[ic]->e);
          if(Evt.clus[ic]->TrueE > 0){
            hRecDivTrueE->Fill( Evt.clus[ic]->TrueE, Evt.clus[ic]->e / Evt.clus[ic]->TrueE);
            if(Evt.clus[ic]->ncells == 1)  hRecDivTrueEOneCell->Fill( Evt.clus[ic]->TrueE, Evt.clus[ic]->e / Evt.clus[ic]->TrueE);
            else if(Evt.clus[ic]->ncells == 2)  hRecDivTrueETwoCell->Fill( Evt.clus[ic]->TrueE, Evt.clus[ic]->e / Evt.clus[ic]->TrueE);
            else if(Evt.clus[ic]->ncells == 3)  hRecDivTrueEThreeCell->Fill( Evt.clus[ic]->TrueE, Evt.clus[ic]->e / Evt.clus[ic]->TrueE);
          }
          hNCellVsEGammasNLTrue_TrueE->Fill(Evt.clus[ic]->TrueE, Evt.clus[ic]->ncells);
          hNCellVsEGammasNLTrue_RecE->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
          hNCellVsETrueGamma->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
        }
        else if( fabs(Evt.clus[ic]->ClusterPDG) == 11){
          hNCellVsETrueElec->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
        }
        else {
          hNCellVsETrueHadr->Fill(Evt.clus[ic]->e, Evt.clus[ic]->ncells);
        }
      }

    }
  }
}

//-----------------------------------------------------
// apply NonLinearity
//-----------------------------------------------------
float DataTree::applyNL(float e, int numcells){

  //--------------- pp 13TeV -------------------------
  // if (period == 0){

    if(isMC){
      e /= FunctionNL_OfficialTB_100MeV_MC_V2(e);
      // fine tuning based on gaussian fits on PCMEMC in pPb5TeV
      // e /= FunctionNL_kSDM(e, 0.987912, -2.94105, -0.273207) ;  // nom. B
      if(period == 1){
        // e /= FunctionNL_kSDM(e, 1.02451, -3.49297, -0.420027); // lowB
        // e /= FunctionNL_kSDM(e, 0.986634, -4.12191, -0.321714);
        // if(period == 1) e *= 1.035;
      } else if(period == 0){
        // e /= FunctionNL_kSDM(e, 0.987912, -2.94105, -0.273207);
      } else if(period == 2){
        // e /= FunctionNL_kSDM(e, 0.987912, -2.94105, -0.273207);
        // e /= 0.98;
      }
      if(numcells == 1){ // additional fine tuning for 1 cell clusters
        e /= FunctionNL_kSDM(e,0, -0.002069903, -0.00669839);
        // e /= 0.995;
        e /= 0.97; // shift estimated from file:///home/joshua/PCG_Software/EMCal_NCellEffi/13TeVNomB_Wide/Pi0Tagging_13TeV_nom_03_18_TrueVsRec_NoFT/pdf/TrueVsRecE.pdf
      }
    } else {
      e /= FunctionNL_OfficialTB_100MeV_Data_V2(e);
      if(numcells == 1){
        e /= 0.97; // estimated from 1cell data/MC from file rootFiles/Pi0Tagging_13TeV_nom_03_18_TrueVsRec_1cellFT
      }
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
bool DataTree::IsClusAcceptedByThreshold(unsigned int i, float agg, float thr){
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
//--- Fill Mass hists
//--------------------------------------
void DataTree::FillMassHists(TString name){
  TFile f(name);
  //-------------- 13TeV ---------------
  if(period == 0){
    if(isMC)grMass = (TGraph*) f.Get("MassPos13TeV_nomB_MC");
    else grMass = (TGraph*) f.Get("MassPos13TeV_nomB_data");
    fMassPos = (TF1*) f.Get("MassPos13TeVNomB_func");
    //-------------- 13TeV ---------------
  } else if(period == 1){
    if(isMC)grMass = (TGraph*) f.Get("MassPos13TeV_MC");
    else grMass = (TGraph*) f.Get("MassPos13TeV_data");
    fMassPos = (TF1*) f.Get("MassPos13TeVLowB_func");
    //-------------- 8TeV ---------------
  } else if (period == 2){
    if(isMC)grMass = (TGraph*) f.Get("MassPos8TeV_MC");
    else grMass = (TGraph*) f.Get("MassPos8TeV_data");
    fMassPos = (TF1*) f.Get("MassPos8TeV_func");
  }
}

//--------------------------------------
//--- Write everything to file
//--------------------------------------
void DataTree::WriteToFile(){
  gSystem->Exec("mkdir rootFiles");

  TString filename = "rootFiles/Pi0Tagging_13TeV_nom_04_26_WithTRD_WithBorderCells_1cellFT";
  // TString filename = "rootFiles/Pi0Tagging_13TeV_nom_03_18_Iso02_TBNL_noScale_1cellFT_RejectBehindTRDSupport";
  if(period ==1) filename = "rootFiles/Pi0Tagging_13TeV_low_03_16_Iso02_TBNL_noScale";
  else if(period ==2) filename = "rootFiles/Pi0Tagging_8TeV_03_12_NoTRD_NoShaper_TBNL";
  std::ifstream ftest(filename+".root");

  TFile fout(filename+".root", (ftest.good() == 0) ? "Recreate" : "Update");
  fout.cd();
  hInvMassVsPt->Write(Form("hInvMassVsPt_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtGG->Write(Form("hInvMassVsPtGG_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtGC->Write(Form("hInvMassVsPtGC_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtCC->Write(Form("hInvMassVsPtCC_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtHadrons->Write(Form("hInvMassVsPtHadrons_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hInvMassVsPtBack->Write(Form("hInvMassVsPtBack_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  // high / low gamma pt
  // hInvMassVsHighGammaPt->Write(Form("hInvMassVsGammaPtHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtLow->Write(Form("hInvMassVsPtLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtHigh->Write(Form("hInvMassVsPtHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtBoth->Write(Form("hInvMassVsPtBoth_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  // seperated to electrosn/gammas
  hInvMassVsPtLowGamma->Write(Form("hInvMassVsPtGammaLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtLowElec->Write(Form("hInvMassVsPtElecLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtHighGamma->Write(Form("hInvMassVsPtGammaHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtHighElec->Write(Form("hInvMassVsPtElecHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtBothGamma->Write(Form("hInvMassVsPtGammaBoth_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtBothElec->Write(Form("hInvMassVsPtElecBoth_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hInvMassVsPtPOnly1cell->Write(Form("hInvMassVsPtPOnly1cell_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  CorrelationLowHigh->Write(Form("CorrelationLowHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  // background
  hInvMassVsHighGammaPtBack->Write(Form("hInvMassVsGammaPtBackHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsLowGammaPtBack->Write(Form("hInvMassVsGammaPtBackLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsBothGammaPtBack->Write(Form("hInvMassVsGammaPtBackBoth_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hInvMassVsPtHigh1cell->Write(Form("hInvMassVsPtHigh1cell_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtHigh2cell->Write(Form("hInvMassVsPtHigh2cell_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtHigh3cell->Write(Form("hInvMassVsPtHigh3cell_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassVsPtDoubleCount->Write(Form("hInvMassVsPtDoubleCount_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hInvMassCheck->Write(Form("hInvMassCheck%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammas->Write(Form("hNCellVsEGammasRight_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNL->Write(Form("hNCellVsEGammasNLRight_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrue->Write(Form("hNCellVsETrueGammasNLRight_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElec->Write(Form("hNCellVsEGammasNLTrueElecRight_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasNLTrue_TrueE->Write(Form("hNCellVsEGammasNLTrue_TrueE_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrue_RecE->Write(Form("hNCellVsEGammasNLTrue_RecE_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasLeft->Write(Form("hNCellVsEGammasLeft_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLLeft->Write(Form("hNCellVsEGammasNLLeft_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueLeft->Write(Form("hNCellVsETrueGammasNLLeft_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElecLeft->Write(Form("hNCellVsEGammasNLTrueElecLeft_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasWide->Write(Form("hNCellVsEGammasWide_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLWide->Write(Form("hNCellVsEGammasNLWide_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueWide->Write(Form("hNCellVsETrueGammasNLWide_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElecWide->Write(Form("hNCellVsEGammasNLTrueElecWide_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasSideBand->Write(Form("hNCellVsEGammasSideBand_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLSideBand->Write(Form("hNCellVsEGammasNLSideBand_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueSideBand->Write(Form("hNCellVsETrueGammasNLSideBand_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElecSideBand->Write(Form("hNCellVsEGammasNLTrueElecSideBand_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasOnlyHighClus->Write(Form("hNCellVsEGammasHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLOnlyHighClus->Write(Form("hNCellVsEGammasNLHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueOnlyHighClus->Write(Form("hNCellVsETrueGammasNLHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElecOnlyHighClus->Write(Form("hNCellVsEGammasNLTrueElecHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasOnlyLowClus->Write(Form("hNCellVsEGammasLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLOnlyLowClus->Write(Form("hNCellVsEGammasNLLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueOnlyLowClus->Write(Form("hNCellVsETrueGammasNLLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElecOnlyLowClus->Write(Form("hNCellVsEGammasNLTrueElecLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasSideBandOnlyHighClus->Write(Form("hNCellVsEGammasSideBandHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLSideBandOnlyHighClus->Write(Form("hNCellVsEGammasNLSideBandHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueSideBandOnlyHighClus->Write(Form("hNCellVsETrueGammasNLSideBandHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElecSideBandOnlyHighClus->Write(Form("hNCellVsEGammasNLTrueElecSideBandHigh_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsEGammasSideBandOnlyLowClus->Write(Form("hNCellVsEGammasSideBandLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLSideBandOnlyLowClus->Write(Form("hNCellVsEGammasNLSideBandLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueSideBandOnlyLowClus->Write(Form("hNCellVsETrueGammasNLSideBandLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEGammasNLTrueElecSideBandOnlyLowClus->Write(Form("hNCellVsEGammasNLTrueElecSideBandLow_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hHitmap->Write(Form("hHitmap_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hHitmapGammas->Write(Form("hHitmapGammas_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hHitmapElec->Write(Form("hHitmapElec_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hHitmapHadrons->Write(Form("hHitmapHadrons_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hClusE->Write(Form("hClusE_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hClusENL->Write(Form("hClusENL_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEelec->Write(Form("hNCellVsEelec_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEelecNL->Write(Form("hNCellVsEelecNL_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsEelecNLTrue->Write(Form("hNCellVsTrueEelecNL_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hEOverPElec->Write(Form("hEOverPElec_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hEOverPElecTrueElec->Write(Form("hEOverPElecTrueElec_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETMNL->Write(Form("hNCellVsETMNL_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsENoExotic50->Write(Form("hNCellVsENoExotic50_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsENoExotic75->Write(Form("hNCellVsENoExotic75_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETMNLS500A100->Write(Form("hNCellVsETMNLS500A100_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETMNLS500A105->Write(Form("hNCellVsETMNLS500A105_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETMNLS525A105->Write(Form("hNCellVsETMNLS525A105_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETMNLS550A110->Write(Form("hNCellVsETMNLS550A110_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETMNLS510A102->Write(Form("hNCellVsETMNLS510A102_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETMNLS500A110->Write(Form("hNCellVsETMNLS500A110_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsETrueGamma->Write(Form("hNCellVsETrueGamma_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETrueElec ->Write(Form("hNCellVsETrueElec_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hNCellVsETrueHadr ->Write(Form("hNCellVsETrueHadr_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hRecVsTrueE->Write(Form("hRecVsTrueE_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hRecDivTrueE->Write(Form("hRecDivTrueE_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hRecDivTrueEOneCell->Write(Form("hRecDivTrueEOneCell_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hRecDivTrueETwoCell->Write(Form("hRecDivTrueETwoCell_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  hRecDivTrueEThreeCell->Write(Form("hRecDivTrueEThreeCell_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);

  hNCellVsENL->Write(Form("hNCellVsENL_%s", (isMC == 1) ? "MC" : "data"), TObject::kOverwrite);
  fout.Close();

  gSystem->Exec(Form("cp MakeNCellEffiPureGammas.cxx %s_Code.cxx", filename.Data()));

}



//--------------------------------------
//--- Start the thing :)
//--------------------------------------
void MakeNCellEffiPureGammas(bool isMC = 1, TString period = "13TeVNomB", bool light = false, int skip = 0){
  TStopwatch watch;
  watch.Start();
  unsigned int nEvt = 1e7;
  if(isMC){
    DataTree myAna(true, period, light, skip);
    // myAna.SetDoNoTRD(1);
    // myAna.SetRejectBehindTRDSupport(true);
    myAna.RejectBorderCells(false);
    myAna.Process(nEvt*2.5);
    myAna.WriteToFile();
  } else {
    DataTree myAna2(false, period, light, skip);
    // myAna2.SetDoNoTRD(1);
    // myAna2.SetRejectBehindTRDSupport(true);
    myAna2.RejectBorderCells(false);
    myAna2.Process(nEvt);
    myAna2.WriteToFile();
  }
  watch.Stop();
  watch.Print();

}
