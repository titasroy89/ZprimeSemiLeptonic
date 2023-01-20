#include <iostream>
#include <memory>
#include <fstream>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/common/include/PrintingModules.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/LeptonScaleFactors.h>
#include <UHH2/common/include/PSWeights.h>

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSystematicsHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicPDFHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicMulticlassNNHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicCHSMatchHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>
#include <UHH2/ZprimeSemiLeptonic/include/ElecTriggerSF.h>

#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <UHH2/HOTVR/include/HadronicTop.h>
#include <UHH2/HOTVR/include/HOTVRScaleFactor.h>
#include <UHH2/HOTVR/include/HOTVRIds.h>

#include "UHH2/common/include/NeuralNetworkBase.hpp"

using namespace std;
using namespace uhh2;


class NeuralNetworkModule: public NeuralNetworkBase {







public:
  explicit NeuralNetworkModule(uhh2::Context&, const std::string & ModelName, const std::string& ConfigName);
  virtual void CreateInputs(uhh2::Event & event) override;
protected:

uhh2::Event::Handle<float> h_InvMET;
uhh2::Event::Handle<float> h_InvS11;
uhh2::Event::Handle<float> h_InvS12;
uhh2::Event::Handle<float> h_InvS13;
uhh2::Event::Handle<float> h_InvS22;
uhh2::Event::Handle<float> h_InvS23;
uhh2::Event::Handle<float> h_InvS33;
uhh2::Event::Handle<float> h_InvST;
uhh2::Event::Handle<float> h_Invmass_jet;
uhh2::Event::Handle<float> h_Invmass_jet1;
uhh2::Event::Handle<float> h_Invmass_jet2;
uhh2::Event::Handle<float> h_Invmass_jet3;
uhh2::Event::Handle<float> h_Invpt_ele;
uhh2::Event::Handle<float> h_Invpt_jet;
uhh2::Event::Handle<float> h_Invpt_jet1;
uhh2::Event::Handle<float> h_Invpt_jet2;
uhh2::Event::Handle<float> h_Invpt_jet3;
uhh2::Event::Handle<float> h_dR_ele_jet;
uhh2::Event::Handle<float> h_deepjetbscore_jet;
uhh2::Event::Handle<float> h_deepjetbscore_jet1;
uhh2::Event::Handle<float> h_deepjetbscore_jet2;
uhh2::Event::Handle<float> h_deepjetbscore_jet3;;
uhh2::Event::Handle<float> h_dphi_jet1_MET;

uhh2::Event::Handle<float> h_dphi_ele_MET ;
uhh2::Event::Handle<float> h_dphi_ele_jet1;

uhh2::Event::Handle<float> h_eta_jet;
uhh2::Event::Handle<float> h_eta_jet1;

uhh2::Event::Handle<float> h_eta_jet2;
uhh2::Event::Handle<float> h_eta_jet3;
uhh2::Event::Handle<float> h_eta_ele;
uhh2::Event::Handle<float> h_phi_jet;
uhh2::Event::Handle<float> h_phi_jet1;
uhh2::Event::Handle<float> h_phi_jet2;
uhh2::Event::Handle<float> h_phi_jet3;
uhh2::Event::Handle<float> h_phi_ele;
uhh2::Event::Handle<float> h_ptrel_ele_jet;
uhh2::Event::Handle<float> h_reliso_ele;
uhh2::Event::Handle<float> h_weight;

};

NeuralNetworkModule::NeuralNetworkModule(Context& ctx, const std::string & ModelName, const std::string& ConfigName): NeuralNetworkBase(ctx, ModelName, ConfigName){
  bool debug = true;
  bool debug2 = true;
  h_InvMET = ctx.get_handle<float> ("InvMET");
  h_InvS11 = ctx.get_handle<float> ("InvS11");
  h_InvS12 = ctx.get_handle<float> ("InvS12");
  h_InvS13 = ctx.get_handle<float> ("InvS13");
  h_InvS22 = ctx.get_handle<float> ("InvS22");
  h_InvS23 = ctx.get_handle<float> ("InvS23");
  h_InvS33 = ctx.get_handle<float> ("InvS33");
  h_InvST = ctx.get_handle<float> ("InvST");

  h_Invmass_jet = ctx.get_handle<float> ("Invmass_jet");
  h_Invmass_jet1 = ctx.get_handle<float> ("Invmass_jet1");
  h_Invmass_jet2 = ctx.get_handle<float> ("Invmass_jet2");
  h_Invmass_jet3 = ctx.get_handle<float> ("Invmass_jet3");
  h_Invpt_ele = ctx.get_handle<float> ("Invpt_ele");
  h_Invpt_jet = ctx.get_handle<float> ("Invpt_jet");
  h_Invpt_jet1 = ctx.get_handle<float> ("Invpt_jet2");
  h_Invpt_jet2 = ctx.get_handle<float> ("Invpt_jet2");
  h_Invpt_jet3 = ctx.get_handle<float> ("Invpt_jet3");  
  h_dR_ele_jet = ctx.get_handle<float> ("dR_ele_jet");
  h_deepjetbscore_jet = ctx.get_handle<float> ("deepjetbscore_jet");
  h_deepjetbscore_jet1 = ctx.get_handle<float> ("deepjetbscore_jet1");
  h_deepjetbscore_jet2 = ctx.get_handle<float> ("deepjetbscore_jet2");
  h_deepjetbscore_jet3 = ctx.get_handle<float> ("deepjetbscore_jet3");
  h_dphi_jet1_MET = ctx.get_handle<float> ("dphi_jet1_MET");
  h_dphi_ele_MET = ctx.get_handle<float> ("dphi_ele_MET");
  h_dphi_ele_jet1 = ctx.get_handle<float> ("dphi_ele_jet1");
  h_eta_jet = ctx.get_handle<float> ("eta_jet");
  h_eta_jet1 = ctx.get_handle<float> ("eta_jet1");
  h_eta_jet2 = ctx.get_handle<float> ("eta_jet2");
  h_eta_jet3 = ctx.get_handle<float> ("eta_jet3");
  h_eta_ele = ctx.get_handle<float> ("eta_ele");

  h_phi_jet = ctx.get_handle<float> ("phi_jet");
  h_phi_jet1 = ctx.get_handle<float> ("phi_jet1");
  h_phi_jet2 = ctx.get_handle<float> ("phi_jet2");
  h_phi_jet3 = ctx.get_handle<float> ("phi_jet3");
  h_phi_ele = ctx.get_handle<float> ("phi_ele");
  h_ptrel_ele_jet = ctx.get_handle<float> ("ptrel_ele_jet");
  h_reliso_ele = ctx.get_handle<float> ("reliso_ele");
  h_weight = ctx.get_handle<float> ("weight");
  //h_CHSjets_matched = ctx.get_handle<std::vector<Jet>>("CHS_matched");
  if (debug) cout << "done with setting variables" << endl;

}

void NeuralNetworkModule::CreateInputs(Event & event){
  bool debug=true;
  bool debug2=true;
  NNInputs.clear();
  NNoutputs.clear();
  if (debug) cout << "cleared NN inputs and outputs" << endl;
  string varname[37];
  string scal[37];
  string mean[37];
  string std[37];
  double mean_val[37];
  double std_val[37];
   
  if (debug) cout << "about to get info from norm.txt" << endl;
  //ifstream normfile ("/nfs/dust/cms/user/titasroy/MLInputs_mu_UL_Nov17/StandardScaler/Merged__NOMINAL/classes_2_QCD_TTbarSemi_1+TTOther+TTbarSemi_2/NormInfo.txt", ios::in);
   ifstream normfile ("/nfs/dust/cms/user/titasroy/MLInputs_ele_HLTinv_24June/MinMaxScaler/Merged__NOMINAL/classes_2_QCD_TTbar+ST+WJets+DYJets+Diboson/NormInfo.txt", ios::in);
  //ifstream normfile ("/nfs/dust/cms/user/titasroy/MLInputs_mu_HLTinv_24June/MinMaxScaler/Merged__NOMINAL/classes_2_QCD_TTbar+ST+WJets+DYJets+Diboson/NormInfo.txt", ios::in);
  if (normfile.is_open()){

        for(int i = 0; i < 37; ++i)
        {       
	   normfile >> varname[i] >> scal[i] >> mean[i] >> std[i];
           mean_val[i] = std::stod(mean[i]);
           std_val[i] = std::stod(std[i]);
 
        }   
        
    normfile.close();
  }
  if (debug) cout << "about to get tensor flow values" << endl;

  NNInputs.push_back( tensorflow::Tensor(tensorflow::DT_FLOAT, {1, 37}));
  
  vector<uhh2::Event::Handle<float>> inputs ={h_InvMET,h_InvS11,h_InvS12,h_InvS13,h_InvS22,h_InvS23,h_InvS33,h_InvST,h_Invmass_jet,h_Invmass_jet1,h_Invmass_jet2,h_Invmass_jet3,h_Invpt_ele,h_Invpt_jet,h_Invpt_jet1,h_Invpt_jet2,h_Invpt_jet3,h_dR_ele_jet,h_deepjetbscore_jet,h_deepjetbscore_jet1,h_deepjetbscore_jet2,h_deepjetbscore_jet3,h_dphi_ele_MET,h_dphi_ele_jet1,h_dphi_jet1_MET,h_eta_ele,h_eta_jet,h_eta_jet1,h_eta_jet2,h_eta_jet3,h_phi_jet,h_phi_jet1,h_phi_jet2,h_phi_jet3,h_phi_ele,h_ptrel_ele_jet,h_reliso_ele}; 
  

 for(int i = 0; i < 37; ++i){
    if (debug) cout << i<< endl;
     NNInputs.at(0).tensor<float, 2>()(0,i)  = (event.get(inputs.at(i))   - mean_val[i]) / (std_val[i]);
   // if (debug) cout << "values:"<< (event.get(inputs.at(i))   - mean_val[i]) / (std_val[i]) << endl;
  }


 if (debug) cout << "got all tensor float values" << endl;
 if (NNInputs.size()!=LayerInputs.size()) throw logic_error("NeuralNetworkModule.cxx: Create a number of inputs diffetent wrt. LayerInputs.size()="+to_string(LayerInputs.size()));
 if (debug) cout << "done with getting inputs"<<endl;
}
class ZprimeAnalysisModule_NNimplement : public ModuleBASE {

public:
  explicit ZprimeAnalysisModule_NNimplement(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:
  bool debug;
  bool debug2;

  // Cleaners
  std::unique_ptr<MuonCleaner>     muon_cleaner_low, muon_cleaner_high;
  std::unique_ptr<ElectronCleaner> electron_cleaner_low, electron_cleaner_high;
  // Scale Factors -- Systematics
  unique_ptr<AnalysisModule> sf_muon_iso_low, sf_muon_id_low, sf_muon_id_high, sf_muon_trigger_low, sf_muon_trigger_high;
  unique_ptr<AnalysisModule> sf_muon_iso_low_dummy, sf_muon_id_dummy, sf_muon_trigger_dummy;
  unique_ptr<AnalysisModule> sf_ele_id_low, sf_ele_id_high, sf_ele_reco;
  unique_ptr<AnalysisModule> sf_ele_id_dummy, sf_ele_reco_dummy;
  unique_ptr<MuonRecoSF> sf_muon_reco;
  unique_ptr<AnalysisModule> sf_ele_trigger;
  unique_ptr<AnalysisModule> sf_btagging;
  //AnalysisModule
  unique_ptr<AnalysisModule> LumiWeight_module, PUWeight_module, TopPtReweight_module, MCScale_module;
  unique_ptr<AnalysisModule> NLOCorrections_module;
  unique_ptr<PSWeights> ps_weights;

   // Top tagging
  unique_ptr<HOTVRTopTagger> TopTaggerHOTVR;
  unique_ptr<AnalysisModule> hadronic_top;
  unique_ptr<AnalysisModule> sf_toptag;
  unique_ptr<DeepAK8TopTagger> TopTaggerDeepAK8;

  // TopTags veto
  unique_ptr<Selection> TopTagVetoSelection;

  // Mass reconstruction
  unique_ptr<ZprimeCandidateBuilder> CandidateBuilder;

  // Chi2 discriminator
  unique_ptr<ZprimeChi2Discriminator> Chi2DiscriminatorZprime;
  unique_ptr<ZprimeCorrectMatchDiscriminator> CorrectMatchDiscriminatorZprime;

  // Selections
  unique_ptr<Selection> Chi2_selection, TTbarMatchable_selection, Chi2CandidateMatched_selection, ZprimeTopTag_selection;
  unique_ptr<Selection>Met_selection;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<Selection> sel_1btag, sel_2btag;
  std::unique_ptr<Selection> HEM_selection;
  unique_ptr<Selection> ThetaStar_selection_bin1, ThetaStar_selection_bin2, ThetaStar_selection_bin3;
  unique_ptr<Selection> SignSplit; 
  // NN variables handles
 // unique_ptr<Variables_NN> Variables_module;

  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;   
  Event::Handle<float> h_weight;

  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

  // Lumi hists
  std::unique_ptr<Hists> lumihists_Weights_Init, lumihists_Weights_PU, lumihists_Weights_Lumi, lumihists_Weights_TopPt, lumihists_Weights_MCScale, lumihists_Chi2;

  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

  // DNN multiclass output hist
  //std::unique_ptr<Hists> h_MulticlassNN_output;




  
 // Event::Handle<float> h_weight;
    Event::Handle<float> h_InvMET;
    Event::Handle<float> h_InvS11;
    Event::Handle<float> h_InvS12;
    Event::Handle<float> h_InvS13;
    Event::Handle<float> h_InvS22;
    Event::Handle<float> h_InvS23;
    Event::Handle<float> h_InvS33;
    Event::Handle<float> h_InvST;
    Event::Handle<float> h_Invmass_jet;
    Event::Handle<float> h_Invmass_jet1;
    Event::Handle<float> h_Invmass_jet2; 
    Event::Handle<float> h_Invmass_jet3;
    Event::Handle<float> h_Invpt_jet;
    Event::Handle<float> h_Invpt_jet1;
    Event::Handle<float> h_Invpt_jet2;
    Event::Handle<float> h_Invpt_jet3;
    Event::Handle<float> h_dR_mu_jet;
    Event::Handle<float> h_deepjetbscore_jet;
    Event::Handle<float> h_deepjetbscore_jet1;
    Event::Handle<float> h_deepjetbscore_jet2;
    Event::Handle<float> h_deepjetbscore_jet3;;
    Event::Handle<float> h_dphi_jet1_MET;

    Event::Handle<float> h_dphi_ele_MET ;
    Event::Handle<float> h_dphi_ele_jet1;
    Event::Handle<float> h_eta_jet;
    Event::Handle<float> h_eta_jet1;

    Event::Handle<float> h_eta_jet2;
    Event::Handle<float> h_eta_jet3;
    Event::Handle<float> h_eta_mu;
    Event::Handle<float> h_phi_jet;
    Event::Handle<float> h_phi_jet1;
    Event::Handle<float> h_phi_jet2;
    Event::Handle<float> h_phi_jet3;
    Event::Handle<float> h_phi_mu;
    Event::Handle<float> h_ptrel_mu_jet;
    Event::Handle<float> h_reliso_mu;

    Event::Handle<std::vector<tensorflow::Tensor> > h_NNoutput;
    Event::Handle<double> h_NNoutput0;
    Event::Handle<double> h_NNoutput1;

  std::unique_ptr<NeuralNetworkModule> NNModule;

   // Configuration
    bool isMC, ishotvr, isdeepAK8;
  string Sys_PU, Prefiring_direction, Sys_TopPt_a, Sys_TopPt_b;
  TString sample;
  int runnr_oldtriggers = 299368;

  bool isUL16preVFP, isUL16postVFP, isUL17, isUL18;
  bool isMuon, isElectron;
  bool isPhoton;
  TString year;

  TH2F *ratio_hist_muon;
  TH2F *ratio_hist_ele;
 // uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

};


void ZprimeAnalysisModule_NNimplement::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_Skimming";
    //book_HFolder(mytag, new TTbarLJHistsSkimming(ctx,mytag));
    mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicHists(ctx,mytag));
  }
}

void ZprimeAnalysisModule_NNimplement::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag + "_Skimming";
 // HFolder(mytag)->fill(event);
  mytag = tag+"_General";
  HFolder(mytag)->fill(event);
}


ZprimeAnalysisModule_NNimplement::ZprimeAnalysisModule_NNimplement(uhh2::Context& ctx){
  debug = false;
  debug2= false;
  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }
  if (debug) cout <<"in analysis module"<<endl;
  
   isMC = (ctx.get("dataset_type") == "MC");
  ishotvr = (ctx.get("is_hotvr") == "true");
  isdeepAK8 = (ctx.get("is_deepAK8") == "true");
  TString mode = "hotvr";
  if(isdeepAK8) mode = "deepAK8";
  string tmp = ctx.get("dataset_version");
  sample = tmp;
  isUL16preVFP  = (ctx.get("dataset_version").find("UL16preVFP")  != std::string::npos);
  isUL16postVFP = (ctx.get("dataset_version").find("UL16postVFP") != std::string::npos);
  isUL17        = (ctx.get("dataset_version").find("UL17")        != std::string::npos);
  isUL18        = (ctx.get("dataset_version").find("UL18")        != std::string::npos);
  if(isUL16preVFP) year = "UL16preVFP";
  if(isUL16postVFP) year = "UL16postVFP";
  if(isUL17) year = "UL17";
  if(isUL18) year = "UL18";

  isPhoton = (ctx.get("dataset_version").find("SinglePhoton") != std::string::npos);
  ElectronId eleID_low  = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80);
  ElectronId eleID_high = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80);
  MuonId     muID_low   = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId     muID_high  = MuonID(Muon::CutBasedIdGlobalHighPt);

  double electron_pt_low;
  if(isUL17){
    electron_pt_low = 38.;
  }
  else{
    electron_pt_low = 35.;
  }
  double muon_pt_low(30.);
  double electron_pt_high(120.);
  double muon_pt_high(55.);

  const MuonId muonID_low(AndId<Muon>(PtEtaCut(muon_pt_low, 2.4), muID_low));
  const ElectronId electronID_low(AndId<Electron>(PtEtaSCCut(electron_pt_low, 2.5), eleID_low));
  const MuonId muonID_high(AndId<Muon>(PtEtaCut(muon_pt_high, 2.4), muID_high));
  const ElectronId electronID_high(AndId<Electron>(PtEtaSCCut(electron_pt_high, 2.5), eleID_high));

  muon_cleaner_low.reset(new MuonCleaner(muonID_low));
  electron_cleaner_low.reset(new ElectronCleaner(electronID_low));
  muon_cleaner_high.reset(new MuonCleaner(muonID_high));
  electron_cleaner_high.reset(new ElectronCleaner(electronID_high));

  // Important selection values
  double chi2_max(30.);
  double MET_cut, HT_lep_cut;
  string trigger_mu_A, trigger_mu_B, trigger_mu_C, trigger_mu_D, trigger_mu_E, trigger_mu_F;
  string trigger_ele_A, trigger_ele_B;
  string trigger_ph_A;
  isMuon = false; isElectron = false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;

  if(isMuon){//semileptonic muon channel
    if(isUL17){
      trigger_mu_A = "HLT_IsoMu27_v*";
    }
    else{
      trigger_mu_A = "HLT_IsoMu24_v*";
    }
    trigger_mu_B = "HLT_IsoTkMu24_v*";
    trigger_mu_C = "HLT_Mu50_v*";
    trigger_mu_D = "HLT_TkMu50_v*";
    trigger_mu_E = "HLT_OldMu100_v*";
    trigger_mu_F = "HLT_TkMu100_v*";
    MET_cut = 50;
    HT_lep_cut = 150;

  }
  if(isElectron){
    trigger_ele_B = "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*";
    if(isUL16preVFP || isUL16postVFP){
      trigger_ele_A = "HLT_Ele27_WPTight_Gsf_v*";
    }
    if(isUL17){
      trigger_ele_A = "HLT_Ele35_WPTight_Gsf_v*";
    }
    if(isUL18){
      trigger_ele_A = "HLT_Ele32_WPTight_Gsf_v*";
    }
    if(isUL16preVFP || isUL16postVFP){
      trigger_ph_A = "HLT_Photon175_v*";
    }
    else{
      trigger_ph_A = "HLT_Photon200_v*";
    }
    MET_cut = 120;
    HT_lep_cut = 0;
  }


  const TopJetId toptagID = AndId<TopJet>(HOTVRTopTag(0.8, 140.0, 220.0, 50.0), Tau32Groomed(0.56));

  Sys_PU = ctx.get("Sys_PU");
  Prefiring_direction = ctx.get("Sys_prefiring");
  Sys_TopPt_a = ctx.get("Systematic_TopPt_a");
  Sys_TopPt_b = ctx.get("Systematic_TopPt_b");

  BTag::algo btag_algo = BTag::DEEPJET;
  BTag::wp btag_wp = BTag::WP_MEDIUM;
  JetId id_btag = BTag(btag_algo, btag_wp);

  double a_toppt = 0.0615;
  double b_toppt = -0.0005; // par b TopPt Reweighting

  // Modules
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
 // TopPtReweight_module.reset(new TopPtReweighting(ctx, a_toppt, b_toppt, Sys_TopPt_a, Sys_TopPt_b, ""));
  MCScale_module.reset(new MCScaleVariation(ctx));
  hadronic_top.reset(new HadronicTop(ctx));
  NLOCorrections_module.reset(new NLOCorrections(ctx));
  ps_weights.reset(new PSWeights(ctx));


  // b-tagging SFs
  sf_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET, "CHS_matched"));

  // set lepton scale factors: see UHH2/common/include/LeptonScaleFactors.h
  sf_muon_iso_low.reset(new uhh2::MuonIsoScaleFactors(ctx, Muon::Selector::PFIsoTight, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_low.reset(new uhh2::MuonIdScaleFactors(ctx, Muon::Selector::CutBasedIdTight, true));
  sf_muon_id_high.reset(new uhh2::MuonIdScaleFactors(ctx, Muon::Selector::CutBasedIdGlobalHighPt, true));
  sf_muon_trigger_low.reset(new uhh2::MuonTriggerScaleFactors(ctx, false, true));
  sf_muon_trigger_high.reset(new uhh2::MuonTriggerScaleFactors(ctx, true, false));
  sf_muon_reco.reset(new MuonRecoSF(ctx));
  sf_ele_id_low.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_iso_V2_wp80, true));
  sf_ele_id_high.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_noIso_V2_wp80, true));
  sf_ele_reco.reset(new uhh2::ElectronRecoScaleFactors(ctx, false, true));

  sf_ele_trigger.reset( new uhh2::ElecTriggerSF(ctx, "central", "eta_ptbins", year) );
  // dummies (needed to aviod set value errors)
  sf_muon_iso_low_dummy.reset(new uhh2::MuonIsoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_id_dummy.reset(new uhh2::MuonIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_muon_trigger_dummy.reset(new uhh2::MuonTriggerScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
  sf_ele_id_dummy.reset(new uhh2::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_ele_reco_dummy.reset(new uhh2::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));


  met_sel.reset(new METCut  (MET_cut   , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(HT_lep_cut, uhh2::infinity));

  Chi2_selection.reset(new Chi2Cut(ctx, 0., chi2_max));
  TTbarMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  Chi2CandidateMatched_selection.reset(new Chi2CandidateMatchedSelection(ctx));
  ZprimeTopTag_selection.reset(new ZprimeTopTagSelection(ctx));

  HEM_selection.reset(new HEMSelection(ctx)); // HEM issue in 2018, veto on leptons and jets

  //Variables_module.reset(new Variables_NN(ctx, mode)); // variables for NN

  // Selections on scattering angle theta star
  double theta_bin1(0.5);
  ThetaStar_selection_bin1.reset(new ThetaStarSelection(ctx, theta_bin1));
  double theta_bin2(0.7);
  ThetaStar_selection_bin2.reset(new ThetaStarSelection(ctx, theta_bin2));
  double theta_bin3(0.9);
  ThetaStar_selection_bin3.reset(new ThetaStarSelection(ctx, theta_bin3));

  // Top Taggers
  TopTaggerHOTVR.reset(new HOTVRTopTagger(ctx));
  TopTaggerDeepAK8.reset(new DeepAK8TopTagger(ctx));

  // TopTags veto
  TopTagVetoSelection.reset(new TopTag_VetoSelection(ctx, mode));

  // Zprime candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // Zprime discriminators
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx));
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx));
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");

  h_weight = ctx.declare_event_output<float> ("weight");

  sel_1btag.reset(new NJetSelection(1, -1, id_btag));
  sel_2btag.reset(new NJetSelection(2,-1, id_btag));
  
  vector<string> histogram_tags = {"DNN_output0","DNN_output1"};
  book_histograms(ctx, histogram_tags); 
  h_InvMET = ctx.get_handle<float> ("InvMET");
  h_InvS11 = ctx.get_handle<float> ("InvS11");
  h_InvS12 = ctx.get_handle<float> ("InvS12");
  h_InvS13 = ctx.get_handle<float> ("InvS13");
  h_InvS22 = ctx.get_handle<float> ("InvS22");
  h_InvS23 = ctx.get_handle<float> ("InvS23");
  h_InvS33 = ctx.get_handle<float> ("InvS33");
  h_InvST = ctx.get_handle<float> ("InvST");

  h_Invmass_jet = ctx.get_handle<float> ("Invmass_jet");
  h_Invmass_jet1 = ctx.get_handle<float> ("Invmass_jet1");
  h_Invmass_jet2 = ctx.get_handle<float> ("Invmass_jet2");
  h_Invmass_jet3 = ctx.get_handle<float> ("Invmass_jet3");
  h_Invpt_jet = ctx.get_handle<float> ("Invpt_jet");
  h_Invpt_jet1 = ctx.get_handle<float> ("Invpt_jet1");
  h_Invpt_jet2 = ctx.get_handle<float> ("Invpt_jet2");
  h_Invpt_jet3 = ctx.get_handle<float> ("Invpt_jet3");
  h_dR_mu_jet = ctx.get_handle<float> ("dR_mu_jet");
  h_deepjetbscore_jet = ctx.get_handle<float> ("deepjetbscore_jet");
  h_deepjetbscore_jet1 = ctx.get_handle<float> ("deepjetbscore_jet1");
  h_deepjetbscore_jet2 = ctx.get_handle<float> ("deepjetbscore_jet2");
  h_deepjetbscore_jet3 = ctx.get_handle<float> ("deepjetbscore_jet3");
  h_dphi_jet1_MET = ctx.get_handle<float> ("dphi_jet1_MET");
  h_dphi_ele_MET = ctx.get_handle<float> ("dphi_ele_MET");
  h_dphi_ele_jet1 = ctx.get_handle<float> ("dphi_ele_jet1");
  h_eta_jet = ctx.get_handle<float> ("eta_jet");
  h_eta_jet1 = ctx.get_handle<float> ("eta_jet1");
  h_eta_jet2 = ctx.get_handle<float> ("eta_jet2");
  h_eta_jet3 = ctx.get_handle<float> ("eta_jet3");
  h_eta_mu = ctx.get_handle<float> ("eta_mu");

  h_phi_jet = ctx.get_handle<float> ("phi_jet");
  h_phi_jet1 = ctx.get_handle<float> ("phi_jet1");
  h_phi_jet2 = ctx.get_handle<float> ("phi_jet2");
  h_phi_jet3 = ctx.get_handle<float> ("phi_jet3");
  h_phi_mu = ctx.get_handle<float> ("phi_mu");
  h_ptrel_mu_jet = ctx.get_handle<float> ("ptrel_mu_jet");
  h_reliso_mu = ctx.get_handle<float> ("reliso_mu");
   
  h_weight = ctx.get_handle<float> ("weight");
  if(isMC){
    TString sample_name = "";
    vector<string> names  = {"ST", "WJets", "DY", "QCD","ALP_ttbar_signal", "ALP_ttbar_interference", "HscalarToTTTo1L1Nu2J_m365_w36p5_res", "HscalarToTTTo1L1Nu2J_m400_w40p0_res", "HscalarToTTTo1L1Nu2J_m500_w50p0_res", "HscalarToTTTo1L1Nu2J_m600_w60p0_res", "HscalarToTTTo1L1Nu2J_m800_w80p0_res", "HscalarToTTTo1L1Nu2J_m1000_w100p0_res", "HscalarToTTTo1L1Nu2J_m365_w36p5_int", "HscalarToTTTo1L1Nu2J_m400_w40p0_int", "HscalarToTTTo1L1Nu2J_m500_w50p0_int", "HscalarToTTTo1L1Nu2J_m600_w60p0_int", "HscalarToTTTo1L1Nu2J_m800_w80p0_int", "HscalarToTTTo1L1Nu2J_m1000_w100p0_int", "HpseudoToTTTo1L1Nu2J_m365_w36p5_res", "HpseudoToTTTo1L1Nu2J_m400_w40p0_res", "HpseudoToTTTo1L1Nu2J_m500_w50p0_res", "HpseudoToTTTo1L1Nu2J_m600_w60p0_res", "HpseudoToTTTo1L1Nu2J_m800_w80p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_res", "HpseudoToTTTo1L1Nu2J_m365_w36p5_int", "HpseudoToTTTo1L1Nu2J_m400_w40p0_int", "HpseudoToTTTo1L1Nu2J_m500_w50p0_int", "HpseudoToTTTo1L1Nu2J_m600_w60p0_int", "HpseudoToTTTo1L1Nu2J_m800_w80p0_int", "HpseudoToTTTo1L1Nu2J_m1000_w100p0_int", "HscalarToTTTo1L1Nu2J_m365_w91p25_res", "HscalarToTTTo1L1Nu2J_m400_w100p0_res", "HscalarToTTTo1L1Nu2J_m500_w125p0_res", "HscalarToTTTo1L1Nu2J_m600_w150p0_res", "HscalarToTTTo1L1Nu2J_m800_w200p0_res", "HscalarToTTTo1L1Nu2J_m1000_w250p0_res", "HscalarToTTTo1L1Nu2J_m365_w91p25_int", "HscalarToTTTo1L1Nu2J_m400_w100p0_int", "HscalarToTTTo1L1Nu2J_m500_w125p0_int", "HscalarToTTTo1L1Nu2J_m600_w150p0_int", "HscalarToTTTo1L1Nu2J_m800_w200p0_int", "HscalarToTTTo1L1Nu2J_m1000_w250p0_int", "HpseudoToTTTo1L1Nu2J_m365_w91p25_res", "HpseudoToTTTo1L1Nu2J_m400_w100p0_res", "HpseudoToTTTo1L1Nu2J_m500_w125p0_res", "HpseudoToTTTo1L1Nu2J_m600_w150p0_res", "HpseudoToTTTo1L1Nu2J_m800_w200p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_res", "HpseudoToTTTo1L1Nu2J_m365_w91p25_int", "HpseudoToTTTo1L1Nu2J_m400_w100p0_int", "HpseudoToTTTo1L1Nu2J_m500_w125p0_int", "HpseudoToTTTo1L1Nu2J_m600_w150p0_int", "HpseudoToTTTo1L1Nu2J_m800_w200p0_int", "HpseudoToTTTo1L1Nu2J_m1000_w250p0_int", "HscalarToTTTo1L1Nu2J_m365_w9p125_res", "HscalarToTTTo1L1Nu2J_m400_w10p0_res", "HscalarToTTTo1L1Nu2J_m500_w12p5_res", "HscalarToTTTo1L1Nu2J_m600_w15p0_res", "HscalarToTTTo1L1Nu2J_m800_w20p0_res", "HscalarToTTTo1L1Nu2J_m1000_w25p0_res", "HscalarToTTTo1L1Nu2J_m365_w9p125_int", "HscalarToTTTo1L1Nu2J_m400_w10p0_int", "HscalarToTTTo1L1Nu2J_m500_w12p5_int", "HscalarToTTTo1L1Nu2J_m600_w15p0_int", "HscalarToTTTo1L1Nu2J_m800_w20p0_int", "HscalarToTTTo1L1Nu2J_m1000_w25p0_int", "HpseudoToTTTo1L1Nu2J_m365_w9p125_res", "HpseudoToTTTo1L1Nu2J_m400_w10p0_res", "HpseudoToTTTo1L1Nu2J_m500_w12p5_res", "HpseudoToTTTo1L1Nu2J_m600_w15p0_res", "HpseudoToTTTo1L1Nu2J_m800_w20p0_res", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_res", "HpseudoToTTTo1L1Nu2J_m365_w9p125_int", "HpseudoToTTTo1L1Nu2J_m400_w10p0_int", "HpseudoToTTTo1L1Nu2J_m500_w12p5_int", "HpseudoToTTTo1L1Nu2J_m600_w15p0_int", "HpseudoToTTTo1L1Nu2J_m800_w20p0_int", "HpseudoToTTTo1L1Nu2J_m1000_w25p0_int", "RSGluonToTT_M-500", "RSGluonToTT_M-1000", "RSGluonToTT_M-1500", "RSGluonToTT_M-2000", "RSGluonToTT_M-2500", "RSGluonToTT_M-3000", "RSGluonToTT_M-3500", "RSGluonToTT_M-4000", "RSGluonToTT_M-4500", "RSGluonToTT_M-5000", "RSGluonToTT_M-5500", "RSGluonToTT_M-6000", "ZPrimeToTT_M400_W40", "ZPrimeToTT_M500_W50", "ZPrimeToTT_M600_W60", "ZPrimeToTT_M700_W70", "ZPrimeToTT_M800_W80", "ZPrimeToTT_M900_W90", "ZPrimeToTT_M1000_W100", "ZPrimeToTT_M1200_W120", "ZPrimeToTT_M1400_W140", "ZPrimeToTT_M1600_W160", "ZPrimeToTT_M1800_W180", "ZPrimeToTT_M2000_W200", "ZPrimeToTT_M2500_W250", "ZPrimeToTT_M3000_W300", "ZPrimeToTT_M3500_W350", "ZPrimeToTT_M4000_W400", "ZPrimeToTT_M4500_W450", "ZPrimeToTT_M5000_W500", "ZPrimeToTT_M6000_W600", "ZPrimeToTT_M7000_W700", "ZPrimeToTT_M8000_W800", "ZPrimeToTT_M9000_W900", "ZPrimeToTT_M400_W120", "ZPrimeToTT_M500_W150", "ZPrimeToTT_M600_W180", "ZPrimeToTT_M700_W210", "ZPrimeToTT_M800_W240", "ZPrimeToTT_M900_W270", "ZPrimeToTT_M1000_W300", "ZPrimeToTT_M1200_W360", "ZPrimeToTT_M1400_W420", "ZPrimeToTT_M1600_W480", "ZPrimeToTT_M1800_W540", "ZPrimeToTT_M2000_W600", "ZPrimeToTT_M2500_W750", "ZPrimeToTT_M3000_W900", "ZPrimeToTT_M3500_W1050", "ZPrimeToTT_M4000_W1200", "ZPrimeToTT_M4500_W1350", "ZPrimeToTT_M5000_W1500", "ZPrimeToTT_M6000_W1800", "ZPrimeToTT_M7000_W2100", "ZPrimeToTT_M8000_W2400", "ZPrimeToTT_M9000_W2700", "ZPrimeToTT_M400_W4", "ZPrimeToTT_M500_W5", "ZPrimeToTT_M600_W6", "ZPrimeToTT_M700_W7", "ZPrimeToTT_M800_W8", "ZPrimeToTT_M900_W9", "ZPrimeToTT_M1000_W10", "ZPrimeToTT_M1200_W12", "ZPrimeToTT_M1400_W14", "ZPrimeToTT_M1600_W16", "ZPrimeToTT_M1800_W18", "ZPrimeToTT_M2000_W20", "ZPrimeToTT_M2500_W25", "ZPrimeToTT_M3000_W30", "ZPrimeToTT_M3500_W35", "ZPrimeToTT_M4000_W40", "ZPrimeToTT_M4500_W45", "ZPrimeToTT_M5000_W50", "ZPrimeToTT_M6000_W60", "ZPrimeToTT_M7000_W70", "ZPrimeToTT_M8000_W80", "ZPrimeToTT_M9000_W90"};

    for(unsigned int i=0; i<names.size(); i++){
      if( ctx.get("dataset_version").find(names.at(i)) != std::string::npos ) sample_name = names.at(i);
    }
    if( (ctx.get("dataset_version").find("TTToHadronic") != std::string::npos) || (ctx.get("dataset_version").find("TTToSemiLeptonic") != std::string::npos) || (ctx.get("dataset_version").find("TTTo2L2Nu") != std::string::npos) ) sample_name = "TTbar";
    if( (ctx.get("dataset_version").find("WW") != std::string::npos) || (ctx.get("dataset_version").find("ZZ") != std::string::npos) || (ctx.get("dataset_version").find("WZ") != std::string::npos) ) sample_name = "Diboson";

    if(isMuon){
      TFile* f_btag2Dsf_muon = new TFile("/nfs/dust/cms/user/deleokse/RunII_106_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/data/customBtagSF_muon_"+year+".root");
      ratio_hist_muon = (TH2F*)f_btag2Dsf_muon->Get("N_Jets_vs_HT_" + sample_name);
      ratio_hist_muon->SetDirectory(0);
    }
    else if(!isMuon){
      TFile* f_btag2Dsf_ele = new TFile("/nfs/dust/cms/user/deleokse/RunII_106_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/data/customBtagSF_electron_"+year+".root");
      ratio_hist_ele = (TH2F*)f_btag2Dsf_ele->Get("N_Jets_vs_HT_" + sample_name);
      ratio_hist_ele->SetDirectory(0);
    }
  }
  //h_MulticlassNN_output.reset(new ZprimeSemiLeptonicMulticlassNNHists(ctx, "MulticlassNN"));

  if (debug) cout << "about to get NN"<<endl;
  h_NNoutput = ctx.get_handle<std::vector<tensorflow::Tensor>>("NNoutput");
  h_NNoutput0 = ctx.declare_event_output<double>("NNoutput0");
  h_NNoutput1 = ctx.declare_event_output<double>("NNoutput1");
  if (debug) cout << "declared  NN"<<endl;
  //NNModule.reset( new NeuralNetworkModule(ctx, "/nfs/dust/cms/user/titasroy/uhh2-106X/CMSSW_10_6_28/src/UHH2/model_Nov30.pb", "/nfs/dust/cms/user/titasroy/uhh2-106X/CMSSW_10_6_28/src/UHH2/model_Nov30.config.pbtxt"));
  NNModule.reset( new NeuralNetworkModule(ctx, "/nfs/dust/cms/user/titasroy/uhh2-106X/CMSSW_10_6_28/src/UHH2/model_ele_July21.pb", "/nfs/dust/cms/user/titasroy/uhh2-106X/CMSSW_10_6_28/src/UHH2/model_ele_July21.config.pbtxt"));
  if (debug) cout << "loaded  NN model"<<endl;

}
  

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimeAnalysisModule_NNimplement::process(uhh2::Event& event){

  event.set(h_NNoutput0, 0);
  event.set(h_NNoutput1, 0);
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  event.set(h_chi2,-100);
  event.set(h_weight,-100);
  if(ishotvr){
    TopTaggerHOTVR->process(event);
    hadronic_top->process(event);
  }else if(isdeepAK8){
    TopTaggerDeepAK8->process(event);
  }
  if(debug) cout<<"Top Tagger ok"<<endl;


  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284);
    }
  LumiWeight_module->process(event);
  if(debug)  cout<<"LumiWeight ok"<<endl;
  PUWeight_module->process(event);
  //TopPtReweight_module->process(event);
  MCScale_module->process(event);
  if (isMC) {
     if (Prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
     else if (Prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
     else if (Prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
  }
  ps_weights->process(event);
  

  double muon_pt_high(55.);
  bool muon_is_low = false;
  bool muon_is_high = false;

  if(isMuon){
    vector<Muon>* muons = event.muons;
    for(unsigned int i=0; i<muons->size(); i++){
      if(event.muons->at(i).pt()<=muon_pt_high){
        muon_is_low = false;
      }else{
        muon_is_high = true;
      }
    }
  }
  sort_by_pt<Muon>(*event.muons);
  double electron_pt_high(120.);
  bool ele_is_low = false;
  bool ele_is_high = false;

  if(isElectron){
    vector<Electron>* electrons = event.electrons;
    for(unsigned int i=0; i<electrons->size(); i++){
      if(event.electrons->at(i).pt()<=electron_pt_high){
        ele_is_low = false;
      }else{
        ele_is_high = true;
      }
    }
  }
  sort_by_pt<Electron>(*event.electrons);
  if(isMuon){
    sf_ele_id_dummy->process(event);
   }
  if(isElectron){
    if(ele_is_low){
      sf_ele_id_low->process(event);
    }
    else if(ele_is_high){
      sf_ele_id_high->process(event);
    }
   // fill_histograms(event, "IdEle_SF");
   }
   if(isMuon){
	if(muon_is_high){
		muon_cleaner_high->process(event);
	}
   }
   if(isElectron){
	if(ele_is_high){
		electron_cleaner_high->process(event);
	}
   }
   if(isMuon){
     if(muon_is_low){
      sf_muon_iso_low->process(event);
     }
      else if(muon_is_high){
    	  sf_muon_iso_low_dummy->process(event);
     }
   // fill_histograms(event, "IsoMuon_SF");
   }
  if(isElectron){
    sf_muon_iso_low_dummy->process(event);
   }
  if(isMuon){
    if(muon_is_low){
      sf_muon_id_low->process(event);
    }
    else if(muon_is_high){
      sf_muon_id_high->process(event);
    }
   // fill_histograms(event, "IdMuon_SF");
  }
  if(isElectron){
    sf_muon_id_dummy->process(event);
  }

  if(isMuon){
    sf_ele_reco_dummy->process(event);
   }
  if(isElectron){
    sf_ele_reco->process(event);
   // fill_histograms(event, "RecoEle_SF");
  }
   sf_muon_reco->process(event);
  //fill_histograms(event, "MuonReco_SF");

  if(isMuon){
    if(muon_is_low){
      sf_muon_trigger_low->process(event);
    }
    if(muon_is_high){
      sf_muon_trigger_high->process(event);
    }
    //fill_histograms(event, "TriggerMuon_SF");
   }
  if(isElectron){
     sf_muon_trigger_dummy->process(event);
  }
  if(!met_sel->passes(event)) return false;
  if(isMuon){
    if(!htlep_sel->passes(event)) return false;
  }
  sel_1btag->passes(event);
  sf_btagging->process(event);
  if(isMC && isMuon){
     float custom_sf;

     vector<Jet>* jets = event.jets;
     int Njets = jets->size();
     double st_jets = 0.;
     for(const auto & jet : *jets) st_jets += jet.pt();
     custom_sf = ratio_hist_muon->GetBinContent( ratio_hist_muon->GetXaxis()->FindBin(Njets), ratio_hist_muon->GetYaxis()->FindBin(st_jets) );

     event.weight *= custom_sf;
  }
  if(isMC && !isMuon){
     float custom_sf;

     vector<Jet>* jets = event.jets;
     int Njets = jets->size();
     double st_jets = 0.;
     for(const auto & jet : *jets) st_jets += jet.pt();
     custom_sf = ratio_hist_ele->GetBinContent( ratio_hist_ele->GetXaxis()->FindBin(Njets), ratio_hist_ele->GetYaxis()->FindBin(st_jets) );

     event.weight *= custom_sf;
  }
  NLOCorrections_module->process(event);
  sf_ele_trigger->process(event);
  CandidateBuilder->process(event);
  Chi2DiscriminatorZprime->process(event);
  CorrectMatchDiscriminatorZprime->process(event);
  //Variables_module->process(event);
  //fill_histograms(event, "NNInputsBeforeReweight");
    



   //Neural network outputs
  if (debug) cout << "starting NN getoutputs" << endl;
  NNModule->process(event);
  if (debug) cout << "done with NNmodule processing events" <<endl;
  std::vector<tensorflow::Tensor> NNoutputs = NNModule->GetOutputs();
  //if (debug) cout << (double)(NNoutputs[0].tensor<float, 2>()(0,0)), (double)(NNoutputs[0].tensor<float, 2>()(0,1)) <<endl;
  event.set(h_NNoutput0, (double)(NNoutputs[0].tensor<float, 2>()(0,0)));
  event.set(h_NNoutput1, (double)(NNoutputs[0].tensor<float, 2>()(0,1)));
  event.set(h_NNoutput, NNoutputs);

  double out0 = (double)(NNoutputs[0].tensor<float, 2>()(0,0));
  double out1 = (double)(NNoutputs[0].tensor<float, 2>()(0,1));
  if (debug2) cout << "output 1"<< out0 << endl;
  if (debug2) cout << "output 2" << out1 << endl;
  vector<double> out_event = {out0, out1};
 // if (out1> out0){
//	if (debug2) cout<< "choosing ttbar"<<endl;
 //  }
  double max_score = 0.0;
  for ( int i = 0; i < 2; i++ ) {
    if ( out_event[i] > max_score) {
    	max_score = out_event[i];
    }
  }

  if (debug2) cout << "Fill NN vars" << endl;
  if (debug2) cout << "max_score: "<< max_score<< endl;
  if (debug2) cout << "out0: " <<out0<<endl;
  if (debug2) cout << "out1: " <<out1<<endl;
  
  if( out0 == max_score ){
      if (debug2) cout << "about to fill histograms for QCD"<< endl;
  // fill_histograms(event, "DNN_output0");
  // if(is_aftercuts){
   //if (debug) cout <<"checking Chi2 in DNN"<<endl;
   	if(!Chi2_selection->passes(event)) return false;
  	 fill_histograms(event, "DNN_output0");
	}
  //else{
//	if (debug2) cout << "about to fill histograms for TTbar"<< endl;
//	if(!Chi2_selection->passes(event)) return false;
//	fill_histograms(event, "DNN_output1");
 //  }
   
   //if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output0");
  //  if(TopJetBtagSubjet_selection->passes(event) && ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output0");
    
   // }
// }
  //if (debug2) cout << "about to check score"<< endl;
  if( out1 == max_score ){
    if (debug2) cout << "about to fill histograms for DNN_output1"<< endl;
  //  fill_histograms(event, "DNN_output1");
  //  if(is_aftercuts){
     if(!Chi2_selection->passes(event)) return false;
  //  if (debug2) cout << "about to fill 2nd node"<<endl;
     fill_histograms(event, "DNN_output1");
  //  if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output1");
    // if(TopJetBtagSubjet_selection->passes(event) && ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output1");

     } 
 // }
  if (debug) cout << "done with filling DNN variables" << endl;
    
  


  return true;
 
}
UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_NNimplement)
