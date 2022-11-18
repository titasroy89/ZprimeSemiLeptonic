#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"

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
#include "UHH2/common/include/LuminosityHists.h"
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>

//#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include <UHH2/HOTVR/include/HadronicTop.h>
#include <UHH2/HOTVR/include/HOTVRScaleFactor.h>
#include <UHH2/HOTVR/include/HOTVRIds.h>
using namespace std;
using namespace uhh2;

#include "UHH2/common/include/NeuralNetworkBase.hpp"



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
uhh2::Event::Handle<float> h_Invpt_jet;
uhh2::Event::Handle<float> h_Invpt_jet1;
uhh2::Event::Handle<float> h_Invpt_jet2;
uhh2::Event::Handle<float> h_Invpt_jet3;
uhh2::Event::Handle<float> h_dR_mu_jet;
uhh2::Event::Handle<float> h_deepjetbscore_jet;
uhh2::Event::Handle<float> h_deepjetbscore_jet1;
uhh2::Event::Handle<float> h_deepjetbscore_jet2;
uhh2::Event::Handle<float> h_deepjetbscore_jet3;;
uhh2::Event::Handle<float> h_dphi_jet1_MET;
uhh2::Event::Handle<float> h_dphi_mu_MET ;
uhh2::Event::Handle<float> h_eta_jet;
uhh2::Event::Handle<float> h_eta_jet1;
uhh2::Event::Handle<float> h_eta_jet2;
uhh2::Event::Handle<float> h_eta_jet3;
uhh2::Event::Handle<float> h_eta_mu;
uhh2::Event::Handle<float> h_phi_jet;
uhh2::Event::Handle<float> h_phi_jet1;
uhh2::Event::Handle<float> h_phi_jet2;
uhh2::Event::Handle<float> h_phi_jet3;
uhh2::Event::Handle<float> h_phi_mu;
uhh2::Event::Handle<float> h_ptrel_mu_jet;
uhh2::Event::Handle<float> h_reliso_mu;
//uhh2::Event::Handle<float> h_weight;

};

NeuralNetworkModule::NeuralNetworkModule(Context& ctx, const std::string & ModelName, const std::string& ConfigName): NeuralNetworkBase(ctx, ModelName, ConfigName){
  bool debug = true;

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
  h_dphi_mu_MET = ctx.get_handle<float> ("dphi_mu_MET");
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
  //h_weight = ctx.get_handle<float> ("weight");
 
  if (debug) cout << "done declaring vars in NeuralNetworkModule" << endl;

}

void NeuralNetworkModule::CreateInputs(Event & event){
  bool debug=true;
  NNInputs.clear();
  NNoutputs.clear();
  if (debug) cout << "cleared NN inputs and outputs" << endl;
  string varname[34];
  string scal[34];
  string mean[34];
  string std[34];
  double mean_val[34];
  double std_val[34];
   
  if (debug) cout << "about to get info from norm.txt" << endl;
  ifstream normfile ("/nfs/dust/cms/user/titasroy/MLInputs_mu_UL_Nov17/StandardScaler/Merged__NOMINAL/classes_2_QCD_TTbarSemi_1+TTOther+TTbarSemi_2/NormInfo.txt", ios::in);
  if (normfile.is_open()){
        

        for(int i = 0; i < 34; ++i)
        {       
	   normfile >> varname[i] >> scal[i] >> mean[i] >> std[i];
           mean_val[i] = std::stod(mean[i]);
           std_val[i] = std::stod(std[i]);
 
        }   
        
    normfile.close();
  }
  if (debug) cout << "about to get tensor flow values" << endl;

  NNInputs.push_back( tensorflow::Tensor(tensorflow::DT_FLOAT, {1, 34}));
  
  vector<uhh2::Event::Handle<float>> inputs ={ h_InvMET,h_InvS11,h_InvS12,h_InvS13,h_InvS22,h_InvS23,h_InvS33,h_InvST,h_Invmass_jet,h_Invmass_jet1,h_Invmass_jet2,h_Invmass_jet3,h_Invpt_jet,h_Invpt_jet1,h_Invpt_jet2,h_Invpt_jet3,h_dR_mu_jet,h_deepjetbscore_jet,h_deepjetbscore_jet1,h_deepjetbscore_jet2,h_deepjetbscore_jet3,h_dphi_jet1_MET,h_dphi_mu_MET,h_eta_jet,h_eta_jet1,h_eta_jet2,h_eta_jet3,h_eta_mu,h_phi_jet,h_phi_jet1,h_phi_jet2,h_phi_jet3,h_phi_mu,h_ptrel_mu_jet,h_reliso_mu}; 
  

 for(int i = 0; i < 34; ++i){
    NNInputs.at(0).tensor<float, 2>()(0,i)  = (event.get(inputs.at(i))   - mean_val[i]) / (std_val[i]);
    if (debug) cout << "values:"<< (event.get(inputs.at(i))   - mean_val[i]) / (std_val[i]) << endl;
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
  std::unique_ptr<MuonCleaner>     muon_cleaner_low, muon_cleaner_high;
  std::unique_ptr<ElectronCleaner> electron_cleaner_low, electron_cleaner_high;
  // Scale Factors -- Systematics
  unique_ptr<AnalysisModule> sf_muon_iso_low, sf_muon_id_low, sf_muon_id_high, sf_muon_trigger_low, sf_muon_trigger_high;
  unique_ptr<AnalysisModule> sf_muon_iso_low_dummy, sf_muon_id_dummy, sf_muon_trigger_dummy;
  unique_ptr<AnalysisModule> sf_ele_id_low, sf_ele_id_high, sf_ele_reco;
  unique_ptr<AnalysisModule> sf_ele_id_dummy, sf_ele_reco_dummy;
  unique_ptr<MuonRecoSF> sf_muon_reco;
  unique_ptr<AnalysisModule> sf_btagging;
  //AnalysisModule
  unique_ptr<AnalysisModule> LumiWeight_module, PUWeight_module, TopPtReweight_module, MCScale_module;
  unique_ptr<AnalysisModule> Corrections_module;

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
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<Selection> sel_1btag, sel_2btag;
  std::unique_ptr<Selection> HEM_selection;
  unique_ptr<Selection> ThetaStar_selection_bin1, ThetaStar_selection_bin2, ThetaStar_selection_bin3;

  // NN variables handles
  //unique_ptr<Variables_NN> Variables_module;

  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;   Event::Handle<float> h_weight;

  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

  // Lumi hists
  std::unique_ptr<Hists> lumihists_Weights_Init, lumihists_Weights_PU, lumihists_Weights_Lumi, lumihists_Weights_TopPt, lumihists_Weights_MCScale, lumihists_Chi2;

  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

  // DNN multiclass output hist
  std::unique_ptr<Hists> h_MulticlassNN_output;




  
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

    Event::Handle<float> h_dphi_mu_MET ;
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
  debug = true;
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
  if (debug) cout <<"done def"<<endl;
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
  if (debug) cout <<"about to modules"<<endl;
  // Modules
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
  TopPtReweight_module.reset(new TopPtReweighting(ctx, a_toppt, b_toppt, Sys_TopPt_a, Sys_TopPt_b, ""));
  MCScale_module.reset(new MCScaleVariation(ctx));
  hadronic_top.reset(new HadronicTop(ctx));
  Corrections_module.reset(new NLOCorrections(ctx));
  if (debug) cout <<"about to btag"<<endl;
  // b-tagging SFs
  sf_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET, "CHS_matched"));

  // set lepton scale factors: see UHH2/common/include/LeptonScaleFactors.h
  //sf_muon_iso_low.reset(new uhh2::MuonIsoScaleFactors(ctx, Muon::Selector::PFIsoTight, Muon::Selector::CutBasedIdTight, true));
  //sf_muon_id_low.reset(new uhh2::MuonIdScaleFactors(ctx, Muon::Selector::CutBasedIdTight, true));
  //sf_muon_id_high.reset(new uhh2::MuonIdScaleFactors(ctx, Muon::Selector::CutBasedIdGlobalHighPt, true));
  //sf_muon_trigger_low.reset(new uhh2::MuonTriggerScaleFactors(ctx, false, true));
  //sf_muon_trigger_high.reset(new uhh2::MuonTriggerScaleFactors(ctx, true, false));
 // sf_muon_reco.reset(new MuonRecoSF(ctx));
 // sf_ele_id_low.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_iso_V2_wp80, true));
 // sf_ele_id_high.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::tag::mvaEleID_Fall17_noIso_V2_wp80, true));
 // sf_ele_reco.reset(new uhh2::ElectronRecoScaleFactors(ctx, false, true));
 // sf_muon_iso_low_dummy.reset(new uhh2::MuonIsoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
 // sf_muon_id_dummy.reset(new uhh2::MuonIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
 // sf_muon_trigger_dummy.reset(new uhh2::MuonTriggerScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));
 // sf_ele_id_dummy.reset(new uhh2::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
 // sf_ele_reco_dummy.reset(new uhh2::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
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
  h_Invpt_jet1 = ctx.get_handle<float> ("Invpt_jet");
  h_Invpt_jet1 = ctx.get_handle<float> ("Invpt_jet1");
  h_Invpt_jet2 = ctx.get_handle<float> ("Invpt_jet2");
  h_Invpt_jet3 = ctx.get_handle<float> ("Invpt_jet3");
  h_dR_mu_jet = ctx.get_handle<float> ("dR_mu_jet");
  h_deepjetbscore_jet = ctx.get_handle<float> ("deepjetbscore_jet");
  h_deepjetbscore_jet1 = ctx.get_handle<float> ("deepjetbscore_jet1");
  h_deepjetbscore_jet2 = ctx.get_handle<float> ("deepjetbscore_jet2");
  h_deepjetbscore_jet3 = ctx.get_handle<float> ("deepjetbscore_jet3");
  h_dphi_jet1_MET = ctx.get_handle<float> ("dphi_jet1_MET");
  h_dphi_mu_MET = ctx.get_handle<float> ("dphi_mu_MET");
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
   
  //h_weight = ctx.get_handle<float> ("weight");
 


  if (debug) cout << "about to get NN"<<endl;
  h_NNoutput = ctx.get_handle<std::vector<tensorflow::Tensor>>("NNoutput");
  h_NNoutput0 = ctx.declare_event_output<double>("NNoutput0");
  h_NNoutput1 = ctx.declare_event_output<double>("NNoutput1");
  if (debug) cout << "declared  NN"<<endl;
  NNModule.reset( new NeuralNetworkModule(ctx, "/nfs/dust/cms/user/titasroy/uhh2-106X/CMSSW_10_6_28/src/UHH2/model_Nov17_34.pb", "/nfs/dust/cms/user/titasroy/uhh2-106X/CMSSW_10_6_28/src/UHH2/model_Nov17_34.config.pbtxt"));
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
  }
  else if(isdeepAK8){
    TopTaggerDeepAK8->process(event);
  }

  fill_histograms(event, "Weights_Init");
  lumihists_Weights_Init->fill(event);

  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284); // calculated following instructions ar https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
  }
  //fill_histograms(event, "Weights_HEM");
  PUWeight_module->process(event);
  if(debug) cout << "PUWeight: ok" << endl;
  fill_histograms(event, "Weights_PU");
  lumihists_Weights_PU->fill(event);
  LumiWeight_module->process(event);
  if(debug) cout << "LumiWeight: ok" << endl;
  fill_histograms(event, "Weights_Lumi");
  lumihists_Weights_Lumi->fill(event);
  MCScale_module->process(event);
  if(debug) cout << "MCScale: ok" << endl;
  fill_histograms(event, "Weights_MCScale");
  lumihists_Weights_MCScale->fill(event);
  if (isMC) {
    if (Prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
    else if (Prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
    else if (Prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
  }
  fill_histograms(event, "Weights_Prefiring");
  ps_weights->process(event);
  if(debug) cout << "MCScale: ok" << endl;
  fill_histograms(event, "Weights_PS");
  lumihists_Weights_PS->fill(event);
  double muon_pt_high(55.);
  bool muon_is_low = false;
  bool muon_is_high = false;

  if(isMuon){
    vector<Muon>* muons = event.muons;
    for(unsigned int i=0; i<muons->size(); i++){
      if(event.muons->at(i).pt()<=muon_pt_high){
        muon_is_low = true;
      }
      else{
        muon_is_high = true;
      }
    }
  }
  sort_by_pt<Muon>(*event.muons);
  if(isMuon){
    if(!isEleTriggerMeasurement){ // For ele trigger SF do not veto additional electrons in muon channel
      if(!EleVeto_selection->passes(event)) return false;
    }
    if(muon_is_low){
      if(!NMuon1_selection->passes(event)) return false;
      muon_cleaner_low->process(event);
      if(!NMuon1_selection->passes(event)) return false;
      fill_histograms(event, "Muon1_LowPt");
    }
    if(muon_is_high){
      if(!NMuon1_selection->passes(event)) return false;
      muon_cleaner_high->process(event);
      if(!NMuon1_selection->passes(event)) return false;
      fill_histograms(event, "Muon1_HighPt");
    }
    if( !(muon_is_high || muon_is_low) ) return false;
    fill_histograms(event, "Muon1_Tot");
  }
  double electron_pt_high(120.);
  bool ele_is_low = false;
  bool ele_is_high = false;

  if(isElectron || (isMuon && isEleTriggerMeasurement)){
    vector<Electron>* electrons = event.electrons;
    for(unsigned int i=0; i<electrons->size(); i++){
      if(abs(event.electrons->at(i).eta()) > 1.44 && abs(event.electrons->at(i).eta()) < 1.57) return false; // remove gap electrons in transition region between the barrel and endcaps of ECAL
      if(event.electrons->at(i).pt()<=electron_pt_high){
        ele_is_low = true;
      }
      else{
        ele_is_high = true;
      }
    }
  }
  sort_by_pt<Electron>(*event.electrons);
  if(isMuon && isEleTriggerMeasurement){
    if(ele_is_low){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_low->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "1Mu1Ele_LowPt");
    }
    if(ele_is_high){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_high->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "1Mu1Ele_HighPt");
    }
    if( !(ele_is_high || ele_is_low) ) return false;
    fill_histograms(event, "1Mu1Ele_Tot");
  }

  // Select exactly 1 electron and 0 muons
  if(isElectron){
    if(!MuonVeto_selection->passes(event)) return false;
    if(ele_is_low){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_low->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "Ele1_LowPt");
    }
    if(ele_is_high){
      if(!NEle1_selection->passes(event)) return false;
      electron_cleaner_high->process(event);
      if(!NEle1_selection->passes(event)) return false;
      fill_histograms(event, "Ele1_HighPt");
    }
    if( !(ele_is_high || ele_is_low)) return false;
    fill_histograms(event, "Ele1_Tot");
  }


  // apply electron id scale factors
  if(isMuon && !isEleTriggerMeasurement){
    sf_ele_id_dummy->process(event);
  }
  if(isMuon && isEleTriggerMeasurement){
    if(ele_is_low){
      sf_ele_id_low->process(event);
    }
    else if(ele_is_high){
      sf_ele_id_high->process(event);
    }
  }
  if(isElectron){
    if(ele_is_low){
      sf_ele_id_low->process(event);
    }
    else if(ele_is_high){
      sf_ele_id_high->process(event);
    }
    fill_histograms(event, "IdEle_SF");
  }
  

  // apply muon isolation scale factors (low pT only)
  if(isMuon){
    if(muon_is_low){
      sf_muon_iso_low->process(event);
    }
    else if(muon_is_high){
      sf_muon_iso_low_dummy->process(event);
    }
    fill_histograms(event, "IsoMuon_SF");
  }
  if(isElectron){
    sf_muon_iso_low_dummy->process(event);
  }
  // apply muon id scale factors
  f(isMuon){
    if(muon_is_low){
      sf_muon_id_low->process(event);
    }
    else if(muon_is_high){
      sf_muon_id_high->process(event);
    }
    fill_histograms(event, "IdMuon_SF");
  }
  if(isElectron){
    sf_muon_id_dummy->process(event);
  }
  if(isMuon && !isEleTriggerMeasurement){
    sf_ele_reco_dummy->process(event);
  }
  if(isMuon && isEleTriggerMeasurement){
    sf_ele_reco->process(event);
  }
  if(isElectron){
    sf_ele_reco->process(event);
    fill_histograms(event, "RecoEle_SF");
  }
  sf_muon_reco->process(event);
  if(isMuon){
    if(muon_is_low){
      if(isUL16preVFP || isUL16postVFP){
        if(!(Trigger_mu_A_selection->passes(event) || Trigger_mu_B_selection->passes(event))) return false;
      }
      else{
        if(!Trigger_mu_A_selection->passes(event)) return false;
      }
    }
   if(muon_is_high){
      if(isUL16preVFP || isUL16postVFP){ // 2016
         if(!isMC){ //2016 DATA RunB
 		if( event.run < 274889){
            if(!Trigger_mu_C_selection->passes(event)) return false;
          }else{ // DATA above Run B
            if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_D_selection->passes(event))) return false;
          }
        }else{ // 2016 MC RunB
          float runB_UL16_mu = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
          if( runB_UL16_mu < 0.1429){
            if(!Trigger_mu_C_selection->passes(event)) return false;
          }else{  // 2016 MC above RunB
            if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_D_selection->passes(event))) return false;
          }
        }
      }
      if(isUL17 && !isMC){ //2017 DATA
        if(event.run <= 299329){ //RunB
         if(!Trigger_mu_C_selection->passes(event)) return false;
        }else{
          if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
        }
      }
      if(isUL17 && isMC){ // 2017 MC
        float runB_mu = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if(runB_mu <= 0.1158){
          if(!Trigger_mu_C_selection->passes(event)) return false;
        }else{
          if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
        }
      }
      if(isUL18){ //2018
        if(!(Trigger_mu_C_selection->passes(event) || Trigger_mu_E_selection->passes(event) || Trigger_mu_F_selection->passes(event))) return false;
      }
    }
    fill_histograms(event, "TriggerMuon");
    lumihists_TriggerMuon->fill(event);
  }
    if(isElectron){
	if(ele_is_high){
		if(isMC && (isUL16preVFP || isUL16postVFP || isUL18) ){
        if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;
      }
      if(isMC && isUL17){
        float runB_ele = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if(runB_ele <= 0.1158){ // in RunB (below runnumb 299329) Ele115 does not exist, use Ele35 instead. To apply randomly in MC if random numb < RunB percetage (11.58%, calculated by Christopher Matthies)
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;
        }else{
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event) || Trigger_ph_A_selection->passes(event))) return false;
        }
      }
      if(!isMC){
        if(isUL16preVFP || isUL16postVFP){
		if(isPhoton){
		  if(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)) return false;
                  if(!Trigger_ph_A_selection->passes(event)) return false;
                }else{ // electron stream
                  if(!(Trigger_ele_A_selection->passes(event) ||  Trigger_ele_B_selection->passes(event) )) return false;
          }
        }
        // 2017
        if(isUL17){
          if(event.run <= 299329){
            if(isPhoton){ // photon stream
              if(Trigger_ele_A_selection->passes(event)) return false;
              if(!Trigger_ph_A_selection->passes(event)) return false;
            }else{ // electron stream
              if(!Trigger_ele_A_selection->passes(event)) return false;
            }
          }else{ // above runnumb with Ele115
            if(isPhoton){ // photon stream
              if(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)) return false;
              if(!Trigger_ph_A_selection->passes(event)) return false;
            }else{ // electron stream
              if(! (Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event))) return false;
            }
          }
        }
        if(isUL18){
          if(!(Trigger_ele_A_selection->passes(event) || Trigger_ele_B_selection->passes(event)|| Trigger_ph_A_selection->passes(event))) return false;
        }

      }
    }
    fill_histograms(event, "TriggerEle");
    lumihists_TriggerEle->fill(event);
  }


   if(isMuon){
    if(muon_is_low){
      sf_muon_trigger_low->process(event);
    }
    if(muon_is_high){
      sf_muon_trigger_high->process(event);
    }
    fill_histograms(event, "TriggerMuon_SF");
  }
  if(isElectron){
    sf_muon_trigger_dummy->process(event);
  }

  if(!isEleTriggerMeasurement){
    if((event.muons->size()+event.electrons->size()) != 1) return false; //veto events without leptons or with too many
  }else{
    if( event.muons->size() != 1 && event.electrons->size() != 1) return false; //for ele trigger SF measurement: 1 ele and 1 mu
  }
  if(debug) cout<<"N leptons ok: Nelectrons="<<event.electrons->size()<<" Nmuons="<<event.muons->size()<<endl;

  if(isMuon && !isEleTriggerMeasurement && muon_is_high){
    if(!TwoDCut_selection->passes(event)) return false;
  }
  fill_histograms(event, "TwoDCut_Muon");
  lumihists_TwoDCut_Muon->fill(event);
  if(isElectron && ele_is_high){
    if(!TwoDCut_selection->passes(event)) return false;
  }
  fill_histograms(event, "TwoDCut_Ele");
  lumihists_TwoDCut_Ele->fill(event);

  if(isMuon && isEleTriggerMeasurement && (muon_is_high || ele_is_high)){
    if(!TwoDCut_selection->passes(event)) return false;
  }
  K4PuppiCHS_matching->process(event);
  h_CHSMatchHists->fill(event);

  if(!Jet1_selection->passes(event)) return false;
  if(debug) cout << "Jet1_selection: ok" << endl;
  fill_histograms(event, "Jet1");
  lumihists_Jet1->fill(event);

  if(!Jet2_selection->passes(event)) return false;
  if(debug) cout << "Jet2_selection: is ok" << endl;
  fill_histograms(event, "Jet2");
  lumihists_Jet2->fill(event);
  if(!met_sel->passes(event)) return false;
  if(debug) cout << "MET: ok" << endl;
  fill_histograms(event, "MET");
  lumihists_MET->fill(event);
  if(isMuon){
    if(!htlep_sel->passes(event)) return false;
    fill_histograms(event, "HTlep");
    lumihists_HTlep->fill(event);
    if(debug) cout << "HTlep: ok" << endl;
  }
  h_CHSMatchHists_beforeBTagSF->fill(event);
  fill_histograms(event, "BeforeBtagSF");
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
  h_CHSMatchHists_after2DBTagSF->fill(event);
  fill_histograms(event, "AfterCustomBtagSF");

  if(!AK4PuppiCHS_BTagging->passes(event)) return false;
  fill_histograms(event, "Btags1");
  h_CHSMatchHists_afterBTag->fill(event);

  NLOCorrections_module->process(event);
  if(!isEleTriggerMeasurement) sf_ele_trigger->process(event);
  CandidateBuilder->process(event);
  fill_histograms(event, "TTbarCandidate");
  if(debug) cout << "CandidateBuilder: ok" << endl;
  CorrectMatchDiscriminatorZprime->process(event);
  fill_histograms(event, "CorrectMatchDiscriminator");
  if(debug) cout << "CorrectMatchDiscriminatorZprime: ok" << endl;
  Chi2DiscriminatorZprime->process(event);
  fill_histograms(event, "Chi2Discriminator");
  if(debug) cout << "Chi2DiscriminatorZprime: ok" << endl;
  sort_by_pt<Jet>(*event.jets);
               



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
  if (debug) cout << "output 1"<< out0 << endl;
  if (debug) cout << "output 2" << out1 << endl;
  vector<double> out_event = {out0, out1};

  double max_score = 0.0;
  for ( int i = 0; i < 2; i++ ) {
    if ( out_event[i] > max_score) {
    max_score = out_event[i];
    }
  }

  if (debug) cout << "Fill NN vars" << endl;
  if (debug) cout << "max_score: "<< max_score<< endl;
  if (debug) cout << "out0: " <<out0<<endl;
  if (debug) cout << "out1: " <<out1<<endl;

  if( out0 == max_score ){
   if (debug) cout << "about to fill histograms for DNN_output0"<< endl;
  // fill_histograms(event, "DNN_output0");
  // if(is_aftercuts){
   //if (debug) cout <<"checking Chi2 in DNN"<<endl;
   if(!Chi2_selection->passes(event)) return false;
   fill_histograms(event, "DNN_output0");
   
   //if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output0");
  //  if(TopJetBtagSubjet_selection->passes(event) && ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output0");
    
   // }
 }

  if( out1 == max_score ){
    if (debug) cout << "about to fill histograms for DNN_output1"<< endl;
  //  fill_histograms(event, "DNN_output1");
  //  if(is_aftercuts){
    if(!Chi2_selection->passes(event)) return false;
    fill_histograms(event, "DNN_output1");
  //  if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output1");
    // if(TopJetBtagSubjet_selection->passes(event) && ZprimeTopTag_selection->passes(event)) fill_histograms(event, "DNN_output1");

    // } 
  }
  if (debug) cout << "done with filling DNN variables" << endl;
    
  


  return true;
 
}
UHH2_REGISTER_ANALYSIS_MODULE(ZprimeAnalysisModule_NNimplement)
