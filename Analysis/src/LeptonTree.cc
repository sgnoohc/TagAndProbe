#include "TagAndProbe/Analysis/interface/LeptonTree.h"
LeptonTree lepton_tree_obj;
namespace lepton_tree {
	const float &met() { return lepton_tree_obj.met(); }
	const float &metPhi() { return lepton_tree_obj.metPhi(); }
	const int &event() { return lepton_tree_obj.event(); }
	const int &lumi() { return lepton_tree_obj.lumi(); }
	const int &run() { return lepton_tree_obj.run(); }
	const bool &filt_csc() { return lepton_tree_obj.filt_csc(); }
	const bool &filt_hbhe() { return lepton_tree_obj.filt_hbhe(); }
	const bool &filt_hcallaser() { return lepton_tree_obj.filt_hcallaser(); }
	const bool &filt_ecaltp() { return lepton_tree_obj.filt_ecaltp(); }
	const bool &filt_trkfail() { return lepton_tree_obj.filt_trkfail(); }
	const bool &filt_eebadsc() { return lepton_tree_obj.filt_eebadsc(); }
	const bool &is_real_data() { return lepton_tree_obj.is_real_data(); }
	const float &scale1fb() { return lepton_tree_obj.scale1fb(); }
	const float &xsec() { return lepton_tree_obj.xsec(); }
	const float &kfactor() { return lepton_tree_obj.kfactor(); }
	const float &gen_met() { return lepton_tree_obj.gen_met(); }
	const float &gen_met_phi() { return lepton_tree_obj.gen_met_phi(); }
	const float &njets() { return lepton_tree_obj.njets(); }
	const float &ht() { return lepton_tree_obj.ht(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets() { return lepton_tree_obj.jets(); }
	const vector<float> &jets_disc() { return lepton_tree_obj.jets_disc(); }
	const TString &sample() { return lepton_tree_obj.sample(); }
	const int &nFOs() { return lepton_tree_obj.nFOs(); }
	const int &nvtx() { return lepton_tree_obj.nvtx(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &p4() { return lepton_tree_obj.p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &tag_p4() { return lepton_tree_obj.tag_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &dilep_p4() { return lepton_tree_obj.dilep_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mc_p4() { return lepton_tree_obj.mc_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mc_motherp4() { return lepton_tree_obj.mc_motherp4(); }
	const int &id() { return lepton_tree_obj.id(); }
	const int &idx() { return lepton_tree_obj.idx(); }
	const float &d0() { return lepton_tree_obj.d0(); }
	const float &dZ() { return lepton_tree_obj.dZ(); }
	const float &d0_err() { return lepton_tree_obj.d0_err(); }
	const int &motherID() { return lepton_tree_obj.motherID(); }
	const int &mc_id() { return lepton_tree_obj.mc_id(); }
	const float &iso() { return lepton_tree_obj.iso(); }
	const bool &passes_id() { return lepton_tree_obj.passes_id(); }
	const bool &passes_id_ptrel() { return lepton_tree_obj.passes_id_ptrel(); }
	const bool &passes_id_miniiso() { return lepton_tree_obj.passes_id_miniiso(); }
	const bool &passes_id_newminiiso() { return lepton_tree_obj.passes_id_newminiiso(); }
	const bool &FO() { return lepton_tree_obj.FO(); }
	const bool &FO_ptrel() { return lepton_tree_obj.FO_ptrel(); }
	const bool &FO_miniiso() { return lepton_tree_obj.FO_miniiso(); }
	const bool &FO_newminiiso() { return lepton_tree_obj.FO_newminiiso(); }
	const bool &FO_NoIso() { return lepton_tree_obj.FO_NoIso(); }
	const float &ip3d() { return lepton_tree_obj.ip3d(); }
	const float &ip3derr() { return lepton_tree_obj.ip3derr(); }
	const int &type() { return lepton_tree_obj.type(); }
	const float &mt() { return lepton_tree_obj.mt(); }
	const float &ptrelv0() { return lepton_tree_obj.ptrelv0(); }
	const float &ptrelv1() { return lepton_tree_obj.ptrelv1(); }
	const float &miniiso() { return lepton_tree_obj.miniiso(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &jet_close_lep() { return lepton_tree_obj.jet_close_lep(); }
	const int &tag_charge() { return lepton_tree_obj.tag_charge(); }
	const float &dilep_mass() { return lepton_tree_obj.dilep_mass(); }
	const float &el_sigmaIEtaIEta_full5x5() { return lepton_tree_obj.el_sigmaIEtaIEta_full5x5(); }
	const float &el_etaSC() { return lepton_tree_obj.el_etaSC(); }
	const float &el_dEtaIn() { return lepton_tree_obj.el_dEtaIn(); }
	const float &el_dPhiIn() { return lepton_tree_obj.el_dPhiIn(); }
	const float &el_hOverE() { return lepton_tree_obj.el_hOverE(); }
	const float &el_ecalEnergy() { return lepton_tree_obj.el_ecalEnergy(); }
	const float &el_eOverPIn() { return lepton_tree_obj.el_eOverPIn(); }
	const bool &el_conv_vtx_flag() { return lepton_tree_obj.el_conv_vtx_flag(); }
	const int &el_exp_innerlayers() { return lepton_tree_obj.el_exp_innerlayers(); }
	const int &el_charge() { return lepton_tree_obj.el_charge(); }
	const int &el_sccharge() { return lepton_tree_obj.el_sccharge(); }
	const int &el_ckf_charge() { return lepton_tree_obj.el_ckf_charge(); }
	const bool &el_threeChargeAgree() { return lepton_tree_obj.el_threeChargeAgree(); }
	const int &mu_pid_PFMuon() { return lepton_tree_obj.mu_pid_PFMuon(); }
	const float &mu_gfit_chi2() { return lepton_tree_obj.mu_gfit_chi2(); }
	const float &mu_gfit_ndof() { return lepton_tree_obj.mu_gfit_ndof(); }
	const int &mu_gfit_validSTAHits() { return lepton_tree_obj.mu_gfit_validSTAHits(); }
	const int &mu_numberOfMatchedStations() { return lepton_tree_obj.mu_numberOfMatchedStations(); }
	const int &mu_validPixelHits() { return lepton_tree_obj.mu_validPixelHits(); }
	const int &mu_nlayers() { return lepton_tree_obj.mu_nlayers(); }
}
