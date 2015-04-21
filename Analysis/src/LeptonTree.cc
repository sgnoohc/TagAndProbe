#include "TagAndProbe/Analysis/interface/LeptonTree.h"
LeptonTree lepton_tree_obj;
namespace lepton_tree {
	const float &evt_pfmet() { return lepton_tree_obj.evt_pfmet(); }
	const float &evt_pfmetPhi() { return lepton_tree_obj.evt_pfmetPhi(); }
	const int &evt_event() { return lepton_tree_obj.evt_event(); }
	const int &evt_lumiBlock() { return lepton_tree_obj.evt_lumiBlock(); }
	const int &evt_run() { return lepton_tree_obj.evt_run(); }
	const bool &filt_csc() { return lepton_tree_obj.filt_csc(); }
	const bool &filt_hbhe() { return lepton_tree_obj.filt_hbhe(); }
	const bool &filt_hcallaser() { return lepton_tree_obj.filt_hcallaser(); }
	const bool &filt_ecaltp() { return lepton_tree_obj.filt_ecaltp(); }
	const bool &filt_trkfail() { return lepton_tree_obj.filt_trkfail(); }
	const bool &filt_eebadsc() { return lepton_tree_obj.filt_eebadsc(); }
	const bool &evt_isRealData() { return lepton_tree_obj.evt_isRealData(); }
	const float &scale1fb() { return lepton_tree_obj.scale1fb(); }
	const float &evt_xsec_incl() { return lepton_tree_obj.evt_xsec_incl(); }
	const float &evt_kfactor() { return lepton_tree_obj.evt_kfactor(); }
	const float &gen_met() { return lepton_tree_obj.gen_met(); }
	const float &gen_metPhi() { return lepton_tree_obj.gen_metPhi(); }
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
	const float &dxyPV() { return lepton_tree_obj.dxyPV(); }
	const float &dZ() { return lepton_tree_obj.dZ(); }
	const float &dxyPV_err() { return lepton_tree_obj.dxyPV_err(); }
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
	const float &sigmaIEtaIEta_full5x5() { return lepton_tree_obj.sigmaIEtaIEta_full5x5(); }
	const float &etaSC() { return lepton_tree_obj.etaSC(); }
	const float &dEtaIn() { return lepton_tree_obj.dEtaIn(); }
	const float &dPhiIn() { return lepton_tree_obj.dPhiIn(); }
	const float &hOverE() { return lepton_tree_obj.hOverE(); }
	const float &ecalEnergy() { return lepton_tree_obj.ecalEnergy(); }
	const float &eOverPIn() { return lepton_tree_obj.eOverPIn(); }
	const bool &conv_vtx_flag() { return lepton_tree_obj.conv_vtx_flag(); }
	const int &exp_innerlayers() { return lepton_tree_obj.exp_innerlayers(); }
	const int &charge() { return lepton_tree_obj.charge(); }
	const int &sccharge() { return lepton_tree_obj.sccharge(); }
	const int &ckf_charge() { return lepton_tree_obj.ckf_charge(); }
	const bool &threeChargeAgree() { return lepton_tree_obj.threeChargeAgree(); }
	const int &pid_PFMuon() { return lepton_tree_obj.pid_PFMuon(); }
	const float &gfit_chi2() { return lepton_tree_obj.gfit_chi2(); }
	const float &gfit_ndof() { return lepton_tree_obj.gfit_ndof(); }
	const int &gfit_validSTAHits() { return lepton_tree_obj.gfit_validSTAHits(); }
	const int &numberOfMatchedStations() { return lepton_tree_obj.numberOfMatchedStations(); }
	const int &validPixelHits() { return lepton_tree_obj.validPixelHits(); }
	const int &nlayers() { return lepton_tree_obj.nlayers(); }
}
