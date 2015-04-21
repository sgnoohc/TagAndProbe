// -*- C++ -*-
#ifndef LeptonTree_H
#define LeptonTree_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class LeptonTree {
private: 
protected: 
	unsigned int index;
	float	evt_pfmet_;
	TBranch *evt_pfmet_branch;
	bool evt_pfmet_isLoaded;
	float	evt_pfmetPhi_;
	TBranch *evt_pfmetPhi_branch;
	bool evt_pfmetPhi_isLoaded;
	int	evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	int	evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
	int	evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	bool	filt_csc_;
	TBranch *filt_csc_branch;
	bool filt_csc_isLoaded;
	bool	filt_hbhe_;
	TBranch *filt_hbhe_branch;
	bool filt_hbhe_isLoaded;
	bool	filt_hcallaser_;
	TBranch *filt_hcallaser_branch;
	bool filt_hcallaser_isLoaded;
	bool	filt_ecaltp_;
	TBranch *filt_ecaltp_branch;
	bool filt_ecaltp_isLoaded;
	bool	filt_trkfail_;
	TBranch *filt_trkfail_branch;
	bool filt_trkfail_isLoaded;
	bool	filt_eebadsc_;
	TBranch *filt_eebadsc_branch;
	bool filt_eebadsc_isLoaded;
	bool	evt_isRealData_;
	TBranch *evt_isRealData_branch;
	bool evt_isRealData_isLoaded;
	float	scale1fb_;
	TBranch *scale1fb_branch;
	bool scale1fb_isLoaded;
	float	evt_xsec_incl_;
	TBranch *evt_xsec_incl_branch;
	bool evt_xsec_incl_isLoaded;
	float	evt_kfactor_;
	TBranch *evt_kfactor_branch;
	bool evt_kfactor_isLoaded;
	float	gen_met_;
	TBranch *gen_met_branch;
	bool gen_met_isLoaded;
	float	gen_metPhi_;
	TBranch *gen_metPhi_branch;
	bool gen_metPhi_isLoaded;
	float	njets_;
	TBranch *njets_branch;
	bool njets_isLoaded;
	float	ht_;
	TBranch *ht_branch;
	bool ht_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *jets_;
	TBranch *jets_branch;
	bool jets_isLoaded;
	vector<float> *jets_disc_;
	TBranch *jets_disc_branch;
	bool jets_disc_isLoaded;
	TString *sample_;
	TBranch *sample_branch;
	bool sample_isLoaded;
	int	nFOs_;
	TBranch *nFOs_branch;
	bool nFOs_isLoaded;
	int	nvtx_;
	TBranch *nvtx_branch;
	bool nvtx_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *p4_;
	TBranch *p4_branch;
	bool p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *tag_p4_;
	TBranch *tag_p4_branch;
	bool tag_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *dilep_p4_;
	TBranch *dilep_p4_branch;
	bool dilep_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mc_p4_;
	TBranch *mc_p4_branch;
	bool mc_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mc_motherp4_;
	TBranch *mc_motherp4_branch;
	bool mc_motherp4_isLoaded;
	int	id_;
	TBranch *id_branch;
	bool id_isLoaded;
	int	idx_;
	TBranch *idx_branch;
	bool idx_isLoaded;
	float	dxyPV_;
	TBranch *dxyPV_branch;
	bool dxyPV_isLoaded;
	float	dZ_;
	TBranch *dZ_branch;
	bool dZ_isLoaded;
	float	dxyPV_err_;
	TBranch *dxyPV_err_branch;
	bool dxyPV_err_isLoaded;
	int	motherID_;
	TBranch *motherID_branch;
	bool motherID_isLoaded;
	int	mc_id_;
	TBranch *mc_id_branch;
	bool mc_id_isLoaded;
	float	iso_;
	TBranch *iso_branch;
	bool iso_isLoaded;
	bool	passes_id_;
	TBranch *passes_id_branch;
	bool passes_id_isLoaded;
	bool	passes_id_ptrel_;
	TBranch *passes_id_ptrel_branch;
	bool passes_id_ptrel_isLoaded;
	bool	passes_id_miniiso_;
	TBranch *passes_id_miniiso_branch;
	bool passes_id_miniiso_isLoaded;
	bool	passes_id_newminiiso_;
	TBranch *passes_id_newminiiso_branch;
	bool passes_id_newminiiso_isLoaded;
	bool	FO_;
	TBranch *FO_branch;
	bool FO_isLoaded;
	bool	FO_ptrel_;
	TBranch *FO_ptrel_branch;
	bool FO_ptrel_isLoaded;
	bool	FO_miniiso_;
	TBranch *FO_miniiso_branch;
	bool FO_miniiso_isLoaded;
	bool	FO_newminiiso_;
	TBranch *FO_newminiiso_branch;
	bool FO_newminiiso_isLoaded;
	bool	FO_NoIso_;
	TBranch *FO_NoIso_branch;
	bool FO_NoIso_isLoaded;
	float	ip3d_;
	TBranch *ip3d_branch;
	bool ip3d_isLoaded;
	float	ip3derr_;
	TBranch *ip3derr_branch;
	bool ip3derr_isLoaded;
	int	type_;
	TBranch *type_branch;
	bool type_isLoaded;
	float	mt_;
	TBranch *mt_branch;
	bool mt_isLoaded;
	float	ptrelv0_;
	TBranch *ptrelv0_branch;
	bool ptrelv0_isLoaded;
	float	ptrelv1_;
	TBranch *ptrelv1_branch;
	bool ptrelv1_isLoaded;
	float	miniiso_;
	TBranch *miniiso_branch;
	bool miniiso_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet_close_lep_;
	TBranch *jet_close_lep_branch;
	bool jet_close_lep_isLoaded;
	int	tag_charge_;
	TBranch *tag_charge_branch;
	bool tag_charge_isLoaded;
	float	dilep_mass_;
	TBranch *dilep_mass_branch;
	bool dilep_mass_isLoaded;
	float	sigmaIEtaIEta_full5x5_;
	TBranch *sigmaIEtaIEta_full5x5_branch;
	bool sigmaIEtaIEta_full5x5_isLoaded;
	float	etaSC_;
	TBranch *etaSC_branch;
	bool etaSC_isLoaded;
	float	dEtaIn_;
	TBranch *dEtaIn_branch;
	bool dEtaIn_isLoaded;
	float	dPhiIn_;
	TBranch *dPhiIn_branch;
	bool dPhiIn_isLoaded;
	float	hOverE_;
	TBranch *hOverE_branch;
	bool hOverE_isLoaded;
	float	ecalEnergy_;
	TBranch *ecalEnergy_branch;
	bool ecalEnergy_isLoaded;
	float	eOverPIn_;
	TBranch *eOverPIn_branch;
	bool eOverPIn_isLoaded;
	bool	conv_vtx_flag_;
	TBranch *conv_vtx_flag_branch;
	bool conv_vtx_flag_isLoaded;
	int	exp_innerlayers_;
	TBranch *exp_innerlayers_branch;
	bool exp_innerlayers_isLoaded;
	int	charge_;
	TBranch *charge_branch;
	bool charge_isLoaded;
	int	sccharge_;
	TBranch *sccharge_branch;
	bool sccharge_isLoaded;
	int	ckf_charge_;
	TBranch *ckf_charge_branch;
	bool ckf_charge_isLoaded;
	bool	threeChargeAgree_;
	TBranch *threeChargeAgree_branch;
	bool threeChargeAgree_isLoaded;
	int	pid_PFMuon_;
	TBranch *pid_PFMuon_branch;
	bool pid_PFMuon_isLoaded;
	float	gfit_chi2_;
	TBranch *gfit_chi2_branch;
	bool gfit_chi2_isLoaded;
	float	gfit_ndof_;
	TBranch *gfit_ndof_branch;
	bool gfit_ndof_isLoaded;
	int	gfit_validSTAHits_;
	TBranch *gfit_validSTAHits_branch;
	bool gfit_validSTAHits_isLoaded;
	int	numberOfMatchedStations_;
	TBranch *numberOfMatchedStations_branch;
	bool numberOfMatchedStations_isLoaded;
	int	validPixelHits_;
	TBranch *validPixelHits_branch;
	bool validPixelHits_isLoaded;
	int	nlayers_;
	TBranch *nlayers_branch;
	bool nlayers_isLoaded;
public: 
void Init(TTree *tree) {
	jets_branch = 0;
	if (tree->GetBranch("jets") != 0) {
		jets_branch = tree->GetBranch("jets");
		if (jets_branch) {jets_branch->SetAddress(&jets_);}
	}
	p4_branch = 0;
	if (tree->GetBranch("p4") != 0) {
		p4_branch = tree->GetBranch("p4");
		if (p4_branch) {p4_branch->SetAddress(&p4_);}
	}
	tag_p4_branch = 0;
	if (tree->GetBranch("tag_p4") != 0) {
		tag_p4_branch = tree->GetBranch("tag_p4");
		if (tag_p4_branch) {tag_p4_branch->SetAddress(&tag_p4_);}
	}
	dilep_p4_branch = 0;
	if (tree->GetBranch("dilep_p4") != 0) {
		dilep_p4_branch = tree->GetBranch("dilep_p4");
		if (dilep_p4_branch) {dilep_p4_branch->SetAddress(&dilep_p4_);}
	}
	mc_p4_branch = 0;
	if (tree->GetBranch("mc_p4") != 0) {
		mc_p4_branch = tree->GetBranch("mc_p4");
		if (mc_p4_branch) {mc_p4_branch->SetAddress(&mc_p4_);}
	}
	mc_motherp4_branch = 0;
	if (tree->GetBranch("mc_motherp4") != 0) {
		mc_motherp4_branch = tree->GetBranch("mc_motherp4");
		if (mc_motherp4_branch) {mc_motherp4_branch->SetAddress(&mc_motherp4_);}
	}
	jet_close_lep_branch = 0;
	if (tree->GetBranch("jet_close_lep") != 0) {
		jet_close_lep_branch = tree->GetBranch("jet_close_lep");
		if (jet_close_lep_branch) {jet_close_lep_branch->SetAddress(&jet_close_lep_);}
	}
  tree->SetMakeClass(1);
	evt_pfmet_branch = 0;
	if (tree->GetBranch("evt_pfmet") != 0) {
		evt_pfmet_branch = tree->GetBranch("evt_pfmet");
		if (evt_pfmet_branch) {evt_pfmet_branch->SetAddress(&evt_pfmet_);}
	}
	evt_pfmetPhi_branch = 0;
	if (tree->GetBranch("evt_pfmetPhi") != 0) {
		evt_pfmetPhi_branch = tree->GetBranch("evt_pfmetPhi");
		if (evt_pfmetPhi_branch) {evt_pfmetPhi_branch->SetAddress(&evt_pfmetPhi_);}
	}
	evt_event_branch = 0;
	if (tree->GetBranch("evt_event") != 0) {
		evt_event_branch = tree->GetBranch("evt_event");
		if (evt_event_branch) {evt_event_branch->SetAddress(&evt_event_);}
	}
	evt_lumiBlock_branch = 0;
	if (tree->GetBranch("evt_lumiBlock") != 0) {
		evt_lumiBlock_branch = tree->GetBranch("evt_lumiBlock");
		if (evt_lumiBlock_branch) {evt_lumiBlock_branch->SetAddress(&evt_lumiBlock_);}
	}
	evt_run_branch = 0;
	if (tree->GetBranch("evt_run") != 0) {
		evt_run_branch = tree->GetBranch("evt_run");
		if (evt_run_branch) {evt_run_branch->SetAddress(&evt_run_);}
	}
	filt_csc_branch = 0;
	if (tree->GetBranch("filt_csc") != 0) {
		filt_csc_branch = tree->GetBranch("filt_csc");
		if (filt_csc_branch) {filt_csc_branch->SetAddress(&filt_csc_);}
	}
	filt_hbhe_branch = 0;
	if (tree->GetBranch("filt_hbhe") != 0) {
		filt_hbhe_branch = tree->GetBranch("filt_hbhe");
		if (filt_hbhe_branch) {filt_hbhe_branch->SetAddress(&filt_hbhe_);}
	}
	filt_hcallaser_branch = 0;
	if (tree->GetBranch("filt_hcallaser") != 0) {
		filt_hcallaser_branch = tree->GetBranch("filt_hcallaser");
		if (filt_hcallaser_branch) {filt_hcallaser_branch->SetAddress(&filt_hcallaser_);}
	}
	filt_ecaltp_branch = 0;
	if (tree->GetBranch("filt_ecaltp") != 0) {
		filt_ecaltp_branch = tree->GetBranch("filt_ecaltp");
		if (filt_ecaltp_branch) {filt_ecaltp_branch->SetAddress(&filt_ecaltp_);}
	}
	filt_trkfail_branch = 0;
	if (tree->GetBranch("filt_trkfail") != 0) {
		filt_trkfail_branch = tree->GetBranch("filt_trkfail");
		if (filt_trkfail_branch) {filt_trkfail_branch->SetAddress(&filt_trkfail_);}
	}
	filt_eebadsc_branch = 0;
	if (tree->GetBranch("filt_eebadsc") != 0) {
		filt_eebadsc_branch = tree->GetBranch("filt_eebadsc");
		if (filt_eebadsc_branch) {filt_eebadsc_branch->SetAddress(&filt_eebadsc_);}
	}
	evt_isRealData_branch = 0;
	if (tree->GetBranch("evt_isRealData") != 0) {
		evt_isRealData_branch = tree->GetBranch("evt_isRealData");
		if (evt_isRealData_branch) {evt_isRealData_branch->SetAddress(&evt_isRealData_);}
	}
	scale1fb_branch = 0;
	if (tree->GetBranch("scale1fb") != 0) {
		scale1fb_branch = tree->GetBranch("scale1fb");
		if (scale1fb_branch) {scale1fb_branch->SetAddress(&scale1fb_);}
	}
	evt_xsec_incl_branch = 0;
	if (tree->GetBranch("evt_xsec_incl") != 0) {
		evt_xsec_incl_branch = tree->GetBranch("evt_xsec_incl");
		if (evt_xsec_incl_branch) {evt_xsec_incl_branch->SetAddress(&evt_xsec_incl_);}
	}
	evt_kfactor_branch = 0;
	if (tree->GetBranch("evt_kfactor") != 0) {
		evt_kfactor_branch = tree->GetBranch("evt_kfactor");
		if (evt_kfactor_branch) {evt_kfactor_branch->SetAddress(&evt_kfactor_);}
	}
	gen_met_branch = 0;
	if (tree->GetBranch("gen_met") != 0) {
		gen_met_branch = tree->GetBranch("gen_met");
		if (gen_met_branch) {gen_met_branch->SetAddress(&gen_met_);}
	}
	gen_metPhi_branch = 0;
	if (tree->GetBranch("gen_metPhi") != 0) {
		gen_metPhi_branch = tree->GetBranch("gen_metPhi");
		if (gen_metPhi_branch) {gen_metPhi_branch->SetAddress(&gen_metPhi_);}
	}
	njets_branch = 0;
	if (tree->GetBranch("njets") != 0) {
		njets_branch = tree->GetBranch("njets");
		if (njets_branch) {njets_branch->SetAddress(&njets_);}
	}
	ht_branch = 0;
	if (tree->GetBranch("ht") != 0) {
		ht_branch = tree->GetBranch("ht");
		if (ht_branch) {ht_branch->SetAddress(&ht_);}
	}
	jets_disc_branch = 0;
	if (tree->GetBranch("jets_disc") != 0) {
		jets_disc_branch = tree->GetBranch("jets_disc");
		if (jets_disc_branch) {jets_disc_branch->SetAddress(&jets_disc_);}
	}
	sample_branch = 0;
	if (tree->GetBranch("sample") != 0) {
		sample_branch = tree->GetBranch("sample");
		if (sample_branch) {sample_branch->SetAddress(&sample_);}
	}
	nFOs_branch = 0;
	if (tree->GetBranch("nFOs") != 0) {
		nFOs_branch = tree->GetBranch("nFOs");
		if (nFOs_branch) {nFOs_branch->SetAddress(&nFOs_);}
	}
	nvtx_branch = 0;
	if (tree->GetBranch("nvtx") != 0) {
		nvtx_branch = tree->GetBranch("nvtx");
		if (nvtx_branch) {nvtx_branch->SetAddress(&nvtx_);}
	}
	id_branch = 0;
	if (tree->GetBranch("id") != 0) {
		id_branch = tree->GetBranch("id");
		if (id_branch) {id_branch->SetAddress(&id_);}
	}
	idx_branch = 0;
	if (tree->GetBranch("idx") != 0) {
		idx_branch = tree->GetBranch("idx");
		if (idx_branch) {idx_branch->SetAddress(&idx_);}
	}
	dxyPV_branch = 0;
	if (tree->GetBranch("dxyPV") != 0) {
		dxyPV_branch = tree->GetBranch("dxyPV");
		if (dxyPV_branch) {dxyPV_branch->SetAddress(&dxyPV_);}
	}
	dZ_branch = 0;
	if (tree->GetBranch("dZ") != 0) {
		dZ_branch = tree->GetBranch("dZ");
		if (dZ_branch) {dZ_branch->SetAddress(&dZ_);}
	}
	dxyPV_err_branch = 0;
	if (tree->GetBranch("dxyPV_err") != 0) {
		dxyPV_err_branch = tree->GetBranch("dxyPV_err");
		if (dxyPV_err_branch) {dxyPV_err_branch->SetAddress(&dxyPV_err_);}
	}
	motherID_branch = 0;
	if (tree->GetBranch("motherID") != 0) {
		motherID_branch = tree->GetBranch("motherID");
		if (motherID_branch) {motherID_branch->SetAddress(&motherID_);}
	}
	mc_id_branch = 0;
	if (tree->GetBranch("mc_id") != 0) {
		mc_id_branch = tree->GetBranch("mc_id");
		if (mc_id_branch) {mc_id_branch->SetAddress(&mc_id_);}
	}
	iso_branch = 0;
	if (tree->GetBranch("iso") != 0) {
		iso_branch = tree->GetBranch("iso");
		if (iso_branch) {iso_branch->SetAddress(&iso_);}
	}
	passes_id_branch = 0;
	if (tree->GetBranch("passes_id") != 0) {
		passes_id_branch = tree->GetBranch("passes_id");
		if (passes_id_branch) {passes_id_branch->SetAddress(&passes_id_);}
	}
	passes_id_ptrel_branch = 0;
	if (tree->GetBranch("passes_id_ptrel") != 0) {
		passes_id_ptrel_branch = tree->GetBranch("passes_id_ptrel");
		if (passes_id_ptrel_branch) {passes_id_ptrel_branch->SetAddress(&passes_id_ptrel_);}
	}
	passes_id_miniiso_branch = 0;
	if (tree->GetBranch("passes_id_miniiso") != 0) {
		passes_id_miniiso_branch = tree->GetBranch("passes_id_miniiso");
		if (passes_id_miniiso_branch) {passes_id_miniiso_branch->SetAddress(&passes_id_miniiso_);}
	}
	passes_id_newminiiso_branch = 0;
	if (tree->GetBranch("passes_id_newminiiso") != 0) {
		passes_id_newminiiso_branch = tree->GetBranch("passes_id_newminiiso");
		if (passes_id_newminiiso_branch) {passes_id_newminiiso_branch->SetAddress(&passes_id_newminiiso_);}
	}
	FO_branch = 0;
	if (tree->GetBranch("FO") != 0) {
		FO_branch = tree->GetBranch("FO");
		if (FO_branch) {FO_branch->SetAddress(&FO_);}
	}
	FO_ptrel_branch = 0;
	if (tree->GetBranch("FO_ptrel") != 0) {
		FO_ptrel_branch = tree->GetBranch("FO_ptrel");
		if (FO_ptrel_branch) {FO_ptrel_branch->SetAddress(&FO_ptrel_);}
	}
	FO_miniiso_branch = 0;
	if (tree->GetBranch("FO_miniiso") != 0) {
		FO_miniiso_branch = tree->GetBranch("FO_miniiso");
		if (FO_miniiso_branch) {FO_miniiso_branch->SetAddress(&FO_miniiso_);}
	}
	FO_newminiiso_branch = 0;
	if (tree->GetBranch("FO_newminiiso") != 0) {
		FO_newminiiso_branch = tree->GetBranch("FO_newminiiso");
		if (FO_newminiiso_branch) {FO_newminiiso_branch->SetAddress(&FO_newminiiso_);}
	}
	FO_NoIso_branch = 0;
	if (tree->GetBranch("FO_NoIso") != 0) {
		FO_NoIso_branch = tree->GetBranch("FO_NoIso");
		if (FO_NoIso_branch) {FO_NoIso_branch->SetAddress(&FO_NoIso_);}
	}
	ip3d_branch = 0;
	if (tree->GetBranch("ip3d") != 0) {
		ip3d_branch = tree->GetBranch("ip3d");
		if (ip3d_branch) {ip3d_branch->SetAddress(&ip3d_);}
	}
	ip3derr_branch = 0;
	if (tree->GetBranch("ip3derr") != 0) {
		ip3derr_branch = tree->GetBranch("ip3derr");
		if (ip3derr_branch) {ip3derr_branch->SetAddress(&ip3derr_);}
	}
	type_branch = 0;
	if (tree->GetBranch("type") != 0) {
		type_branch = tree->GetBranch("type");
		if (type_branch) {type_branch->SetAddress(&type_);}
	}
	mt_branch = 0;
	if (tree->GetBranch("mt") != 0) {
		mt_branch = tree->GetBranch("mt");
		if (mt_branch) {mt_branch->SetAddress(&mt_);}
	}
	ptrelv0_branch = 0;
	if (tree->GetBranch("ptrelv0") != 0) {
		ptrelv0_branch = tree->GetBranch("ptrelv0");
		if (ptrelv0_branch) {ptrelv0_branch->SetAddress(&ptrelv0_);}
	}
	ptrelv1_branch = 0;
	if (tree->GetBranch("ptrelv1") != 0) {
		ptrelv1_branch = tree->GetBranch("ptrelv1");
		if (ptrelv1_branch) {ptrelv1_branch->SetAddress(&ptrelv1_);}
	}
	miniiso_branch = 0;
	if (tree->GetBranch("miniiso") != 0) {
		miniiso_branch = tree->GetBranch("miniiso");
		if (miniiso_branch) {miniiso_branch->SetAddress(&miniiso_);}
	}
	tag_charge_branch = 0;
	if (tree->GetBranch("tag_charge") != 0) {
		tag_charge_branch = tree->GetBranch("tag_charge");
		if (tag_charge_branch) {tag_charge_branch->SetAddress(&tag_charge_);}
	}
	dilep_mass_branch = 0;
	if (tree->GetBranch("dilep_mass") != 0) {
		dilep_mass_branch = tree->GetBranch("dilep_mass");
		if (dilep_mass_branch) {dilep_mass_branch->SetAddress(&dilep_mass_);}
	}
	sigmaIEtaIEta_full5x5_branch = 0;
	if (tree->GetBranch("sigmaIEtaIEta_full5x5") != 0) {
		sigmaIEtaIEta_full5x5_branch = tree->GetBranch("sigmaIEtaIEta_full5x5");
		if (sigmaIEtaIEta_full5x5_branch) {sigmaIEtaIEta_full5x5_branch->SetAddress(&sigmaIEtaIEta_full5x5_);}
	}
	etaSC_branch = 0;
	if (tree->GetBranch("etaSC") != 0) {
		etaSC_branch = tree->GetBranch("etaSC");
		if (etaSC_branch) {etaSC_branch->SetAddress(&etaSC_);}
	}
	dEtaIn_branch = 0;
	if (tree->GetBranch("dEtaIn") != 0) {
		dEtaIn_branch = tree->GetBranch("dEtaIn");
		if (dEtaIn_branch) {dEtaIn_branch->SetAddress(&dEtaIn_);}
	}
	dPhiIn_branch = 0;
	if (tree->GetBranch("dPhiIn") != 0) {
		dPhiIn_branch = tree->GetBranch("dPhiIn");
		if (dPhiIn_branch) {dPhiIn_branch->SetAddress(&dPhiIn_);}
	}
	hOverE_branch = 0;
	if (tree->GetBranch("hOverE") != 0) {
		hOverE_branch = tree->GetBranch("hOverE");
		if (hOverE_branch) {hOverE_branch->SetAddress(&hOverE_);}
	}
	ecalEnergy_branch = 0;
	if (tree->GetBranch("ecalEnergy") != 0) {
		ecalEnergy_branch = tree->GetBranch("ecalEnergy");
		if (ecalEnergy_branch) {ecalEnergy_branch->SetAddress(&ecalEnergy_);}
	}
	eOverPIn_branch = 0;
	if (tree->GetBranch("eOverPIn") != 0) {
		eOverPIn_branch = tree->GetBranch("eOverPIn");
		if (eOverPIn_branch) {eOverPIn_branch->SetAddress(&eOverPIn_);}
	}
	conv_vtx_flag_branch = 0;
	if (tree->GetBranch("conv_vtx_flag") != 0) {
		conv_vtx_flag_branch = tree->GetBranch("conv_vtx_flag");
		if (conv_vtx_flag_branch) {conv_vtx_flag_branch->SetAddress(&conv_vtx_flag_);}
	}
	exp_innerlayers_branch = 0;
	if (tree->GetBranch("exp_innerlayers") != 0) {
		exp_innerlayers_branch = tree->GetBranch("exp_innerlayers");
		if (exp_innerlayers_branch) {exp_innerlayers_branch->SetAddress(&exp_innerlayers_);}
	}
	charge_branch = 0;
	if (tree->GetBranch("charge") != 0) {
		charge_branch = tree->GetBranch("charge");
		if (charge_branch) {charge_branch->SetAddress(&charge_);}
	}
	sccharge_branch = 0;
	if (tree->GetBranch("sccharge") != 0) {
		sccharge_branch = tree->GetBranch("sccharge");
		if (sccharge_branch) {sccharge_branch->SetAddress(&sccharge_);}
	}
	ckf_charge_branch = 0;
	if (tree->GetBranch("ckf_charge") != 0) {
		ckf_charge_branch = tree->GetBranch("ckf_charge");
		if (ckf_charge_branch) {ckf_charge_branch->SetAddress(&ckf_charge_);}
	}
	threeChargeAgree_branch = 0;
	if (tree->GetBranch("threeChargeAgree") != 0) {
		threeChargeAgree_branch = tree->GetBranch("threeChargeAgree");
		if (threeChargeAgree_branch) {threeChargeAgree_branch->SetAddress(&threeChargeAgree_);}
	}
	pid_PFMuon_branch = 0;
	if (tree->GetBranch("pid_PFMuon") != 0) {
		pid_PFMuon_branch = tree->GetBranch("pid_PFMuon");
		if (pid_PFMuon_branch) {pid_PFMuon_branch->SetAddress(&pid_PFMuon_);}
	}
	gfit_chi2_branch = 0;
	if (tree->GetBranch("gfit_chi2") != 0) {
		gfit_chi2_branch = tree->GetBranch("gfit_chi2");
		if (gfit_chi2_branch) {gfit_chi2_branch->SetAddress(&gfit_chi2_);}
	}
	gfit_ndof_branch = 0;
	if (tree->GetBranch("gfit_ndof") != 0) {
		gfit_ndof_branch = tree->GetBranch("gfit_ndof");
		if (gfit_ndof_branch) {gfit_ndof_branch->SetAddress(&gfit_ndof_);}
	}
	gfit_validSTAHits_branch = 0;
	if (tree->GetBranch("gfit_validSTAHits") != 0) {
		gfit_validSTAHits_branch = tree->GetBranch("gfit_validSTAHits");
		if (gfit_validSTAHits_branch) {gfit_validSTAHits_branch->SetAddress(&gfit_validSTAHits_);}
	}
	numberOfMatchedStations_branch = 0;
	if (tree->GetBranch("numberOfMatchedStations") != 0) {
		numberOfMatchedStations_branch = tree->GetBranch("numberOfMatchedStations");
		if (numberOfMatchedStations_branch) {numberOfMatchedStations_branch->SetAddress(&numberOfMatchedStations_);}
	}
	validPixelHits_branch = 0;
	if (tree->GetBranch("validPixelHits") != 0) {
		validPixelHits_branch = tree->GetBranch("validPixelHits");
		if (validPixelHits_branch) {validPixelHits_branch->SetAddress(&validPixelHits_);}
	}
	nlayers_branch = 0;
	if (tree->GetBranch("nlayers") != 0) {
		nlayers_branch = tree->GetBranch("nlayers");
		if (nlayers_branch) {nlayers_branch->SetAddress(&nlayers_);}
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		evt_pfmet_isLoaded = false;
		evt_pfmetPhi_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		filt_csc_isLoaded = false;
		filt_hbhe_isLoaded = false;
		filt_hcallaser_isLoaded = false;
		filt_ecaltp_isLoaded = false;
		filt_trkfail_isLoaded = false;
		filt_eebadsc_isLoaded = false;
		evt_isRealData_isLoaded = false;
		scale1fb_isLoaded = false;
		evt_xsec_incl_isLoaded = false;
		evt_kfactor_isLoaded = false;
		gen_met_isLoaded = false;
		gen_metPhi_isLoaded = false;
		njets_isLoaded = false;
		ht_isLoaded = false;
		jets_isLoaded = false;
		jets_disc_isLoaded = false;
		sample_isLoaded = false;
		nFOs_isLoaded = false;
		nvtx_isLoaded = false;
		p4_isLoaded = false;
		tag_p4_isLoaded = false;
		dilep_p4_isLoaded = false;
		mc_p4_isLoaded = false;
		mc_motherp4_isLoaded = false;
		id_isLoaded = false;
		idx_isLoaded = false;
		dxyPV_isLoaded = false;
		dZ_isLoaded = false;
		dxyPV_err_isLoaded = false;
		motherID_isLoaded = false;
		mc_id_isLoaded = false;
		iso_isLoaded = false;
		passes_id_isLoaded = false;
		passes_id_ptrel_isLoaded = false;
		passes_id_miniiso_isLoaded = false;
		passes_id_newminiiso_isLoaded = false;
		FO_isLoaded = false;
		FO_ptrel_isLoaded = false;
		FO_miniiso_isLoaded = false;
		FO_newminiiso_isLoaded = false;
		FO_NoIso_isLoaded = false;
		ip3d_isLoaded = false;
		ip3derr_isLoaded = false;
		type_isLoaded = false;
		mt_isLoaded = false;
		ptrelv0_isLoaded = false;
		ptrelv1_isLoaded = false;
		miniiso_isLoaded = false;
		jet_close_lep_isLoaded = false;
		tag_charge_isLoaded = false;
		dilep_mass_isLoaded = false;
		sigmaIEtaIEta_full5x5_isLoaded = false;
		etaSC_isLoaded = false;
		dEtaIn_isLoaded = false;
		dPhiIn_isLoaded = false;
		hOverE_isLoaded = false;
		ecalEnergy_isLoaded = false;
		eOverPIn_isLoaded = false;
		conv_vtx_flag_isLoaded = false;
		exp_innerlayers_isLoaded = false;
		charge_isLoaded = false;
		sccharge_isLoaded = false;
		ckf_charge_isLoaded = false;
		threeChargeAgree_isLoaded = false;
		pid_PFMuon_isLoaded = false;
		gfit_chi2_isLoaded = false;
		gfit_ndof_isLoaded = false;
		gfit_validSTAHits_isLoaded = false;
		numberOfMatchedStations_isLoaded = false;
		validPixelHits_isLoaded = false;
		nlayers_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (evt_pfmet_branch != 0) evt_pfmet();
	if (evt_pfmetPhi_branch != 0) evt_pfmetPhi();
	if (evt_event_branch != 0) evt_event();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_run_branch != 0) evt_run();
	if (filt_csc_branch != 0) filt_csc();
	if (filt_hbhe_branch != 0) filt_hbhe();
	if (filt_hcallaser_branch != 0) filt_hcallaser();
	if (filt_ecaltp_branch != 0) filt_ecaltp();
	if (filt_trkfail_branch != 0) filt_trkfail();
	if (filt_eebadsc_branch != 0) filt_eebadsc();
	if (evt_isRealData_branch != 0) evt_isRealData();
	if (scale1fb_branch != 0) scale1fb();
	if (evt_xsec_incl_branch != 0) evt_xsec_incl();
	if (evt_kfactor_branch != 0) evt_kfactor();
	if (gen_met_branch != 0) gen_met();
	if (gen_metPhi_branch != 0) gen_metPhi();
	if (njets_branch != 0) njets();
	if (ht_branch != 0) ht();
	if (jets_branch != 0) jets();
	if (jets_disc_branch != 0) jets_disc();
	if (sample_branch != 0) sample();
	if (nFOs_branch != 0) nFOs();
	if (nvtx_branch != 0) nvtx();
	if (p4_branch != 0) p4();
	if (tag_p4_branch != 0) tag_p4();
	if (dilep_p4_branch != 0) dilep_p4();
	if (mc_p4_branch != 0) mc_p4();
	if (mc_motherp4_branch != 0) mc_motherp4();
	if (id_branch != 0) id();
	if (idx_branch != 0) idx();
	if (dxyPV_branch != 0) dxyPV();
	if (dZ_branch != 0) dZ();
	if (dxyPV_err_branch != 0) dxyPV_err();
	if (motherID_branch != 0) motherID();
	if (mc_id_branch != 0) mc_id();
	if (iso_branch != 0) iso();
	if (passes_id_branch != 0) passes_id();
	if (passes_id_ptrel_branch != 0) passes_id_ptrel();
	if (passes_id_miniiso_branch != 0) passes_id_miniiso();
	if (passes_id_newminiiso_branch != 0) passes_id_newminiiso();
	if (FO_branch != 0) FO();
	if (FO_ptrel_branch != 0) FO_ptrel();
	if (FO_miniiso_branch != 0) FO_miniiso();
	if (FO_newminiiso_branch != 0) FO_newminiiso();
	if (FO_NoIso_branch != 0) FO_NoIso();
	if (ip3d_branch != 0) ip3d();
	if (ip3derr_branch != 0) ip3derr();
	if (type_branch != 0) type();
	if (mt_branch != 0) mt();
	if (ptrelv0_branch != 0) ptrelv0();
	if (ptrelv1_branch != 0) ptrelv1();
	if (miniiso_branch != 0) miniiso();
	if (jet_close_lep_branch != 0) jet_close_lep();
	if (tag_charge_branch != 0) tag_charge();
	if (dilep_mass_branch != 0) dilep_mass();
	if (sigmaIEtaIEta_full5x5_branch != 0) sigmaIEtaIEta_full5x5();
	if (etaSC_branch != 0) etaSC();
	if (dEtaIn_branch != 0) dEtaIn();
	if (dPhiIn_branch != 0) dPhiIn();
	if (hOverE_branch != 0) hOverE();
	if (ecalEnergy_branch != 0) ecalEnergy();
	if (eOverPIn_branch != 0) eOverPIn();
	if (conv_vtx_flag_branch != 0) conv_vtx_flag();
	if (exp_innerlayers_branch != 0) exp_innerlayers();
	if (charge_branch != 0) charge();
	if (sccharge_branch != 0) sccharge();
	if (ckf_charge_branch != 0) ckf_charge();
	if (threeChargeAgree_branch != 0) threeChargeAgree();
	if (pid_PFMuon_branch != 0) pid_PFMuon();
	if (gfit_chi2_branch != 0) gfit_chi2();
	if (gfit_ndof_branch != 0) gfit_ndof();
	if (gfit_validSTAHits_branch != 0) gfit_validSTAHits();
	if (numberOfMatchedStations_branch != 0) numberOfMatchedStations();
	if (validPixelHits_branch != 0) validPixelHits();
	if (nlayers_branch != 0) nlayers();
}

	float &evt_pfmet()
	{
		if (not evt_pfmet_isLoaded) {
			if (evt_pfmet_branch != 0) {
				evt_pfmet_branch->GetEntry(index);
			} else { 
				printf("branch evt_pfmet_branch does not exist!\n");
				exit(1);
			}
			evt_pfmet_isLoaded = true;
		}
		return evt_pfmet_;
	}
	float &evt_pfmetPhi()
	{
		if (not evt_pfmetPhi_isLoaded) {
			if (evt_pfmetPhi_branch != 0) {
				evt_pfmetPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_pfmetPhi_branch does not exist!\n");
				exit(1);
			}
			evt_pfmetPhi_isLoaded = true;
		}
		return evt_pfmetPhi_;
	}
	int &evt_event()
	{
		if (not evt_event_isLoaded) {
			if (evt_event_branch != 0) {
				evt_event_branch->GetEntry(index);
			} else { 
				printf("branch evt_event_branch does not exist!\n");
				exit(1);
			}
			evt_event_isLoaded = true;
		}
		return evt_event_;
	}
	int &evt_lumiBlock()
	{
		if (not evt_lumiBlock_isLoaded) {
			if (evt_lumiBlock_branch != 0) {
				evt_lumiBlock_branch->GetEntry(index);
			} else { 
				printf("branch evt_lumiBlock_branch does not exist!\n");
				exit(1);
			}
			evt_lumiBlock_isLoaded = true;
		}
		return evt_lumiBlock_;
	}
	int &evt_run()
	{
		if (not evt_run_isLoaded) {
			if (evt_run_branch != 0) {
				evt_run_branch->GetEntry(index);
			} else { 
				printf("branch evt_run_branch does not exist!\n");
				exit(1);
			}
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	bool &	filt_csc()
	{
		if (not filt_csc_isLoaded) {
			if (filt_csc_branch != 0) {
				filt_csc_branch->GetEntry(index);
			} else { 
				printf("branch filt_csc_branch does not exist!\n");
				exit(1);
			}
			filt_csc_isLoaded = true;
		}
		return filt_csc_;
	}
	bool &	filt_hbhe()
	{
		if (not filt_hbhe_isLoaded) {
			if (filt_hbhe_branch != 0) {
				filt_hbhe_branch->GetEntry(index);
			} else { 
				printf("branch filt_hbhe_branch does not exist!\n");
				exit(1);
			}
			filt_hbhe_isLoaded = true;
		}
		return filt_hbhe_;
	}
	bool &	filt_hcallaser()
	{
		if (not filt_hcallaser_isLoaded) {
			if (filt_hcallaser_branch != 0) {
				filt_hcallaser_branch->GetEntry(index);
			} else { 
				printf("branch filt_hcallaser_branch does not exist!\n");
				exit(1);
			}
			filt_hcallaser_isLoaded = true;
		}
		return filt_hcallaser_;
	}
	bool &	filt_ecaltp()
	{
		if (not filt_ecaltp_isLoaded) {
			if (filt_ecaltp_branch != 0) {
				filt_ecaltp_branch->GetEntry(index);
			} else { 
				printf("branch filt_ecaltp_branch does not exist!\n");
				exit(1);
			}
			filt_ecaltp_isLoaded = true;
		}
		return filt_ecaltp_;
	}
	bool &	filt_trkfail()
	{
		if (not filt_trkfail_isLoaded) {
			if (filt_trkfail_branch != 0) {
				filt_trkfail_branch->GetEntry(index);
			} else { 
				printf("branch filt_trkfail_branch does not exist!\n");
				exit(1);
			}
			filt_trkfail_isLoaded = true;
		}
		return filt_trkfail_;
	}
	bool &	filt_eebadsc()
	{
		if (not filt_eebadsc_isLoaded) {
			if (filt_eebadsc_branch != 0) {
				filt_eebadsc_branch->GetEntry(index);
			} else { 
				printf("branch filt_eebadsc_branch does not exist!\n");
				exit(1);
			}
			filt_eebadsc_isLoaded = true;
		}
		return filt_eebadsc_;
	}
	bool &	evt_isRealData()
	{
		if (not evt_isRealData_isLoaded) {
			if (evt_isRealData_branch != 0) {
				evt_isRealData_branch->GetEntry(index);
			} else { 
				printf("branch evt_isRealData_branch does not exist!\n");
				exit(1);
			}
			evt_isRealData_isLoaded = true;
		}
		return evt_isRealData_;
	}
	float &scale1fb()
	{
		if (not scale1fb_isLoaded) {
			if (scale1fb_branch != 0) {
				scale1fb_branch->GetEntry(index);
			} else { 
				printf("branch scale1fb_branch does not exist!\n");
				exit(1);
			}
			scale1fb_isLoaded = true;
		}
		return scale1fb_;
	}
	float &evt_xsec_incl()
	{
		if (not evt_xsec_incl_isLoaded) {
			if (evt_xsec_incl_branch != 0) {
				evt_xsec_incl_branch->GetEntry(index);
			} else { 
				printf("branch evt_xsec_incl_branch does not exist!\n");
				exit(1);
			}
			evt_xsec_incl_isLoaded = true;
		}
		return evt_xsec_incl_;
	}
	float &evt_kfactor()
	{
		if (not evt_kfactor_isLoaded) {
			if (evt_kfactor_branch != 0) {
				evt_kfactor_branch->GetEntry(index);
			} else { 
				printf("branch evt_kfactor_branch does not exist!\n");
				exit(1);
			}
			evt_kfactor_isLoaded = true;
		}
		return evt_kfactor_;
	}
	float &gen_met()
	{
		if (not gen_met_isLoaded) {
			if (gen_met_branch != 0) {
				gen_met_branch->GetEntry(index);
			} else { 
				printf("branch gen_met_branch does not exist!\n");
				exit(1);
			}
			gen_met_isLoaded = true;
		}
		return gen_met_;
	}
	float &gen_metPhi()
	{
		if (not gen_metPhi_isLoaded) {
			if (gen_metPhi_branch != 0) {
				gen_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch gen_metPhi_branch does not exist!\n");
				exit(1);
			}
			gen_metPhi_isLoaded = true;
		}
		return gen_metPhi_;
	}
	float &njets()
	{
		if (not njets_isLoaded) {
			if (njets_branch != 0) {
				njets_branch->GetEntry(index);
			} else { 
				printf("branch njets_branch does not exist!\n");
				exit(1);
			}
			njets_isLoaded = true;
		}
		return njets_;
	}
	float &ht()
	{
		if (not ht_isLoaded) {
			if (ht_branch != 0) {
				ht_branch->GetEntry(index);
			} else { 
				printf("branch ht_branch does not exist!\n");
				exit(1);
			}
			ht_isLoaded = true;
		}
		return ht_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets()
	{
		if (not jets_isLoaded) {
			if (jets_branch != 0) {
				jets_branch->GetEntry(index);
			} else { 
				printf("branch jets_branch does not exist!\n");
				exit(1);
			}
			jets_isLoaded = true;
		}
		return *jets_;
	}
	const vector<float> &jets_disc()
	{
		if (not jets_disc_isLoaded) {
			if (jets_disc_branch != 0) {
				jets_disc_branch->GetEntry(index);
			} else { 
				printf("branch jets_disc_branch does not exist!\n");
				exit(1);
			}
			jets_disc_isLoaded = true;
		}
		return *jets_disc_;
	}
	TString &sample()
	{
		if (not sample_isLoaded) {
			if (sample_branch != 0) {
				sample_branch->GetEntry(index);
			} else { 
				printf("branch sample_branch does not exist!\n");
				exit(1);
			}
			sample_isLoaded = true;
		}
		return *sample_;
	}
	int &nFOs()
	{
		if (not nFOs_isLoaded) {
			if (nFOs_branch != 0) {
				nFOs_branch->GetEntry(index);
			} else { 
				printf("branch nFOs_branch does not exist!\n");
				exit(1);
			}
			nFOs_isLoaded = true;
		}
		return nFOs_;
	}
	int &nvtx()
	{
		if (not nvtx_isLoaded) {
			if (nvtx_branch != 0) {
				nvtx_branch->GetEntry(index);
			} else { 
				printf("branch nvtx_branch does not exist!\n");
				exit(1);
			}
			nvtx_isLoaded = true;
		}
		return nvtx_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &p4()
	{
		if (not p4_isLoaded) {
			if (p4_branch != 0) {
				p4_branch->GetEntry(index);
			} else { 
				printf("branch p4_branch does not exist!\n");
				exit(1);
			}
			p4_isLoaded = true;
		}
		return *p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &tag_p4()
	{
		if (not tag_p4_isLoaded) {
			if (tag_p4_branch != 0) {
				tag_p4_branch->GetEntry(index);
			} else { 
				printf("branch tag_p4_branch does not exist!\n");
				exit(1);
			}
			tag_p4_isLoaded = true;
		}
		return *tag_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &dilep_p4()
	{
		if (not dilep_p4_isLoaded) {
			if (dilep_p4_branch != 0) {
				dilep_p4_branch->GetEntry(index);
			} else { 
				printf("branch dilep_p4_branch does not exist!\n");
				exit(1);
			}
			dilep_p4_isLoaded = true;
		}
		return *dilep_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mc_p4()
	{
		if (not mc_p4_isLoaded) {
			if (mc_p4_branch != 0) {
				mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch mc_p4_branch does not exist!\n");
				exit(1);
			}
			mc_p4_isLoaded = true;
		}
		return *mc_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mc_motherp4()
	{
		if (not mc_motherp4_isLoaded) {
			if (mc_motherp4_branch != 0) {
				mc_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch mc_motherp4_branch does not exist!\n");
				exit(1);
			}
			mc_motherp4_isLoaded = true;
		}
		return *mc_motherp4_;
	}
	int &id()
	{
		if (not id_isLoaded) {
			if (id_branch != 0) {
				id_branch->GetEntry(index);
			} else { 
				printf("branch id_branch does not exist!\n");
				exit(1);
			}
			id_isLoaded = true;
		}
		return id_;
	}
	int &idx()
	{
		if (not idx_isLoaded) {
			if (idx_branch != 0) {
				idx_branch->GetEntry(index);
			} else { 
				printf("branch idx_branch does not exist!\n");
				exit(1);
			}
			idx_isLoaded = true;
		}
		return idx_;
	}
	float &dxyPV()
	{
		if (not dxyPV_isLoaded) {
			if (dxyPV_branch != 0) {
				dxyPV_branch->GetEntry(index);
			} else { 
				printf("branch dxyPV_branch does not exist!\n");
				exit(1);
			}
			dxyPV_isLoaded = true;
		}
		return dxyPV_;
	}
	float &dZ()
	{
		if (not dZ_isLoaded) {
			if (dZ_branch != 0) {
				dZ_branch->GetEntry(index);
			} else { 
				printf("branch dZ_branch does not exist!\n");
				exit(1);
			}
			dZ_isLoaded = true;
		}
		return dZ_;
	}
	float &dxyPV_err()
	{
		if (not dxyPV_err_isLoaded) {
			if (dxyPV_err_branch != 0) {
				dxyPV_err_branch->GetEntry(index);
			} else { 
				printf("branch dxyPV_err_branch does not exist!\n");
				exit(1);
			}
			dxyPV_err_isLoaded = true;
		}
		return dxyPV_err_;
	}
	int &motherID()
	{
		if (not motherID_isLoaded) {
			if (motherID_branch != 0) {
				motherID_branch->GetEntry(index);
			} else { 
				printf("branch motherID_branch does not exist!\n");
				exit(1);
			}
			motherID_isLoaded = true;
		}
		return motherID_;
	}
	int &mc_id()
	{
		if (not mc_id_isLoaded) {
			if (mc_id_branch != 0) {
				mc_id_branch->GetEntry(index);
			} else { 
				printf("branch mc_id_branch does not exist!\n");
				exit(1);
			}
			mc_id_isLoaded = true;
		}
		return mc_id_;
	}
	float &iso()
	{
		if (not iso_isLoaded) {
			if (iso_branch != 0) {
				iso_branch->GetEntry(index);
			} else { 
				printf("branch iso_branch does not exist!\n");
				exit(1);
			}
			iso_isLoaded = true;
		}
		return iso_;
	}
	bool &	passes_id()
	{
		if (not passes_id_isLoaded) {
			if (passes_id_branch != 0) {
				passes_id_branch->GetEntry(index);
			} else { 
				printf("branch passes_id_branch does not exist!\n");
				exit(1);
			}
			passes_id_isLoaded = true;
		}
		return passes_id_;
	}
	bool &	passes_id_ptrel()
	{
		if (not passes_id_ptrel_isLoaded) {
			if (passes_id_ptrel_branch != 0) {
				passes_id_ptrel_branch->GetEntry(index);
			} else { 
				printf("branch passes_id_ptrel_branch does not exist!\n");
				exit(1);
			}
			passes_id_ptrel_isLoaded = true;
		}
		return passes_id_ptrel_;
	}
	bool &	passes_id_miniiso()
	{
		if (not passes_id_miniiso_isLoaded) {
			if (passes_id_miniiso_branch != 0) {
				passes_id_miniiso_branch->GetEntry(index);
			} else { 
				printf("branch passes_id_miniiso_branch does not exist!\n");
				exit(1);
			}
			passes_id_miniiso_isLoaded = true;
		}
		return passes_id_miniiso_;
	}
	bool &	passes_id_newminiiso()
	{
		if (not passes_id_newminiiso_isLoaded) {
			if (passes_id_newminiiso_branch != 0) {
				passes_id_newminiiso_branch->GetEntry(index);
			} else { 
				printf("branch passes_id_newminiiso_branch does not exist!\n");
				exit(1);
			}
			passes_id_newminiiso_isLoaded = true;
		}
		return passes_id_newminiiso_;
	}
	bool &	FO()
	{
		if (not FO_isLoaded) {
			if (FO_branch != 0) {
				FO_branch->GetEntry(index);
			} else { 
				printf("branch FO_branch does not exist!\n");
				exit(1);
			}
			FO_isLoaded = true;
		}
		return FO_;
	}
	bool &	FO_ptrel()
	{
		if (not FO_ptrel_isLoaded) {
			if (FO_ptrel_branch != 0) {
				FO_ptrel_branch->GetEntry(index);
			} else { 
				printf("branch FO_ptrel_branch does not exist!\n");
				exit(1);
			}
			FO_ptrel_isLoaded = true;
		}
		return FO_ptrel_;
	}
	bool &	FO_miniiso()
	{
		if (not FO_miniiso_isLoaded) {
			if (FO_miniiso_branch != 0) {
				FO_miniiso_branch->GetEntry(index);
			} else { 
				printf("branch FO_miniiso_branch does not exist!\n");
				exit(1);
			}
			FO_miniiso_isLoaded = true;
		}
		return FO_miniiso_;
	}
	bool &	FO_newminiiso()
	{
		if (not FO_newminiiso_isLoaded) {
			if (FO_newminiiso_branch != 0) {
				FO_newminiiso_branch->GetEntry(index);
			} else { 
				printf("branch FO_newminiiso_branch does not exist!\n");
				exit(1);
			}
			FO_newminiiso_isLoaded = true;
		}
		return FO_newminiiso_;
	}
	bool &	FO_NoIso()
	{
		if (not FO_NoIso_isLoaded) {
			if (FO_NoIso_branch != 0) {
				FO_NoIso_branch->GetEntry(index);
			} else { 
				printf("branch FO_NoIso_branch does not exist!\n");
				exit(1);
			}
			FO_NoIso_isLoaded = true;
		}
		return FO_NoIso_;
	}
	float &ip3d()
	{
		if (not ip3d_isLoaded) {
			if (ip3d_branch != 0) {
				ip3d_branch->GetEntry(index);
			} else { 
				printf("branch ip3d_branch does not exist!\n");
				exit(1);
			}
			ip3d_isLoaded = true;
		}
		return ip3d_;
	}
	float &ip3derr()
	{
		if (not ip3derr_isLoaded) {
			if (ip3derr_branch != 0) {
				ip3derr_branch->GetEntry(index);
			} else { 
				printf("branch ip3derr_branch does not exist!\n");
				exit(1);
			}
			ip3derr_isLoaded = true;
		}
		return ip3derr_;
	}
	int &type()
	{
		if (not type_isLoaded) {
			if (type_branch != 0) {
				type_branch->GetEntry(index);
			} else { 
				printf("branch type_branch does not exist!\n");
				exit(1);
			}
			type_isLoaded = true;
		}
		return type_;
	}
	float &mt()
	{
		if (not mt_isLoaded) {
			if (mt_branch != 0) {
				mt_branch->GetEntry(index);
			} else { 
				printf("branch mt_branch does not exist!\n");
				exit(1);
			}
			mt_isLoaded = true;
		}
		return mt_;
	}
	float &ptrelv0()
	{
		if (not ptrelv0_isLoaded) {
			if (ptrelv0_branch != 0) {
				ptrelv0_branch->GetEntry(index);
			} else { 
				printf("branch ptrelv0_branch does not exist!\n");
				exit(1);
			}
			ptrelv0_isLoaded = true;
		}
		return ptrelv0_;
	}
	float &ptrelv1()
	{
		if (not ptrelv1_isLoaded) {
			if (ptrelv1_branch != 0) {
				ptrelv1_branch->GetEntry(index);
			} else { 
				printf("branch ptrelv1_branch does not exist!\n");
				exit(1);
			}
			ptrelv1_isLoaded = true;
		}
		return ptrelv1_;
	}
	float &miniiso()
	{
		if (not miniiso_isLoaded) {
			if (miniiso_branch != 0) {
				miniiso_branch->GetEntry(index);
			} else { 
				printf("branch miniiso_branch does not exist!\n");
				exit(1);
			}
			miniiso_isLoaded = true;
		}
		return miniiso_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &jet_close_lep()
	{
		if (not jet_close_lep_isLoaded) {
			if (jet_close_lep_branch != 0) {
				jet_close_lep_branch->GetEntry(index);
			} else { 
				printf("branch jet_close_lep_branch does not exist!\n");
				exit(1);
			}
			jet_close_lep_isLoaded = true;
		}
		return *jet_close_lep_;
	}
	int &tag_charge()
	{
		if (not tag_charge_isLoaded) {
			if (tag_charge_branch != 0) {
				tag_charge_branch->GetEntry(index);
			} else { 
				printf("branch tag_charge_branch does not exist!\n");
				exit(1);
			}
			tag_charge_isLoaded = true;
		}
		return tag_charge_;
	}
	float &dilep_mass()
	{
		if (not dilep_mass_isLoaded) {
			if (dilep_mass_branch != 0) {
				dilep_mass_branch->GetEntry(index);
			} else { 
				printf("branch dilep_mass_branch does not exist!\n");
				exit(1);
			}
			dilep_mass_isLoaded = true;
		}
		return dilep_mass_;
	}
	float &sigmaIEtaIEta_full5x5()
	{
		if (not sigmaIEtaIEta_full5x5_isLoaded) {
			if (sigmaIEtaIEta_full5x5_branch != 0) {
				sigmaIEtaIEta_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch sigmaIEtaIEta_full5x5_branch does not exist!\n");
				exit(1);
			}
			sigmaIEtaIEta_full5x5_isLoaded = true;
		}
		return sigmaIEtaIEta_full5x5_;
	}
	float &etaSC()
	{
		if (not etaSC_isLoaded) {
			if (etaSC_branch != 0) {
				etaSC_branch->GetEntry(index);
			} else { 
				printf("branch etaSC_branch does not exist!\n");
				exit(1);
			}
			etaSC_isLoaded = true;
		}
		return etaSC_;
	}
	float &dEtaIn()
	{
		if (not dEtaIn_isLoaded) {
			if (dEtaIn_branch != 0) {
				dEtaIn_branch->GetEntry(index);
			} else { 
				printf("branch dEtaIn_branch does not exist!\n");
				exit(1);
			}
			dEtaIn_isLoaded = true;
		}
		return dEtaIn_;
	}
	float &dPhiIn()
	{
		if (not dPhiIn_isLoaded) {
			if (dPhiIn_branch != 0) {
				dPhiIn_branch->GetEntry(index);
			} else { 
				printf("branch dPhiIn_branch does not exist!\n");
				exit(1);
			}
			dPhiIn_isLoaded = true;
		}
		return dPhiIn_;
	}
	float &hOverE()
	{
		if (not hOverE_isLoaded) {
			if (hOverE_branch != 0) {
				hOverE_branch->GetEntry(index);
			} else { 
				printf("branch hOverE_branch does not exist!\n");
				exit(1);
			}
			hOverE_isLoaded = true;
		}
		return hOverE_;
	}
	float &ecalEnergy()
	{
		if (not ecalEnergy_isLoaded) {
			if (ecalEnergy_branch != 0) {
				ecalEnergy_branch->GetEntry(index);
			} else { 
				printf("branch ecalEnergy_branch does not exist!\n");
				exit(1);
			}
			ecalEnergy_isLoaded = true;
		}
		return ecalEnergy_;
	}
	float &eOverPIn()
	{
		if (not eOverPIn_isLoaded) {
			if (eOverPIn_branch != 0) {
				eOverPIn_branch->GetEntry(index);
			} else { 
				printf("branch eOverPIn_branch does not exist!\n");
				exit(1);
			}
			eOverPIn_isLoaded = true;
		}
		return eOverPIn_;
	}
	bool &	conv_vtx_flag()
	{
		if (not conv_vtx_flag_isLoaded) {
			if (conv_vtx_flag_branch != 0) {
				conv_vtx_flag_branch->GetEntry(index);
			} else { 
				printf("branch conv_vtx_flag_branch does not exist!\n");
				exit(1);
			}
			conv_vtx_flag_isLoaded = true;
		}
		return conv_vtx_flag_;
	}
	int &exp_innerlayers()
	{
		if (not exp_innerlayers_isLoaded) {
			if (exp_innerlayers_branch != 0) {
				exp_innerlayers_branch->GetEntry(index);
			} else { 
				printf("branch exp_innerlayers_branch does not exist!\n");
				exit(1);
			}
			exp_innerlayers_isLoaded = true;
		}
		return exp_innerlayers_;
	}
	int &charge()
	{
		if (not charge_isLoaded) {
			if (charge_branch != 0) {
				charge_branch->GetEntry(index);
			} else { 
				printf("branch charge_branch does not exist!\n");
				exit(1);
			}
			charge_isLoaded = true;
		}
		return charge_;
	}
	int &sccharge()
	{
		if (not sccharge_isLoaded) {
			if (sccharge_branch != 0) {
				sccharge_branch->GetEntry(index);
			} else { 
				printf("branch sccharge_branch does not exist!\n");
				exit(1);
			}
			sccharge_isLoaded = true;
		}
		return sccharge_;
	}
	int &ckf_charge()
	{
		if (not ckf_charge_isLoaded) {
			if (ckf_charge_branch != 0) {
				ckf_charge_branch->GetEntry(index);
			} else { 
				printf("branch ckf_charge_branch does not exist!\n");
				exit(1);
			}
			ckf_charge_isLoaded = true;
		}
		return ckf_charge_;
	}
	bool &	threeChargeAgree()
	{
		if (not threeChargeAgree_isLoaded) {
			if (threeChargeAgree_branch != 0) {
				threeChargeAgree_branch->GetEntry(index);
			} else { 
				printf("branch threeChargeAgree_branch does not exist!\n");
				exit(1);
			}
			threeChargeAgree_isLoaded = true;
		}
		return threeChargeAgree_;
	}
	int &pid_PFMuon()
	{
		if (not pid_PFMuon_isLoaded) {
			if (pid_PFMuon_branch != 0) {
				pid_PFMuon_branch->GetEntry(index);
			} else { 
				printf("branch pid_PFMuon_branch does not exist!\n");
				exit(1);
			}
			pid_PFMuon_isLoaded = true;
		}
		return pid_PFMuon_;
	}
	float &gfit_chi2()
	{
		if (not gfit_chi2_isLoaded) {
			if (gfit_chi2_branch != 0) {
				gfit_chi2_branch->GetEntry(index);
			} else { 
				printf("branch gfit_chi2_branch does not exist!\n");
				exit(1);
			}
			gfit_chi2_isLoaded = true;
		}
		return gfit_chi2_;
	}
	float &gfit_ndof()
	{
		if (not gfit_ndof_isLoaded) {
			if (gfit_ndof_branch != 0) {
				gfit_ndof_branch->GetEntry(index);
			} else { 
				printf("branch gfit_ndof_branch does not exist!\n");
				exit(1);
			}
			gfit_ndof_isLoaded = true;
		}
		return gfit_ndof_;
	}
	int &gfit_validSTAHits()
	{
		if (not gfit_validSTAHits_isLoaded) {
			if (gfit_validSTAHits_branch != 0) {
				gfit_validSTAHits_branch->GetEntry(index);
			} else { 
				printf("branch gfit_validSTAHits_branch does not exist!\n");
				exit(1);
			}
			gfit_validSTAHits_isLoaded = true;
		}
		return gfit_validSTAHits_;
	}
	int &numberOfMatchedStations()
	{
		if (not numberOfMatchedStations_isLoaded) {
			if (numberOfMatchedStations_branch != 0) {
				numberOfMatchedStations_branch->GetEntry(index);
			} else { 
				printf("branch numberOfMatchedStations_branch does not exist!\n");
				exit(1);
			}
			numberOfMatchedStations_isLoaded = true;
		}
		return numberOfMatchedStations_;
	}
	int &validPixelHits()
	{
		if (not validPixelHits_isLoaded) {
			if (validPixelHits_branch != 0) {
				validPixelHits_branch->GetEntry(index);
			} else { 
				printf("branch validPixelHits_branch does not exist!\n");
				exit(1);
			}
			validPixelHits_isLoaded = true;
		}
		return validPixelHits_;
	}
	int &nlayers()
	{
		if (not nlayers_isLoaded) {
			if (nlayers_branch != 0) {
				nlayers_branch->GetEntry(index);
			} else { 
				printf("branch nlayers_branch does not exist!\n");
				exit(1);
			}
			nlayers_isLoaded = true;
		}
		return nlayers_;
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern LeptonTree lepton_tree_obj;
#endif

namespace lepton_tree {
	const float &evt_pfmet();
	const float &evt_pfmetPhi();
	const int &evt_event();
	const int &evt_lumiBlock();
	const int &evt_run();
	const bool &filt_csc();
	const bool &filt_hbhe();
	const bool &filt_hcallaser();
	const bool &filt_ecaltp();
	const bool &filt_trkfail();
	const bool &filt_eebadsc();
	const bool &evt_isRealData();
	const float &scale1fb();
	const float &evt_xsec_incl();
	const float &evt_kfactor();
	const float &gen_met();
	const float &gen_metPhi();
	const float &njets();
	const float &ht();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets();
	const vector<float> &jets_disc();
	const TString &sample();
	const int &nFOs();
	const int &nvtx();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &tag_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &dilep_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mc_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mc_motherp4();
	const int &id();
	const int &idx();
	const float &dxyPV();
	const float &dZ();
	const float &dxyPV_err();
	const int &motherID();
	const int &mc_id();
	const float &iso();
	const bool &passes_id();
	const bool &passes_id_ptrel();
	const bool &passes_id_miniiso();
	const bool &passes_id_newminiiso();
	const bool &FO();
	const bool &FO_ptrel();
	const bool &FO_miniiso();
	const bool &FO_newminiiso();
	const bool &FO_NoIso();
	const float &ip3d();
	const float &ip3derr();
	const int &type();
	const float &mt();
	const float &ptrelv0();
	const float &ptrelv1();
	const float &miniiso();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &jet_close_lep();
	const int &tag_charge();
	const float &dilep_mass();
	const float &sigmaIEtaIEta_full5x5();
	const float &etaSC();
	const float &dEtaIn();
	const float &dPhiIn();
	const float &hOverE();
	const float &ecalEnergy();
	const float &eOverPIn();
	const bool &conv_vtx_flag();
	const int &exp_innerlayers();
	const int &charge();
	const int &sccharge();
	const int &ckf_charge();
	const bool &threeChargeAgree();
	const int &pid_PFMuon();
	const float &gfit_chi2();
	const float &gfit_ndof();
	const int &gfit_validSTAHits();
	const int &numberOfMatchedStations();
	const int &validPixelHits();
	const int &nlayers();
}
#endif
