#include "TagAndProbe/Analysis/interface/Measurement.h"

#include <iostream>

using namespace tnp;
using namespace std;

// Tools
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "TagAndProbe/Analysis/interface/LeptonTree.h"
#include "TagAndProbe/Analysis/interface/LeptonSelections.h"

namespace tnp
{
    Lepton::value_type GetLeptonFromString(const std::string& lepton_name)
    {
        if (lt::string_lower(lepton_name) == "muon"    ) {return Lepton::Muon;    }
        if (lt::string_lower(lepton_name) == "electron") {return Lepton::Electron;}
        throw std::invalid_argument("[tnp::GetLeptonFromString]: ERROR - invalid value!"); 
    }

    std::string GetStringFromLepton(const Lepton::value_type lepton_type)
    {
        if (lepton_type == Lepton::Muon    ) return "muon";
        if (lepton_type == Lepton::Electron) return "electron";
        throw std::invalid_argument(Form("[tnp::GetStringFromLepton]: ERROR - invalid value! %u", lepton_type)); 
    }

    Selection::value_type GetSelectionFromString(const std::string& sel_name)
    {
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaGsfElectron"    )) {return Selection::EGammaGsfElectron;    } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaPFElectron"     )) {return Selection::EGammaPFElectron;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaPFElectronChIso")) {return Selection::EGammaPFElectronChIso;} 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOG"      )) {return Selection::EGammaMediumPOG;      } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOGnoDEta")) {return Selection::EGammaMediumPOGnoDEta;} 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOGnoIP"  )) {return Selection::EGammaMediumPOGnoIP;  } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOGnoDPhi")) {return Selection::EGammaMediumPOGnoDPhi;} 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOGnoSieie")) {return Selection::EGammaMediumPOGnoSieie;} 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOGnoConv")) {return Selection::EGammaMediumPOGnoConv;} 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOGIso"   )) {return Selection::EGammaMediumPOGIso;   } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumSTOP"     )) {return Selection::EGammaMediumSTOP;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumSTOPIso"  )) {return Selection::EGammaMediumSTOPIso;  } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaVetoHAD"        )) {return Selection::EGammaVetoHAD;        } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaVetoHADIso"     )) {return Selection::EGammaVetoHADIso;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaTightSS"        )) {return Selection::EGammaTightSS;        }  
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaTightSSIso"     )) {return Selection::EGammaTightSSIso;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaTightSSMVA"     )) {return Selection::EGammaTightSSMVA;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaTightSSMVA3ch"  )) {return Selection::EGammaTightSSMVA3ch;  } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaTightSSMVA3chIP")) {return Selection::EGammaTightSSMVA3chIP;} 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaTightSSMVA3chIPconv")) {return Selection::EGammaTightSSMVA3chIPconv;} 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenID"       )) {return Selection::MuTightWPDenID;       } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenIso"      )) {return Selection::MuTightWPDenIso;      } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenBoth"     )) {return Selection::MuTightWPDenBoth;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPNum"         )) {return Selection::MuTightWPNum;         } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuPFDen"              )) {return Selection::MuPFDen;              } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuPFChIso"            )) {return Selection::MuPFChIso;            } 
        throw std::invalid_argument("[tnp::GetSelectionFromString]: ERROR - invalid value!"); 
    }

    std::string GetStringFromSelection(const Selection::value_type sel_type)
    {
        if (sel_type == Selection::EGammaGsfElectron    ) return "EGammaGsfElectron";
        if (sel_type == Selection::EGammaPFElectron     ) return "EGammaPFElectron";
        if (sel_type == Selection::EGammaPFElectronChIso) return "EGammaPFElectronChIso";
        if (sel_type == Selection::EGammaMediumPOG      ) return "EGammaMediumPOG";
        if (sel_type == Selection::EGammaMediumPOGnoDEta) return "EGammaMediumPOGnoDEta";
        if (sel_type == Selection::EGammaMediumPOGnoIP  ) return "EGammaMediumPOGnoIP";
        if (sel_type == Selection::EGammaMediumPOGnoDPhi) return "EGammaMediumPOGnoDPhi";
        if (sel_type == Selection::EGammaMediumPOGnoSieie) return "EGammaMediumPOGnoSieie";
        if (sel_type == Selection::EGammaMediumPOGnoConv) return "EGammaMediumPOGnoConv";
        if (sel_type == Selection::EGammaMediumPOGIso   ) return "EGammaMediumPOGIso";
        if (sel_type == Selection::EGammaMediumSTOP     ) return "EGammaMediumSTOP";
        if (sel_type == Selection::EGammaMediumSTOPIso  ) return "EGammaMediumSTOPIso";
        if (sel_type == Selection::EGammaVetoHAD        ) return "EGammaVetoHAD";
        if (sel_type == Selection::EGammaVetoHADIso     ) return "EGammaVetoHADIso";
        if (sel_type == Selection::EGammaTightSS        ) return "EGammaTightSS";
        if (sel_type == Selection::EGammaTightSSIso     ) return "EGammaTightSSIso";
        if (sel_type == Selection::EGammaTightSSMVA     ) return "EGammaTightSSMVA";
        if (sel_type == Selection::EGammaTightSSMVA3ch  ) return "EGammaTightSSMVA3ch";
        if (sel_type == Selection::EGammaTightSSMVA3chIP) return "EGammaTightSSMVA3chIP";
        if (sel_type == Selection::EGammaTightSSMVA3chIPconv) return "EGammaTightSSMVA3chIPconv";
        if (sel_type == Selection::MuTightWPDenID       ) return "MuTightWPDenID";
        if (sel_type == Selection::MuTightWPDenIso      ) return "MuTightWPDenIso";
        if (sel_type == Selection::MuTightWPDenBoth     ) return "MuTightWPDenBoth";
        if (sel_type == Selection::MuTightWPNum         ) return "MuTightWPNum";
        if (sel_type == Selection::MuPFDen              ) return "MuPFDen";
        if (sel_type == Selection::MuPFChIso            ) return "MuPFChIso";
        throw std::invalid_argument("[tnp::GetStringFromSelection]: ERROR - invalid value!"); 
    }

    // passes selection based on above enum
    bool PassesSelection(const Lepton::value_type lepton_type, const Selection::value_type selection, const bool is_data)
    {
        using namespace lepton_tree;

        // --------------------------------------------------------------------------- //
        // electrons
        // --------------------------------------------------------------------------- //

        if (lepton_type == Lepton::Electron)
        {
            // cut values and variables
	    const bool el_is_barrel   = fabs(etaSC()) < 1.479;
            //const bool el_is_endcap   = fabs(etaSC()) > 1.566 && fabs(etaSC()) < 2.5;
            //const bool el_is_crack    = fabs(etaSC()) < 1.566 && fabs(etaSC()) > 1.442;
            const float el_tag_pt      = tag_p4().pt(); 
            const float el_tag_pt_cut  = 30.0;

            // cut decisions 
            const bool el_passes_tag_pt    = (el_tag_pt > el_tag_pt_cut);
            const bool el_passes_tag_eta   = (fabs(tag_p4().eta()) < 2.1); 
            const bool el_passes_tag_trig  = evt_isRealData() ? tag_HLT_Ele22_eta2p1_WPLoose_Gsf() > 0: true; //true; //GZ need update
	    //	    const bool el_GsfElectron_den = ( !el_is_crack && el_passes_pt && el_passes_trig_tag );
	    const bool el_GsfElectron_den = ( el_passes_tag_pt && el_passes_tag_trig && el_passes_tag_eta );
	    const bool el_PFElectron_den  = (isPF() && fabs(dZ()) < 0.1); // additional requirement included in denominator definition. 
	    const bool el_passes_chIso = (AbsTrkIso() / p4().pt() < 0.2);

	    // implement by hand STOP_medium_v2
	    bool passSieie = true, passDeta = true, passDphi = true, passHoverE = true, passOemoop = true, passIP = true, passConv = true, passIso = true;
	    if (el_is_barrel) {
	      if (sigmaIEtaIEta_full5x5()   >= 0.009996) passSieie = false;
	      if (fabs(dEtaIn()    )  >= 0.008925) passDeta = false;
	      if (fabs(dPhiIn()    )  >= 0.035973) passDphi = false;
	      if (hOverE()            >= 0.050537) passHoverE = false;
	      if (fabs( (1.0/ecalEnergy()     ) -  (eOverPIn()/ecalEnergy())) >= 0.091942) passOemoop = false;
	      if (fabs(dxyPV())                       >= 0.012235) passIP = false;
	      if (fabs(dZ())                       >= 0.042020) passIP = false;
	      if (exp_innerlayers()  > 1        ) passConv = false;
	      if (conv_vtx_flag()        ) passConv = false;
	      if (RelIso03EA() > 0.107587) passIso = false;
	    }
	    else if ( fabs(etaSC()) < 2.5 ) {
	      if (sigmaIEtaIEta_full5x5()                   >= 0.030135) passSieie = false;
	      if (fabs(dEtaIn()                          )  >= 0.007429) passDeta = false;
	      if (fabs(dPhiIn()                          )  >= 0.067879) passDphi = false;
	      if (hOverE()                                  >= 0.086782) passHoverE = false;
	      if (fabs( (1.0/ecalEnergy()) -(eOverPIn()/ecalEnergy())) >= 0.100683) passOemoop = false;
	      if (fabs(dxyPV()      )                       >= 0.036719) passIP = false;
	      if (fabs(dZ()       )                       >= 0.138142) passIP = false;
	      if (exp_innerlayers()                         > 1        ) passConv = false;
	      if (conv_vtx_flag()                                      ) passConv = false;
	      if (RelIso03EA() > 0.113254 ) passIso = false;
	    }
	    bool el_passes_STOP_medium_v2_noiso = (passSieie && passDeta && passDphi && passHoverE && passOemoop && passIP && passConv);
	    bool el_passes_STOP_medium_v2_noiso_noDEta  = (passSieie && passDphi && passHoverE && passOemoop && passIP && passConv);
	    bool el_passes_STOP_medium_v2_noiso_noIP    = (passSieie && passDeta && passDphi && passHoverE && passOemoop && passConv);
	    bool el_passes_STOP_medium_v2_noiso_noDPhi  = (passSieie && passDeta && passHoverE && passOemoop && passIP && passConv);
	    bool el_passes_STOP_medium_v2_noiso_noSieie = (passDeta && passDphi && passHoverE && passOemoop && passIP && passConv);
	    bool el_passes_STOP_medium_v2_noiso_noConv  = (passSieie && passDeta && passDphi && passHoverE && passOemoop && passIP );
	    bool el_passes_STOP_medium_v2_noDEta  = el_passes_STOP_medium_v2_noiso_noDEta  && passIso;
	    bool el_passes_STOP_medium_v2_noIP    = el_passes_STOP_medium_v2_noiso_noIP    && passIso;
	    bool el_passes_STOP_medium_v2_noDPhi  = el_passes_STOP_medium_v2_noiso_noDPhi  && passIso;
	    bool el_passes_STOP_medium_v2_noSieie = el_passes_STOP_medium_v2_noiso_noSieie && passIso;
	    bool el_passes_STOP_medium_v2_noConv  = el_passes_STOP_medium_v2_noiso_noConv  && passIso;

	    bool el_passes_STOP_medium_v2_iso =  el_passes_STOP_medium_v2_noiso && miniisoDB() < 0.1;
	    bool el_passes_POG_mediumV2_iso = el_passes_STOP_medium_v2_noiso && passIso;

	    float MVA = mva();
	    bool el_passes_SS_tightMVA = false;
	    if (fabs(etaSC()) < 0.8 &&  MVA > 0.73) el_passes_SS_tightMVA = true;
	    if (fabs(etaSC()) >= 0.8 && fabs(etaSC()) <= 1.479 && MVA > 0.57) el_passes_SS_tightMVA = true;
	    if (fabs(etaSC()) > 1.479 && MVA > 0.05) el_passes_SS_tightMVA = true;
	    

	    bool el_passes_threeChargeAgree = threeChargeAgree();
	    bool el_passes_SS_tightIP = true;
	    if (fabs(dZ()) >= 0.1) el_passes_SS_tightIP = false;
	    if (fabs(ip3d())/ip3derr() >= 4) el_passes_SS_tightIP = false;
	    bool el_passes_SS_conv = true;
	    if (conv_vtx_flag()) el_passes_SS_conv = false;
	    if (exp_innerlayers() > 0) el_passes_SS_conv = false;
	    
	    if (selection == Selection::EGammaGsfElectron    ) {
                if (not el_GsfElectron_den)       {return false;}
	    }
	    if (selection == Selection::EGammaPFElectron     ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_PFElectron_den)        {return false;}
	    }
	    if (selection == Selection::EGammaPFElectronChIso) {
	      if (not el_GsfElectron_den)         {return false;}
	      if (not el_PFElectron_den)          {return false;}
	      if (not el_passes_chIso)            {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOG      ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noiso)    {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOGnoDEta      ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noDEta)    {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOGnoIP      ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noIP)    {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOGnoDPhi      ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noDPhi)    {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOGnoSieie      ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noSieie)    {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOGnoConv      ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noConv)    {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOGIso   ) {
                if (not el_GsfElectron_den)       {return false;}
		if (not el_passes_POG_mediumV2_iso)    {return false;}
	    }
	    if (selection == Selection::EGammaMediumSTOP     ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noiso) {return false;}
	    }
	    if (selection == Selection::EGammaMediumSTOPIso  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_iso) {return false;}
	    }
	    if (selection == Selection::EGammaVetoHAD  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not passes_HAD_veto_noiso_v3()) {return false;}
	    }
	    if (selection == Selection::EGammaVetoHADIso  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not passes_HAD_veto_v3())     {return false;}
	    }	    
	    if (selection == Selection::EGammaTightSS  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not passes_SS_tight_noiso_v3()) {return false;}
	    }
	    if (selection == Selection::EGammaTightSSIso  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not passes_SS_tight_v3())     {return false;}
	    }	  
	    if (selection == Selection::EGammaTightSSMVA  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_SS_tightMVA) {return false;}
	    }  
	    if (selection == Selection::EGammaTightSSMVA3ch  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_SS_tightMVA) {return false;}
                if (not el_passes_threeChargeAgree) {return false;}
	    }  
	    if (selection == Selection::EGammaTightSSMVA3chIP  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_SS_tightMVA) {return false;}
                if (not el_passes_threeChargeAgree) {return false;}
                if (not el_passes_SS_tightIP) {return false;}
	    } 
	    if (selection == Selection::EGammaTightSSMVA3chIPconv  ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_SS_tightMVA) {return false;}
                if (not el_passes_threeChargeAgree) {return false;}
                if (not el_passes_SS_tightIP) {return false;}
                if (not el_passes_SS_conv) {return false;}
	    } 


        }

        // --------------------------------------------------------------------------- //
        // muons
        // --------------------------------------------------------------------------- //

        if (lepton_type == Lepton::Muon)
        {
            // cut values and variables
            const float mu_tag_pt      = tag_p4().pt();
            const float mu_iso         = miniiso();
            const float mu_iso_pog_cut = 0.15;  
            const float mu_tag_pt_cut  = 20.0;

            // cut decisions 
            const bool mu_passes_pt       = (mu_tag_pt > mu_tag_pt_cut);
            const bool mu_passes_trig_tag = (tag_HLT_IsoMu20() > 0) || (tag_HLT_IsoTkMu20() > 0); 
            const bool mu_passes_pog_iso  = (mu_iso < mu_iso_pog_cut); 
            const bool mu_passes_pog_id   = passes_SS_tight_noiso_v3(); 
            const bool mu_passes_PFChIso  = (AbsTrkIso() / p4().pt() < 0.2); //RelIso03() < 0.2; //
            const bool mu_passes_PF       = isPF(); 
            // Muon POG Selections (2012)
            // From: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
            // --------------------------------------------------------------------------- //

            // Isolation
            if (selection == Selection::MuTightWPDenIso)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
                if (not mu_passes_pog_id)   {return false;}
            }
	    
            // ID
            if (selection == Selection::MuTightWPDenID)
            {
                if (not mu_passes_pt)        {return false;}
                if (not mu_passes_trig_tag)  {return false;}
                if (not mu_passes_pog_iso)   {return false;}
            }

            // Both ID and isolation
            if (selection == Selection::MuTightWPDenBoth)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
            }

            // Numerator
            if (selection == Selection::MuTightWPNum)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
                if (not mu_passes_pog_id)   {return false;}
                if (not mu_passes_pog_iso)  {return false;}
            }
	    if (selection == Selection::MuPFDen)
	    {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_PF)       {return false;}
                if (not mu_passes_trig_tag) {return false;}	      
	    }
	    if (selection == Selection::MuPFChIso)
	    {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_PF)       {return false;}
                if (not mu_passes_trig_tag) {return false;}	     
                if (not mu_passes_PFChIso)  {return false;}	      
	    }
            
        }

        // other values are invalid
        if (lepton_type == Lepton::static_size)
        {
            std::cout << "lepton_type is invalid" << std::endl;
            return false;
        }

        // if we got here -- it's selected
        return true;
    }

} // namespace tnp
 
