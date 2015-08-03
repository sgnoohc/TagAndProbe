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
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOG"      )) {return Selection::EGammaMediumPOG;      } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumPOGIso"   )) {return Selection::EGammaMediumPOGIso;   } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumSTOP"     )) {return Selection::EGammaMediumSTOP;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumSTOPIso"  )) {return Selection::EGammaMediumSTOPIso;  } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaVetoHAD"        )) {return Selection::EGammaVetoHAD;        } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaVetoHADIso"     )) {return Selection::EGammaVetoHADIso;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenID"       )) {return Selection::MuTightWPDenID;       } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenIso"      )) {return Selection::MuTightWPDenIso;      } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenBoth"     )) {return Selection::MuTightWPDenBoth;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPNum"         )) {return Selection::MuTightWPNum;         } 
        throw std::invalid_argument("[tnp::GetSelectionFromString]: ERROR - invalid value!"); 
    }

    std::string GetStringFromSelection(const Selection::value_type sel_type)
    {
        if (sel_type == Selection::EGammaGsfElectron    ) return "EGammaGsfElectron";
        if (sel_type == Selection::EGammaMediumPOG      ) return "EGammaMediumPOG";
        if (sel_type == Selection::EGammaMediumPOGIso   ) return "EGammaMediumPOGIso";
        if (sel_type == Selection::EGammaMediumSTOP     ) return "EGammaMediumSTOP";
        if (sel_type == Selection::EGammaMediumSTOPIso  ) return "EGammaMediumSTOPIso";
        if (sel_type == Selection::EGammaVetoHAD        ) return "EGammaVetoHAD";
        if (sel_type == Selection::EGammaVetoHADIso     ) return "EGammaVetoHADIso";
        if (sel_type == Selection::MuTightWPDenID       ) return "MuTightWPDenID";
        if (sel_type == Selection::MuTightWPDenIso      ) return "MuTightWPDenIso";
        if (sel_type == Selection::MuTightWPDenBoth     ) return "MuTightWPDenBoth";
        if (sel_type == Selection::MuTightWPNum         ) return "MuTightWPNum";
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
            const bool el_is_crack    = fabs(etaSC()) < 1.566 && fabs(etaSC()) > 1.442;
            const float el_tag_pt      = tag_p4().pt(); 
            const float el_tag_pt_cut  = 30.0;

            // cut decisions 
            const bool el_passes_pt       = (el_tag_pt > el_tag_pt_cut);
            const bool el_passes_trig_tag = evt_isRealData() ? tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0: true; //true; //GZ need update
	    const bool el_GsfElectron_den = ( !el_is_crack && el_passes_pt && el_passes_trig_tag );

	    // implement by hand STOP_medium_v2
	    bool el_passes_STOP_medium_v2_noiso = true;
	    if (el_is_barrel) {
	      if (sigmaIEtaIEta_full5x5()   >= 0.009996) el_passes_STOP_medium_v2_noiso = false;// full5x5_sigmaIetaIeta   
	      if (fabs(dEtaIn()    )  >= 0.008925) el_passes_STOP_medium_v2_noiso = false;// abs(dEtaIn)                     
	      if (fabs(dPhiIn()    )  >= 0.035973) el_passes_STOP_medium_v2_noiso = false;// abs(dPhiIn)                     
	      if (hOverE()            >= 0.050537) el_passes_STOP_medium_v2_noiso = false;// hOverE                                  
	      if (fabs( (1.0/ecalEnergy()     ) -  (eOverPIn()/ecalEnergy())) >= 0.091942) el_passes_STOP_medium_v2_noiso = false;// ooEmooP    
	      if (fabs(dxyPV())                       >= 0.012235) el_passes_STOP_medium_v2_noiso = false;// abs(d0)                                 
	      if (fabs(dZ())                       >= 0.042020) el_passes_STOP_medium_v2_noiso = false;// abs(dz)                                 
	      if (exp_innerlayers()  > 1        ) el_passes_STOP_medium_v2_noiso = false;// expectedMissingInnerHits 
	      if (conv_vtx_flag()        ) el_passes_STOP_medium_v2_noiso = false;// pass conversion veto    
	    }
	    else if ( fabs(etaSC()) < 2.5 ) {
	      if (sigmaIEtaIEta_full5x5()                   >= 0.030135) el_passes_STOP_medium_v2_noiso = false;// full5x5_sigmaIetaIeta   
	      if (fabs(dEtaIn()                          )  >= 0.007429) el_passes_STOP_medium_v2_noiso = false;// abs(dEtaIn)                      
	      if (fabs(dPhiIn()                          )  >= 0.067879) el_passes_STOP_medium_v2_noiso = false;// abs(dPhiIn)                      
	      if (hOverE()                                  >= 0.086782) el_passes_STOP_medium_v2_noiso = false;// hOverE                                   
	      if (fabs( (1.0/ecalEnergy()) -(eOverPIn()/ecalEnergy())) >= 0.100683) el_passes_STOP_medium_v2_noiso = false;// ooEmooP                       
	      if (fabs(dxyPV()      )                       >= 0.036719) el_passes_STOP_medium_v2_noiso = false;// abs(d0)                                 
	      if (fabs(dZ()       )                       >= 0.138142) el_passes_STOP_medium_v2_noiso = false;// abs(dz)                                  
	      if (exp_innerlayers()                         > 1        ) el_passes_STOP_medium_v2_noiso = false;// expectedMissingInnerHits
	      if (conv_vtx_flag()                                      ) el_passes_STOP_medium_v2_noiso = false;// pass conversion veto    
	    }
	    bool el_passes_STOP_medium_v2_iso =  miniisoDB() < 0.1;
	    bool el_passes_POG_mediumV2_iso = el_passes_STOP_medium_v2_noiso;
	    if (el_is_barrel && RelIso03EA() > 0.107587)  el_passes_POG_mediumV2_iso = false;
	    if (!el_is_barrel && RelIso03EA() > 0.113254)  el_passes_POG_mediumV2_iso = false;
	    
	    if (selection == Selection::EGammaGsfElectron    ) {
                if (not el_GsfElectron_den)       {return false;}
	    }
	    if (selection == Selection::EGammaMediumPOG      ) {
                if (not el_GsfElectron_den)       {return false;}
                if (not el_passes_STOP_medium_v2_noiso)    {return false;}
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
            const float mu_tag_pt_cut  = 30.0;

            // cut decisions 
            const bool mu_passes_pt       = (mu_tag_pt > mu_tag_pt_cut);
            const bool mu_passes_trig_tag = true; //GZ need update
            const bool mu_passes_pog_iso  = (mu_iso < mu_iso_pog_cut); 
            const bool mu_passes_pog_id   = passes_SS_tight_noiso_v3(); 

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
 
