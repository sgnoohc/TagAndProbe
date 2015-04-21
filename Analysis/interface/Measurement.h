#ifndef TNP_MEASUREMENT_H
#define TNP_MEASUREMENT_H

#include <string>

namespace tnp
{
    // simple lepton type
    struct Lepton
    {
        enum value_type 
        {
            Muon,
            Electron,
            static_size
        };
    };

    // simple selection type
    struct Selection
    {
        enum value_type 
        {
            // EGamma POG medium working point
            EGammaMediumWPDenID,
            EGammaMediumWPDenIso,
            EGammaMediumWPDenBoth,
            EGammaMediumWPNum,

            // Muon POG tight working point
            MuTightWPDenID,
            MuTightWPDenIso,
            MuTightWPDenBoth,
            MuTightWPNum,

            static_size
        };
    };

    // get the lepton lepton from a string
    Lepton::value_type GetLeptonFromString(const std::string& lepton_name);

    // get the string from the lepton type 
    std::string GetStringFromLepton(const Lepton::value_type lepton_type);

    // get the selection from a string
    Selection::value_type GetSelectionFromString(const std::string& sel_name);

    // get the string from the Selection 
    std::string GetStringFromSelection(const Selection::value_type sel_type);

    // passes selection based on above enum
    bool PassesSelection(const Lepton::value_type lepton_type, const Selection::value_type selection, const bool is_data);

} // namespace tnp

#endif //TNP_MEASUREMENT_H
