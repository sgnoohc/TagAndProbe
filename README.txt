# TagAndProbe

To set things up on UAF:

   export SCRAM_ARCH=slc6_amd64_gcc481
   source /cvmfs/cms.cern.ch/cmsset_default.sh
   cmsrel CMSSW_7_2_0
   cd CMSSW_7_2_0/src
   cmsenv

   git clone https://github.com/kelleyrw/AnalysisTools.git 
   cd $CMSSW_BASE/src/AnalysisTools
   git checkout tnp_V00-00-00
   cd $CMSSW_BASE/src/

   git clone git@github.com:gzevi/TagAndProbe.git

   scramv1 b -j10

To run the first steps (make dilepMass plots and fit them to 
extrac efficiencies) on a small LeptonTree:
   cd TagAndProbe/Analysis
   tnp_make_plots config/ElectronID_2015test.py
   tnp_eff_plots  config/ElectronID_2015test.py

-------
For detailed instructions, look at Ryan's twiki, skipping the babymaking 
http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/TagAndProbe#Calculate_the_Efficiency
-------

Regarding LeptonTrees:

Useful tool to remake LeptonTree.h/cc files after changing the LeptonTree
    TagAndProbe/Analysis/tools/makeLeptonTreeClassFiles.sh

LeptonTree making is currently based on the SSAnalysis.
To install (in a different directory), follow:
   https://github.com/cmstas/SSAnalysis/blob/master/install.sh
And to get the changes for TagAndProbe
   git checkout gzeviTNP
To make trees, all the code is in FakeRate/babymaker
