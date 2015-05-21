#!/bin/bash

#use: ./makeWebsite.sh ElectronID_EGammaMediumWPDilepTrig DileptonTrigger

#locName = ElectronID_EGammaMediumWPDilepTrig
#wwwName = DileptonTrigger

locName=$1
wwwName=$2

locDir=$CMSSW_BASE/src/TagAndProbe/Analysis/plots/$locName/electron/EGammaMediumWPDenBoth_EGammaMediumWPNum/

wwwDir=~/www/public_html/$wwwName/
mkdir -p $wwwDir/MC
mkdir -p $wwwDir/Data

cp $locDir/compare/index.php $wwwDir/
cp $locDir/compare/index.php $wwwDir/Data
cp $locDir/compare/index.php $wwwDir/MC
cp $locDir/compare/p_eff_pt_eta* $wwwDir/
cp $locDir/compare/eff_tables.txt $wwwDir/
   
cp $locDir/data_*_el_eff/png/p_*_pt*vs*eta* $wwwDir/Data
cp $locDir/dy_*_eff/png/p_*_pt*vs*eta* $wwwDir/MC


