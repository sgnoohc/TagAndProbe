#!/bin/bash

#use: ./makeWebsite.sh ElectronID_EGammaMediumWPDilepTrig DileptonTrigger

#locName = ElectronID_EGammaMediumSTOPIso/electron/EGammaGsfElectron_EGammaMediumSTOPIso/
#wwwName = DileptonTrigger

locName=$1
wwwName=$2

locDir=$CMSSW_BASE/src/TagAndProbe/Analysis/plots/$locName/

wwwDir=~/www/public_html/$wwwName/
mkdir -p $wwwDir/MC
mkdir -p $wwwDir/Data
rm $wwwDir/MC/* $wwwDir/Data/* $wwwDir/*.png

cp $locDir/compare/index.php $wwwDir/
cp $locDir/compare/index.php $wwwDir/Data
cp $locDir/compare/index.php $wwwDir/MC
cp $locDir/compare/p_eff_pt_eta* $wwwDir/
cp $locDir/compare/p_sf*.png $wwwDir/
cp $locDir/compare/eff_tables.txt $wwwDir/
   
cp $locDir/data_*_el_eff/png/p_*_pt*vs*eta* $wwwDir/Data
cp $locDir/dy_*_eff/png/p_*_pt*vs*eta* $wwwDir/MC

cd ~/www/public_html/
rm webtest.zip
zip -r webtest.zip $wwwName/*
scp webtest.zip $USER@lxplus.cern.ch:~/www/public_html/
ssh $USER@lxplus.cern.ch "cd  www/public_html/; unzip webtest.zip;"
cd -


