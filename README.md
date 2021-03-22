*******************************
Introduction to CMSSW: http://cms-sw.github.io


**Work on the uGT emulator** 
---
Focus on L1Trigger/L1TriggerGlobal

-  CMSSW version and branch used for tests:
   ```
   cmsrel CMSSW_11_2_0_pre10
   cd CMSSW_11_2_0_pre10/src/	
   cmsenv
   git cms-init
   git-cms-addpkg L1Trigger/L1TGlobal
   git checkout -b 1120-uGT-emulate-cuts-unconstrainedPt-impactParam	
   ```
-  Instructions to produce test vectors (useful way to test the code):
   ```
   cp L1Trigger/L1TGlobal/test/runGlobalFakeInputProducer.py .
   scram b -j  4
   voms-proxy-init --voms cms
   cmsRun runGlobalFakeInputProducer.py    
   ```
   L1 Trigger menu currently used for tests is: 
   ```
   L1Trigger/L1TGlobal/data/Luminosity/startup/L1Menu_test_mass_3_body_reduced.xml 
   ```
   It contains four algorithms for three-muon events with different conditions.
