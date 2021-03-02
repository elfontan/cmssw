*******************************
Introduction to CMSSW: http://cms-sw.github.io


**Work on the uGT emulator** 
---
Focus on L1Trigger/L1TriggerGlobal

-  CMSSW version used for tests:
   ```
   cmsrel CMSSW_11_2_0_pre10
   ```
-  Instructions to produce test vectors (useful way to test the code):
   ```
   cd CMSSW_11_2_0_pre10/src/
   cmsenv
   git cms-init   
   git-cms-addpkg L1Trigger/L1TGlobal
   git checkout -b 1120-uGT-emulate-cuts-unconstrainedPt-impactParam
   cp L1Trigger/L1TGlobal/test/runGlobalFakeInputProducer.py .
   scram b -j  4
   voms-proxy-init --voms cms
   cmsRun runGlobalFakeInputProducer.py 
   
   ```
   
