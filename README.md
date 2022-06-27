# NuEScatteringReco
Repository with the code I used for my reconstruction and analysis of nu-e elastic scattering events in SBND experiment.
This reposotory includes:

- NuEScatAnalyzer (NuEScatAna_module.cc + NuEScatAnaConfig.fcl + run_NuEScatAna.fcl) = Analyzer module to use with a single Nu-E interaction sample. The imput must be a .root file that contains events with GENIE + GEANT4 + detsim + reco1 simulation steps. Be careful: it is necessary not to drop the SimEnergyDeposit objects in g4 step.

- MuonAnalyzer (MuonAna_module.cc + MuonAnaConfig.fcl + run_MuonAna.fcl) = Analyzer module to use with a single muon interaction sample (this was use in the Kalman Filter direction reconstruction attempt). The imput must be a .root file that contains events with GENIE (prodsingle fhicl files) + GEANT4 + detsim + reco1 simulation steps. Be careful: it is necessary not to drop the SimEnergyDeposit objects in g4 step.

- KalmanFilter (SegmentsForKalman_MuonSample.cc + kalman4D_MuonSample.py + kalman5D_MuonSample.py) = macros for direction reconstruction using Kalman Filter with muon sample. 5D and 4D state vectors algorithms are included.

- ShowerSegmentation (TreeCreatorSegmentsAndPortions.cc)  = macros for shower segmentation and shower portions building. Angle between the segments with respect to the first one and shower portions directions are also obtained. A cut in these variables aims to retain the useful e.m. shower space-points for direction reconstruction.

- ShowerDirectionReco (ThetaRMSCut.cc + DerivEValueCut.cc + MaxEValueCut.cc) = File with the shower direction reconstruction macros (that uses ShowerSegmentation output as the input file), each one with its correpsonding point selection method and PCA.

- ShowerDirectionReco_KalmanFilter (ThetaRMSCut_and_KalmanFilterDirReco.py) = shower reconstruction macro (that uses ShowerSegmentation output as the input file) that uses Theta-RMS cut and reconstructs shower direction with Kalman Filter

