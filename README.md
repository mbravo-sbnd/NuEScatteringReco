# NuEScatteringReco
Repository with the code I used for my reconstruction and analysis of nu-e elastic scattering events in SBND experiment.
This reposotory includes:

- NuEScatAnalyzer (NuEScatAna_module.cc + NuEScatAnaConfig.fcl + run_NuEScatAna.fcl) = Analyzer module to use with a single Nu-E interaction sample. The imput must be a .root file that contains events with GENIE + GEANT4 + detsim + reco1 simulation steps. Be careful: it is necessary not to drop the SimEnergyDeposit objects in g4 step.

- MuonAnalyzer (MuonAna_module.cc + MuonAnaConfig.fcl + run_MuonAna.fcl) = Analyzer module to use with a single muon interaction sample (this was use in the Kalman Filter direction reconstruction attempt). The imput must be a .root file that contains events with GENIE (prodsingle fhicl files) + GEANT4 + detsim + reco1 simulation steps. Be careful: it is necessary not to drop the SimEnergyDeposit objects in g4 step.

- ShowerSegmentation = File with macros used for shower segmentation and angle variables between the segments with respect to the first one. A cut in these variables aims to retain the useful e.m. shower space-points for direction reconstruction.

- KalmanFilter = File with macros for direction reconstruction using Kalman Filter with muon sample. 5D and 4D state vectors algorithms are included.

- ShowerDirectionReco = File with the shower direction reconstruction macros.
  -  KalmanFilter = File with 4D-state-vector-Kalman Filter reconstruction.
  -  PCA = File with PCA reconstruction.

-  MeetingSlides = Slides made for friday meetings

