# NPS2023 Data Production
## Introduction
This is the script for clustering, photon reconstruction, and analysis tree output after waveform analysis and energy calibration for NPS. \
The methods are based on the NPS software package (https://github.com/hhuang-hep/NPS_SOFT.git), which is also used for photon reconstruction in simulation.

## Clustering and photon reconstruction
1. Input pulses 
    - TCaloBlock::AddPulse(Float_t energy, Float_t time)
    - energy: calibration coefficients*amplitudes from waveform fitting
    - time: offset-corrected timing from waveform fitting

2. Selection of pulses for clustering 
    - TCaloEvent::TriggerSim(Double_t FourBlocksThreshold)
    - Select pulses based on a threshold for every 2x2 group of blocks.
    - A pulse will be selected once it belong to a group of energy > threshold
    - Currently (1/9/2026), this threshold is set to 0.2 GeV. Can be set different according to different kinematics.

3. Clustering method and algorithm
    - TCaloEvent::DoClustering(Double_t timemin, Double_t timemax, Float_t BlockThreshold)
    - Pulses with their timing in (timemin, timemax) and energy > BlockThreshold will be considered
    - We are using pulse time in (-3, 3) and BlockThreshold = 0 (no threshold for individual pulse)
    - Clustering algorithm: cellular automaton (Nuclear Instruments and Methods in Physics Research A 362, 478-486, 1995)

4. Photon reconstruction
    - TDVCSEvent::GetPhoton(Int_t clus, Float_t a, Float_t b, Float_t &a_out, Float_t &x_corr, Float_t &y_corr)
    - This function takes information from cluster "clus", corrects cluster positions with shower depth "a", then returns the 4-momenta of reconstructed photons.
    - Parameter "b" is not inuse.
    - a_out, x_corr, and y_corr are set for output branches (see below)

## Output Trees
4 trees are output for analysis. Details of their branches are listed below.
Here is a brief discription:
1. t_prod: production tree of coincident events (i.e., pulses with timing within (-3, 3))
2. t_accdt1: production tree of accidental enents with negative pulse time within (-11, -5)
3. t_accdt2: production tree of accidental enents with positive pulse time within (5, 11)

Structure of the production trees above are the same. Information of all NPS clusters are stored

4. t_pi0sub: Monte-Carlo tree of pi0 contamination for DVCS analysis. 
    - 5000 random decays at pi0 rest fram when finding a pi0 candidate (M_pi0 +- 3 sigma) in 2 cluster events.
    - Energy cut for MC photons: 500 MeV (readout threshold during the experiment)
    - Events with 1 photon in the NPS acceptance are stored.

## Branches of production trees (t_prod, t_accdt1, t_accdt2)
### Branches from T tree
| Branch name | Type | Discription |
|:------|:------|:------|
| g.runnum | double |
| g.evnum | double |
| H.gtr.dp | double | 
| H.gtr.dp | double | 
| H.gtr.ph | double | 
| H.gtr.th | double | 
| H.gtr.px | double | 
| H.gtr.py | double | 
| H.gtr.pz | double | 
| H.react.x | double | 
| H.react.y | double | 
| H.react.z | double | 
| H.hod.beta | double | 
| H.cal.etottracknorm | double | 
| H.cer.npeSum | double |
### NPS related branches (after waveform fit + energy calibrated)
| Branch name | Type | Discription |
|:------|:------|:------|
| NPS.prod.nclust | int | Number of cluster in th event |
| NPS.prod.clusE | `vector<double>` | Cluster energy |
| NPS.prod.clusSize | `vector<int>` | Cluster size |
| NPS.prod.clusX | `vector<double>` | Cluster X position |
| NPS.prod.clusX.corr | `vector<double>` | Cluster X position after shower depth correction |
| NPS.prod.clusY | `vector<double>` | Cluster Y position |
| NPS.prod.clusY.corr | `vector<double>` | Cluster Y position after shower depth correction |
| NPS.prod.clusZ | `vector<double>` | Distance from NPS surface to target center |
| NPS.prod.clusT | `vector<double>` | Cluster timing |
| NPS.prod.clusDepth | `vector<double>` | Shower depth along the trajectories of photons |
| NPS.prod.trk.px | `vector<double>` | Reconstructed photon px |
| NPS.prod.trk.py | `vector<double>` | Reconstructed photon py |
| NPS.prod.trk.pz | `vector<double>` | Reconstructed photon pz |
| NPS.prod.trk.ene | `vector<double>` | Reconstructed photon energy |
- Note: cluster x and y are in the coordinate of NPS surface
### Branches for quick analysis/checking things
| Branch name | Type | Discription |
|:------|:------|:------|
| NPS.prod.M | double | Invariant mass only when there are 2 clusters |
| NPS.prod.Mx2 | double | Missing mass square only when there 1 or 2 cluster(s) |

## Branches of pi0 contamination (t_pi0sub) for DVCS analysis
### Branches from T tree

| Branch name | Type | Discription |
|:------|:------|:------|
| g.runnum | double |
| g.evnum | double |
| H.gtr.dp | double | 
| H.gtr.dp | double | 
| H.gtr.ph | double | 
| H.gtr.th | double | 
| H.gtr.px | double | 
| H.gtr.py | double | 
| H.gtr.pz | double | 
| H.react.x | double | 
| H.react.y | double | 
| H.react.z | double | 
| H.hod.beta | double | 
| H.cal.etottracknorm | double | 
| H.cer.npeSum | double |

### Pi0 comtamination from MC
| Branch name | Type | Discription |
|:------|:------|:------|
| pi0.cont.Xc | `vector<double>` | Impact position of photons on NPS surface |
| pi0.cont.Yc | `vector<double>` | Impact position of photons on NPS surface |
| pi0.cont.trk.px | `vector<double>` | Photon px |
| pi0.cont.trk.py | `vector<double>` | Photon py |
| pi0.cont.trk.pz | `vector<double>` | Photon pz |
| pi0.cont.trk.ene | `vector<double>` | Photon energy |
| pi0.cont.Mx2 | `vector<double>` | Missing mass square |
| pi0.cont.N0 | int | Counts of cases when 0 photon in NPS acceptance |
| pi0.cont.N1 | int | Counts of cases when 1 photon in NPS acceptance |
| pi0.cont.N2 | int | Counts of cases when 2 photons in NPS acceptance |
| pi0.cont.weight | double | Weight = 1/N2 |