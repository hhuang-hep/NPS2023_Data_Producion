#include "/home/hhuang/workspace/group_h/analysis/MyHeader/Analysis.h"
#include "/group/nps/hhuang/software/NPS_SOFT/TDVCSDB.h"

Int_t ndecay = 5000; // number of decays for pi0 contamination simulation

void prodTree(int run_number, int iseg)
{   
    TH1::SetDefaultSumw2();

    // Kinematics dependent variables that are not in DB________________________________________
    Double_t clusTrsH = 0.2; // GeV
    Double_t pi0_sigm = 0; // GeV, for pi0 mass window cut of pi0 contamination subtraction
    // temp for x36_5_3
    if(run_number == 3728) pi0_sigm = 0.004602; 
    else if(run_number == 3729) pi0_sigm = 0.004319; 
    else if(run_number == 3731) pi0_sigm = 0.004398; 
    else if(run_number == 3732) pi0_sigm = 0.004114; 
    else if(run_number == 3733) pi0_sigm = 0.003763; 
    else if(run_number == 3737) pi0_sigm = 0.004604; 
    else if(run_number == 3755) pi0_sigm = 0.004519; 
    else if(run_number == 3756) pi0_sigm = 0.004575; 
    else if(run_number == 3757) pi0_sigm = 0.004489; 
    else if(run_number == 3758) pi0_sigm = 0.004530; 
    else if(run_number == 3759) pi0_sigm = 0.004721; 
    else if(run_number == 3760) pi0_sigm = 0.004238; 
    else if(run_number == 3762) pi0_sigm = 0.004052; 
    else if(run_number == 3767) pi0_sigm = 0.004417;
    //end temp
    if(pi0_sigm == 0){
        cout<<"ERROR: pi0 width not set for run "<<run_number<<"!!!"<<endl;
        return;
    }

    // Timing offsets from Mark
    // ifstream fTimingOffset("/group/nps/mathison/analysis/NPS_Offsets/wfOffsets_3020.txt"); // x36_2
    // ifstream fTimingOffset("/group/nps/mathison/analysis/NPS_Offsets/wfOffsets_4253.txt"); // x60_4b
    // ifstream fTimingOffset("/group/nps/mathison/analysis/NPS_Offsets/wfOffsets_60_4b.txt"); // x60_4b fixed
    ifstream fTimingOffset("/group/nps/mathison/analysis/NPS_Offsets/wfOffsets_3728.txt"); // x36_5_3
    if(!fTimingOffset){
        cout<<"ERROR: Can't find the timing offsets file!!!"<<endl;
        return;
    }
    Double_t timingOffset[1080];
    for(int iblk = 0; iblk < 1080; iblk++) fTimingOffset >> timingOffset[iblk];

    // Pi0 calibration coefficients
    ifstream fpi0coef;
    // fpi0coef.open("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/Result/x60_4b_LH2_wf_v2_cycle0_cycle1/coef_pi0Calib_4253.txt"); // x60_4b wavefor fit
    fpi0coef.open(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/Result/x36_5_3_LH2_wf_cycle0_cycle1/coef_pi0Calib_%d.txt", run_number)); // x36_5_3 waveform fit
    if(!fpi0coef){
        cout<<"ERROR: Can't find the calibration coefficient file!!!"<<endl;
        return;
    }
    Double_t coefPi0[1080]; //pi0 calibration coefficients (GeV/mV) in NPS numbering Scheme
    for(int iblk = 0; iblk < 1080; iblk++) fpi0coef >> coefPi0[iblk];
    
    // Input rootfiles__________________________________________________________
    // TString dataDir = "/cache/hallc/c-nps/analysis/pass2/WF";
    TString dataDir = ".";
    TString filename = Form("nps_production_%d_%d_wf.root", run_number, iseg);

    TFile *infile = TFile::Open(Form("%s/%s", dataDir.Data(), filename.Data()));
    if(!infile || infile->IsZombie()){
        std::cerr<<Form("Error: can not open %s/%s", dataDir.Data(), filename.Data())<<endl;
        return;
    }

    // Output rootfiles__________________________________________________________
    // TString outputDir = "/volatile/hallc/nps/hhuang/farmFile/Production/DVCS";
    TString outputDir = ".";
    TString outfilename = Form("nps_production_%d_%d_wf_calib.root", run_number, iseg);

    TFile *outfile = new TFile(Form("%s/%s", outputDir.Data(), outfilename.Data()), "recreate");

    if (!outfile){ // Check if the output file was created successfully
        std::cerr << "Error creating output file.\n";
        return;
    }

    // T tree__________________________________________________________
    TTree *t_T = (TTree*)infile->Get("T");
    if (!t_T){
        std::cerr << "Error: no T tree found in the file.\n";
        // infile->Close();
        return;
    }
    
    Double_t runNb_T;
    Double_t evtNb_T;
    Double_t H_gtr_dp;
    Double_t H_gtr_ph;
    Double_t H_gtr_th;
    Double_t H_gtr_px;
    Double_t H_gtr_py;
    Double_t H_gtr_pz;
    Double_t H_react_x;
    Double_t H_react_y;
    Double_t H_react_z;
    Double_t H_hod_beta;
    Double_t H_cal_etottracknorm;
    Double_t H_cer_npeSum;

    // Int_t Ndata_NPS_cal_fly_adcCounter;
    // Double_t NPS_cal_fly_adcCounter[1080];
    // Double_t NPS_cal_fly_adcSampPulseAmp[1080];
    // Double_t NPS_cal_fly_adcSampPulseInt[1080];
    // Double_t NPS_cal_fly_adcSampPulseTime[1080];
    // Double_t NPS_cal_fly_adcSampPed[1080];

    // Disabla all and turn on the branches we need
    t_T->SetBranchStatus("*", false);
    t_T->SetBranchStatus("g.runnum", true);
    t_T->SetBranchStatus("g.evnum", true);
    t_T->SetBranchStatus("H.gtr.dp", true);
    t_T->SetBranchStatus("H.gtr.ph", true);
    t_T->SetBranchStatus("H.gtr.th", true);
    t_T->SetBranchStatus("H.gtr.px", true);
    t_T->SetBranchStatus("H.gtr.py", true);
    t_T->SetBranchStatus("H.gtr.pz", true);
    t_T->SetBranchStatus("H.react.x", true);
    t_T->SetBranchStatus("H.react.y", true);
    t_T->SetBranchStatus("H.react.z", true);
    t_T->SetBranchStatus("H.hod.beta", true);
    t_T->SetBranchStatus("H.cal.etottracknorm", true);
    t_T->SetBranchStatus("H.cer.npeSum", true);

    // NPS fADC information
    // t_T->SetBranchStatus("Ndata.NPS.cal.fly.adcCounter", true);
    // t_T->SetBranchStatus("NPS.cal.fly.adcCounter", true);
    // t_T->SetBranchStatus("NPS.cal.fly.adcSampPulseAmp", true);
    // t_T->SetBranchStatus("NPS.cal.fly.adcSampPulseInt", true);
    // t_T->SetBranchStatus("NPS.cal.fly.adcSampPulseTime", true);
    // t_T->SetBranchStatus("NPS.cal.fly.adcSampPed", true);

    //Event Level variables
    t_T->SetBranchAddress("g.runnum", &runNb_T);
    t_T->SetBranchAddress("g.evnum", &evtNb_T); // Global event number of T tree
    // HMS information
    t_T->SetBranchAddress("H.gtr.dp", &H_gtr_dp);
    t_T->SetBranchAddress("H.gtr.ph", &H_gtr_ph);
    t_T->SetBranchAddress("H.gtr.th", &H_gtr_th);
    t_T->SetBranchAddress("H.gtr.px", &H_gtr_px);
    t_T->SetBranchAddress("H.gtr.py", &H_gtr_py);
    t_T->SetBranchAddress("H.gtr.pz", &H_gtr_pz);
    t_T->SetBranchAddress("H.react.x", &H_react_x);
    t_T->SetBranchAddress("H.react.y", &H_react_y);
    t_T->SetBranchAddress("H.react.z", &H_react_z);
    t_T->SetBranchAddress("H.hod.beta", &H_hod_beta);
    t_T->SetBranchAddress("H.cal.etottracknorm", &H_cal_etottracknorm);
    t_T->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);

    // NPS fADC information
    // t_T->SetBranchAddress("Ndata.NPS.cal.fly.adcCounter", &Ndata_NPS_cal_fly_adcCounter);
    // t_T->SetBranchAddress("NPS.cal.fly.adcCounter", &NPS_cal_fly_adcCounter);
    // t_T->SetBranchAddress("NPS.cal.fly.adcSampPulseAmp", &NPS_cal_fly_adcSampPulseAmp);
    // t_T->SetBranchAddress("NPS.cal.fly.adcSampPulseInt", &NPS_cal_fly_adcSampPulseInt);
    // t_T->SetBranchAddress("NPS.cal.fly.adcSampPulseTime", &NPS_cal_fly_adcSampPulseTime);
    // t_T->SetBranchAddress("NPS.cal.fly.adcSampPed", &NPS_cal_fly_adcSampPed);

    // TChain of waveform tree___________________________________________________________
    TTree *t_wf = (TTree*)infile->Get("WF");
    if (!t_wf) {
        std::cerr << "Error: no WF tree found in the file.\n";
        return;
    }

    Double_t evtNb_wf;
    vector<Double_t> *chi2 = nullptr;
    vector<Int_t> *wfnpulse = nullptr;
    vector<Double_t> *wfampl = nullptr;
    vector<Double_t> *wftime = nullptr;
    t_wf->SetBranchAddress("evt", &evtNb_wf); // Global event number of wf tree
    t_wf->SetBranchAddress("chi2", &chi2); // Number of pulses found by TSpectrum in this block
    t_wf->SetBranchAddress("wfnpulse", &wfnpulse); // Number of pulses found by TSpectrum in this block 
    t_wf->SetBranchAddress("wfampl", &wfampl); // Waveform amplitude of each fitted pulse (mV)
    t_wf->SetBranchAddress("wftime", &wftime); // Pulse time of each fitted pulse per block (ns)

    // Use the TVirtualIndex to access the entry corresponding to the on in T tree
    TVirtualIndex *vIdx = t_wf->GetTreeIndex();
    TTreeIndex *idx = dynamic_cast<TTreeIndex *>(vIdx);
    if (!idx)
    {
        std::cerr << "Error: no TTreeIndex found in the output tree.\n";
        infile->Close();
        return;
    }

    // Get the index array
    Long64_t *indexArray = idx->GetIndex();
    Int_t nevt_wf = t_wf->GetEntries();
    Int_t nevt_T = t_T->GetEntries();
    Double_t lastevent = 0.;

    if (nevt_wf != nevt_T){ // Check the number of events in both trees
        std::cerr << "Error: number of events in WF tree (" << nevt_wf << ") does not match number of events in T tree (" << nevt_T << ").\n";
        return;
    }

    // Connect to the database and get kinematics variables________________________________________________
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    Double_t Beam_energy = *db->GetEntry_d("BEAM_param_Energy", run_number); // Beam energy in GeV
    Double_t HMS_mom = *db->GetEntry_d("SIMU_param_HMSmomentum", run_number); // HMS central momentum in GeV/c
    Double_t HMS_angle = *db->GetEntry_d("SIMU_param_HMSangle", run_number); // HMS angle in rad
    Double_t Target_amu = *db->GetEntry_d("TARGET_param_Amu", run_number); // Target amu

    // The angle and distance of calorimeter
    Double_t NPS_dist = *db->GetEntry_d("CALO_geom_Dist", run_number); // NPS distence in cm
    Double_t NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", run_number); // NPS angle in rad
    Int_t *caloMaskBlock = new Int_t[1080]; // Get the mask block information in NPS numbering Scheme
    caloMaskBlock = db->GetEntry_i("CALO_flag_MaskBlock", run_number);

    // Elastic coefficients
    // Double_t *coefElas = new Double_t[1080]; // Get the elastic coefficients (GeV/mV) in NPS numbering Scheme
    // coefElas = db->GetEntry_d("CALO_calib_ElasCoef", run_number);

    cout<<"Start clustering for Run "<<run_number<<", segment "<<iseg<<" with cluster threshold "<<clusTrsH<<" GeV"<<endl;
    cout<<"========== Run information =========="<<endl;
    cout<<"Beam energy: "<<Beam_energy<<" GeV"<<endl;
    cout<<"HMS momentum: "<<HMS_mom<<" GeV/c"<<endl;
    cout<<"HMS angle: "<<HMS_angle*TMath::RadToDeg()<<" degrees"<<endl;
    cout<<"Target amu: "<<Target_amu<<" GeV/c^2"<<endl;
    cout<<"NPS distance: "<<NPS_dist<<" cm"<<endl;
    cout<<"NPS angle: "<<NPS_angle*TMath::RadToDeg()<<" degrees"<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== NPS maske blocks =========="<<endl;
    for(int i = 0; i < 1080; i++) if(caloMaskBlock[i] == 1) cout<<i<<" ";
    cout<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== NPS pi0 calibration coefficients =========="<<endl;
    for(int i = 0; i < 1080; i++) cout<<coefPi0[i]<<" ";
    cout<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== NPS timing offsets =========="<<endl;
    for(int iblk = 0; iblk < 1080; iblk++) cout<<timingOffset[iblk]<<" ";
    cout<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;

    // Initial settings for clustering__________________________________________________
    TLorentzVector beam(0, 0, Beam_energy, Beam_energy);
    TLorentzVector p0(0, 0, 0, m_p);

    TDVCSEvent *ev = new TDVCSEvent(run_number);
    TCaloEvent *caloev = new TCaloEvent(run_number);
    ev->GetGeometry()->SetCaloTheta(-1*NPS_angle); // rad Note: the angle have to be negative here
    ev->GetGeometry()->SetCaloDist(NPS_dist); // cm

    // Ontput Tree for production data___________________________________________________________
    TTree *t_prod = new TTree("t_prod", "t_prod");

    // Variables for NPS information in the tree
    Int_t nclust;
    vector<double> clusE;
    vector<int> clusSize;
    vector<double> clusX;
    vector<double> clusX_corr;
    vector<double> clusY;
    vector<double> clusY_corr;
    vector<double> clusT;
    vector<double> clusDepth;
    vector<double> trkPx;
    vector<double> trkPy;
    vector<double> trkPz;
    vector<double> trkE;
    Double_t M;
    Double_t Mx2;

    //Event Level variables
    t_prod->Branch("g.runnum", &runNb_T);
    t_prod->Branch("g.evnum", &evtNb_T); // Global event number of T tree
    // HMS information
    t_prod->Branch("H.gtr.dp", &H_gtr_dp);
    t_prod->Branch("H.gtr.ph", &H_gtr_ph);
    t_prod->Branch("H.gtr.th", &H_gtr_th);
    t_prod->Branch("H.gtr.px", &H_gtr_px);
    t_prod->Branch("H.gtr.py", &H_gtr_py);
    t_prod->Branch("H.gtr.pz", &H_gtr_pz);
    t_prod->Branch("H.react.x", &H_react_x);
    t_prod->Branch("H.react.y", &H_react_y);
    t_prod->Branch("H.react.z", &H_react_z);
    t_prod->Branch("H.hod.beta", &H_hod_beta);
    t_prod->Branch("H.cal.etottracknorm", &H_cal_etottracknorm);
    t_prod->Branch("H.cer.npeSum", &H_cer_npeSum);
    // // NPS waveform fitting results
    // t_prod->Branch("evt", &evtNb_wf); // Global event number of wf tree
    // t_prod->Branch("chi2", &chi2); // Number of pulses found by TSpectrum in this block
    // t_prod->Branch("wfnpulse", &wfnpulse); // Number of pulses found by TSpectrum in this block 
    // t_prod->Branch("wfampl", &wfampl); // Waveform amplitude of each fitted pulse (mV)
    // t_prod->Branch("wftime", &wftime); // Pulse time of each fitted pulse per block (ns)
    // NPS information
    t_prod->Branch("NPS.prod.nclust", &nclust, "nclust/I");
    t_prod->Branch("NPS.prod.clusE", &clusE);
    t_prod->Branch("NPS.prod.clusSize", &clusSize);
    t_prod->Branch("NPS.prod.clusX", &clusX);
    t_prod->Branch("NPS.prod.clusX.corr", &clusX_corr);
    t_prod->Branch("NPS.prod.clusY", &clusY);
    t_prod->Branch("NPS.prod.clusY.corr", &clusY_corr);
    t_prod->Branch("NPS.prod.clusZ", &NPS_dist);
    t_prod->Branch("NPS.prod.clusT", &clusT);
    t_prod->Branch("NPS.prod.clusDepth", &clusDepth);
    t_prod->Branch("NPS.prod.trk.px", &trkPx);
    t_prod->Branch("NPS.prod.trk.py", &trkPy);
    t_prod->Branch("NPS.prod.trk.pz", &trkPz);
    t_prod->Branch("NPS.prod.trk.ene", &trkE);
    t_prod->Branch("NPS.prod.M", &M, "M/D");
    t_prod->Branch("NPS.prod.Mx2", &Mx2, "Mx2/D");

    // Ontput Tree for accidental data #1___________________________________________________________
    TTree *t_accdt1 = new TTree("t_accdt1", "t_accdt1");
    //Event Level variables
    t_accdt1->Branch("g.runnum", &runNb_T);
    t_accdt1->Branch("g.evnum", &evtNb_T); // Global event number of T tree
    // HMS information
    t_accdt1->Branch("H.gtr.dp", &H_gtr_dp);
    t_accdt1->Branch("H.gtr.ph", &H_gtr_ph);
    t_accdt1->Branch("H.gtr.th", &H_gtr_th);
    t_accdt1->Branch("H.gtr.px", &H_gtr_px);
    t_accdt1->Branch("H.gtr.py", &H_gtr_py);
    t_accdt1->Branch("H.gtr.pz", &H_gtr_pz);
    t_accdt1->Branch("H.react.x", &H_react_x);
    t_accdt1->Branch("H.react.y", &H_react_y);
    t_accdt1->Branch("H.react.z", &H_react_z);
    t_accdt1->Branch("H.hod.beta", &H_hod_beta);
    t_accdt1->Branch("H.cal.etottracknorm", &H_cal_etottracknorm);
    t_accdt1->Branch("H.cer.npeSum", &H_cer_npeSum);
    // NPS information
    t_accdt1->Branch("NPS.prod.nclust", &nclust, "nclust/I");
    t_accdt1->Branch("NPS.prod.clusE", &clusE);
    t_accdt1->Branch("NPS.prod.clusSize", &clusSize);
    t_accdt1->Branch("NPS.prod.clusX", &clusX);
    t_accdt1->Branch("NPS.prod.clusX.corr", &clusX_corr);
    t_accdt1->Branch("NPS.prod.clusY", &clusY);
    t_accdt1->Branch("NPS.prod.clusY.corr", &clusY_corr);
    t_accdt1->Branch("NPS.prod.clusZ", &NPS_dist);
    t_accdt1->Branch("NPS.prod.clusT", &clusT);
    t_accdt1->Branch("NPS.prod.clusDepth", &clusDepth);
    t_accdt1->Branch("NPS.prod.trk.px", &trkPx);
    t_accdt1->Branch("NPS.prod.trk.py", &trkPy);
    t_accdt1->Branch("NPS.prod.trk.pz", &trkPz);
    t_accdt1->Branch("NPS.prod.trk.ene", &trkE);
    t_accdt1->Branch("NPS.prod.M", &M);
    t_accdt1->Branch("NPS.prod.Mx2", &Mx2);

    // Ontput Tree for accidental data #2___________________________________________________________
    TTree *t_accdt2 = new TTree("t_accdt2", "t_accdt2");
    //Event Level variables
    t_accdt2->Branch("g.runnum", &runNb_T);
    t_accdt2->Branch("g.evnum", &evtNb_T); // Global event number of T tree
    // HMS information
    t_accdt2->Branch("H.gtr.dp", &H_gtr_dp);
    t_accdt2->Branch("H.gtr.ph", &H_gtr_ph);
    t_accdt2->Branch("H.gtr.th", &H_gtr_th);
    t_accdt2->Branch("H.gtr.px", &H_gtr_px);
    t_accdt2->Branch("H.gtr.py", &H_gtr_py);
    t_accdt2->Branch("H.gtr.pz", &H_gtr_pz);
    t_accdt2->Branch("H.react.x", &H_react_x);
    t_accdt2->Branch("H.react.y", &H_react_y);
    t_accdt2->Branch("H.react.z", &H_react_z);
    t_accdt2->Branch("H.hod.beta", &H_hod_beta);
    t_accdt2->Branch("H.cal.etottracknorm", &H_cal_etottracknorm);
    t_accdt2->Branch("H.cer.npeSum", &H_cer_npeSum);
    // NPS information
    t_accdt2->Branch("NPS.prod.nclust", &nclust, "nclust/I");
    t_accdt2->Branch("NPS.prod.clusE", &clusE);
    t_accdt2->Branch("NPS.prod.clusSize", &clusSize);
    t_accdt2->Branch("NPS.prod.clusX", &clusX);
    t_accdt2->Branch("NPS.prod.clusX.corr", &clusX_corr);
    t_accdt2->Branch("NPS.prod.clusY", &clusY);
    t_accdt2->Branch("NPS.prod.clusY.corr", &clusY_corr);
    t_accdt2->Branch("NPS.prod.clusZ", &NPS_dist);
    t_accdt2->Branch("NPS.prod.clusT", &clusT);
    t_accdt2->Branch("NPS.prod.clusDepth", &clusDepth);
    t_accdt2->Branch("NPS.prod.trk.px", &trkPx);
    t_accdt2->Branch("NPS.prod.trk.py", &trkPy);
    t_accdt2->Branch("NPS.prod.trk.pz", &trkPz);
    t_accdt2->Branch("NPS.prod.trk.ene", &trkE);
    t_accdt2->Branch("NPS.prod.M", &M);
    t_accdt2->Branch("NPS.prod.Mx2", &Mx2);

    // Ontput Tree for pi0 contamination___________________________________________________________
    TTree *t_pi0sub = new TTree("t_pi0sub", "t_pi0sub");
    //Event Level variables
    t_pi0sub->Branch("g.runnum", &runNb_T);
    t_pi0sub->Branch("g.evnum", &evtNb_T); // Global event number of T tree
    // HMS information
    t_pi0sub->Branch("H.gtr.dp", &H_gtr_dp);
    t_pi0sub->Branch("H.gtr.ph", &H_gtr_ph);
    t_pi0sub->Branch("H.gtr.th", &H_gtr_th);
    t_pi0sub->Branch("H.gtr.px", &H_gtr_px);
    t_pi0sub->Branch("H.gtr.py", &H_gtr_py);
    t_pi0sub->Branch("H.gtr.pz", &H_gtr_pz);
    t_pi0sub->Branch("H.react.x", &H_react_x);
    t_pi0sub->Branch("H.react.y", &H_react_y);
    t_pi0sub->Branch("H.react.z", &H_react_z);
    t_pi0sub->Branch("H.hod.beta", &H_hod_beta);
    t_pi0sub->Branch("H.cal.etottracknorm", &H_cal_etottracknorm);
    t_pi0sub->Branch("H.cer.npeSum", &H_cer_npeSum);
    // photons from pi0 contamination
    Int_t N0, N1, N2;
    Double_t weight;
    vector<double> phXc;
    vector<double> phYc;
    vector<double> phPx;
    vector<double> phPy;
    vector<double> phPz;
    vector<double> phE;
    vector<double> phMx2;
    t_pi0sub->Branch("pi0.cont.Xc", &phXc);
    t_pi0sub->Branch("pi0.cont.Yc", &phYc);
    t_pi0sub->Branch("pi0.cont.trk.px", &phPx);
    t_pi0sub->Branch("pi0.cont.trk.py", &phPy);
    t_pi0sub->Branch("pi0.cont.trk.pz", &phPz);
    t_pi0sub->Branch("pi0.cont.trk.ene", &phE);
    t_pi0sub->Branch("pi0.cont.Mx2", &phMx2);
    t_pi0sub->Branch("pi0.cont.N0", &N0);
    t_pi0sub->Branch("pi0.cont.N1", &N1);
    t_pi0sub->Branch("pi0.cont.N2", &N2);
    t_pi0sub->Branch("pi0.cont.weight", &weight);

    // Output histograms__________________________________________________________
    TH1F *h_reactX = new TH1F("h_reactX", "H.react.x;Reaction X [cm];Counts", 400, -0.4, 0.1);
    TH1F *h_reactY = new TH1F("h_reactY", "H.react.y;Reaction Y [cm];Counts", 400, -0.2, 0.2);
    TH1F *h_reactZ = new TH1F("h_reactZ", "H.react.z;Reaction Z [cm];Counts", 400, -20, 20);

    TH1F *h_th = new TH1F("h_th", "H.gtr.th", 400, -2, 2);
    TH2F *hh_dp_ph = new TH2F("hh_dp_ph", "H.gtr.dp vs. H.gtr.ph;", 100, -0.05, 0.05, 300, -15, 15);

    TH1F *h_etottracknorm = new TH1F("h_etottracknorm", "H.cal.etottracknorm", 1000, 0, 10);
    TH1F *h_npeSum = new TH1F("h_npeSum", "H.cer.npeSum", 1000, 0, 10);
    TH1F *h_beta = new TH1F("h_beta", "H.hod.beta;#beta;Counts", 1500, 0, 1.5);
    TH1F *h_NpsTime = new TH1F("h_NpsTime", "NPS timing with offset correction (ampl. > 10 mV);Pulse time [ns];Number of events", 20000, -100, 100); // When using the waveform tree

    TH1F *h_nclust = new TH1F("h_nclust", "Number of clusters;Number of clusters;Counts", 21, -0.5, 20.5);

    // Event loop
    Int_t nevt = nevt_T;
    // Int_t nevt = 10000; // for testing
    Double_t tmean = 0.;
    for(int it = 0; it < 3; it++){ // Loop over three different time windows
        if(it == 0) tmean = 0.; // Production time window = 0+-3ns
        else if(it == 1) tmean = -8.; // Accidental#1 time window = -8+-3ns
        else if(it == 2) tmean = 8.; // Accidental#1 time window = 8+-3ns

        for(Int_t ievt = 0; ievt < nevt; ievt++){ // Event loop
            t_T->GetEntry(ievt);
            t_wf->GetEntry(indexArray[ievt]); // indexArray[i] is the original-entry number for the i-th smallest evt
            if(ievt%100000==0) cout << "Looking at entry = " << ievt << "  (" << 100.*ievt/nevt_T << "%)" << endl;

            // Initialize for production and accidental trees
            nclust = 0;
            clusE.clear();
            clusSize.clear();
            clusX.clear();
            clusX_corr.clear();
            clusY.clear();
            clusY_corr.clear();
            clusT.clear();
            clusDepth.clear();
            trkPx.clear();
            trkPy.clear();
            trkPz.clear();
            trkE.clear();
            M = -999;
            Mx2 = -999;
            // Initialize for pi0 contamination tree
            N0 = 0;
            N1 = 0;
            N2 = 0;
            weight = 1.;
            phXc.clear();
            phYc.clear();
            phPx.clear();
            phPy.clear();
            phPz.clear();
            phE.clear();
            phMx2.clear();

            if(abs(H_react_x) > 100 || abs(H_react_y) > 100 || abs(H_react_z) > 100){ 
                // Skip photon reconstruction for strange vertex postion of 1e+38
                // Fill the Tree then jump to next event without clustering
                if(it == 0){t_prod->Fill(); continue;}
                else if(it == 1){t_accdt1->Fill(); continue;}
                else if(it == 2){t_accdt2->Fill(); continue;}
            }
            // Fill the histograms before selection
            h_reactX->Fill(H_react_x);
            h_reactY->Fill(H_react_y);
            h_reactZ->Fill(H_react_z);
            h_th->Fill(H_gtr_th);
            hh_dp_ph->Fill(H_gtr_ph, H_gtr_dp);
            h_etottracknorm->Fill(H_cal_etottracknorm);
            h_npeSum->Fill(H_cer_npeSum);
            h_beta->Fill(H_hod_beta);
            

            // Using the waveform tree_____________________________________________
            Int_t flatIdx = 0; // accumulated number of pulses in previous blocks
            for (size_t iblk=0; iblk<wfnpulse->size(); iblk++){
                Int_t nPulse = (*wfnpulse)[iblk];
                Int_t iblk_sim = bnConv_NewToOld(iblk); // block number with simulation numbering scheme
                TCaloBlock *block = caloev->AddBlock(iblk_sim); // add block with simulation numbering scheme
                if (nPulse > 0 && (*chi2)[iblk] > 0 && !caloMaskBlock[iblk]){ // if we want to read the pulse in this block
                    Int_t pulseIdx = flatIdx + 0; // index of the first pulse in this block
                    for (Int_t p = 0; p < nPulse; p++){ // Find the pulse closest to the mean time
                        if(abs((*wftime)[flatIdx + p]+timingOffset[iblk]-tmean) <= abs((*wftime)[pulseIdx]+timingOffset[iblk]-tmean)) pulseIdx = flatIdx + p; 
                    }
                    if(it == 0 && (*wfampl)[pulseIdx] > 10) h_NpsTime->Fill((*wftime)[pulseIdx]+timingOffset[iblk]); // Fill the histogram with timing offset correction

                    if(abs((*wftime)[pulseIdx]+timingOffset[iblk]-tmean) < 3) block->AddPulse((*wfampl)[pulseIdx] * coefPi0[iblk], (*wftime)[pulseIdx]+timingOffset[iblk]-tmean); // Add energy and timing to the block
                } // block selection
                flatIdx += nPulse;
            } // loop over 1080 blocks

            caloev->TriggerSim(clusTrsH);     // Energy threshold of every 2x2 blocks
            caloev->DoClustering(-3, 3); // Time window cut (-3,3) [ns] for clustering
            ev->SetCaloEvent(caloev);
            ev->SetVertex(H_react_x, H_react_y, H_react_z);

            nclust = caloev->GetNbClusters();
            for(int iclus = 0; iclus < nclust; iclus++){
                caloev->GetCluster(iclus)->Analyze(1);

                clusE.push_back(caloev->GetCluster(iclus)->GetEnergy());
                clusSize.push_back(caloev->GetCluster(iclus)->GetClusSize());
                clusX.push_back(caloev->GetCluster(iclus)->GetX());
                clusY.push_back(caloev->GetCluster(iclus)->GetY());

                Float_t a_out, x_corr_out, y_corr_out;
                TLorentzVector photon = ev->GetPhoton(iclus, 7, 0, a_out, x_corr_out, y_corr_out);

                clusDepth.push_back(a_out);
                clusX_corr.push_back(x_corr_out);
                clusY_corr.push_back(y_corr_out);

                trkPx.push_back(photon.Px());
                trkPy.push_back(photon.Py());
                trkPz.push_back(photon.Pz());
                trkE.push_back(photon.E());

                // calculate cluster time
                Double_t sum_et = 0.;
                for(int iblk = 0; iblk < caloev->GetCluster(iclus)->GetClusSize(); iblk++){
                    Double_t blockE = caloev->GetCluster(iclus)->GetBlock(iblk)->GetEnergy(0);
                    Double_t blockT = caloev->GetCluster(iclus)->GetBlock(iblk)->GetTime(0);
                    if(blockE > 0) sum_et+= blockE*blockT;
                    else(cout<<"Warning: block energy <= 0 ("<<blockE<<" GeV)"<<endl);
                }
                clusT.push_back(sum_et/caloev->GetCluster(iclus)->GetEnergy());
            }

            h_nclust->Fill(nclust);

            Double_t H_gtr_e = sqrt(H_gtr_px*H_gtr_px+H_gtr_py*H_gtr_py+H_gtr_pz*H_gtr_pz+0.000511*0.000511);
            TLorentzVector kpvec(H_gtr_px, H_gtr_py, H_gtr_pz, H_gtr_e);
            if(caloev->GetNbClusters() == 1){
                TLorentzVector photon1 = ev->GetPhoton(0);
                Mx2 = (beam + p0 - kpvec - photon1).Mag2();
            } // clus==1

            else if(caloev->GetNbClusters() == 2){
                TLorentzVector photon1 = ev->GetPhoton(0);
                TLorentzVector photon2 = ev->GetPhoton(1);
                M = (photon1 + photon2).M();
                Mx2 = (beam + p0 - kpvec - photon1 - photon2).Mag2();

                // pi0 subtraction__________________________________________________________________
                if(it == 0){
                    // acceptance cut for photons
                    Int_t nCol_inMgnt = 0; // number of columns under the shadow of magnet
                    if(290 < NPS_dist && NPS_dist < 310) nCol_inMgnt = 5; // col0-4 when NPS@300cm
                    else if(340 < NPS_dist && NPS_dist < 360) nCol_inMgnt = 4; // col0-3 when NPS@350cm
                    else if(390 < NPS_dist && NPS_dist < 410) nCol_inMgnt = 2; // col0-1 when NPS@400cm

                    Double_t x1 = caloev->GetCluster(0)->GetX();
                    Double_t y1 = caloev->GetCluster(0)->GetY();
                    Double_t ene1 = caloev->GetCluster(0)->GetEnergy();
                    
                    Double_t x2 = caloev->GetCluster(1)->GetX();
                    Double_t y2 = caloev->GetCluster(1)->GetY();
                    Double_t ene2 = caloev->GetCluster(1)->GetEnergy();

                    Double_t Xcut_min = -2.16*(15-1-nCol_inMgnt);
                    Double_t Xcut_max = 2.16*(15-1);
                    Double_t Ycut_min = -2.16*(18-1);
                    Double_t Ycut_max = 2.16*(18-1);

                    if(abs(M-0.1349766) < pi0_sigm // pi0 mass window cut, pi0_sigm is run dependent
                    && ene1 > 0.5 && ene2 > 0.5
                    && Xcut_min < x1 && x1 < Xcut_max && Ycut_min < y1 && y1 < Ycut_max
                    && Xcut_min < x2 && x2 < Xcut_max && Ycut_min < y2 && y2 < Ycut_max){
                        
                        TLorentzVector pion_tlv = photon1 + photon2;
                        for(int idc = 0; idc < ndecay; idc++){
                            Double_t cosalpha = 2.*gRandom->Rndm()-1.; // uniform polor angle cos(alpha) between -1 and 1
                            Double_t phi2 = 2.*TMath::Pi()*gRandom->Rndm(); // uniform azimuthal angle between 0 and 2*pi
                            // Decay in pion rest frame
                            Double_t e_photon = 0.5*M;
                            Double_t px_photon = e_photon*TMath::Sin(TMath::ACos(cosalpha))*TMath::Cos(phi2);
                            Double_t py_photon = e_photon*TMath::Sin(TMath::ACos(cosalpha))*TMath::Sin(phi2);
                            Double_t pz_photon = e_photon*cosalpha;
                            TLorentzVector photonA(px_photon, py_photon, pz_photon, e_photon);
                            TLorentzVector photonB(-1*px_photon, -1*py_photon, -1*pz_photon, e_photon);
                            // Boost to lab frame
                            photonA.Boost(pion_tlv.BoostVector());
                            photonB.Boost(pion_tlv.BoostVector());
                            Double_t eneA = photonA.E();
                            Double_t eneB = photonB.E();
                            // Impact position on NPS surface
                            Double_t beta = TMath::ATan(photonA.Px() / photonA.Pz());
                            Double_t xcA = NPS_dist*TMath::Tan(beta-NPS_angle)-H_react_z*TMath::Sin(beta)/TMath::Cos(beta-NPS_angle);
                            Double_t ycA = (TMath::Sqrt(NPS_dist*NPS_dist+xcA*xcA)-H_react_z*(TMath::Cos(beta)+TMath::Sin(beta)*TMath::Tan(beta-NPS_angle))) * photonA.Py();
                            ycA = ycA/TMath::Sqrt(photonA.Px()*photonA.Px()+photonA.Pz()*photonA.Pz());

                            beta = TMath::ATan(photonB.Px() / photonB.Pz());
                            Double_t xcB = NPS_dist*TMath::Tan(beta-NPS_angle)-H_react_z*TMath::Sin(beta)/TMath::Cos(beta-NPS_angle);
                            Double_t ycB = (TMath::Sqrt(NPS_dist*NPS_dist+xcB*xcB)-H_react_z*(TMath::Cos(beta)+TMath::Sin(beta)*TMath::Tan(beta-NPS_angle)))*photonB.Py();
                            ycB = ycB/TMath::Sqrt(photonB.Px()*photonB.Px()+photonB.Pz()*photonB.Pz());

                            // if the decayed photon pass the acceptance cut
                            Bool_t passA = false;
                            Bool_t passB = false;

                            if(eneA > 0.5 
                            && Xcut_min < xcA && xcA < Xcut_max 
                            && Ycut_min < ycA && ycA < Ycut_max) passA = true;

                            if(eneB > 0.5 
                            && Xcut_min < xcB && xcB < Xcut_max 
                            && Ycut_min < ycB && ycB < Ycut_max) passB = true;

                            if(!passA && !passB) N0++; // no photon pass the acceptance cut
                            if(passA && passB) N2++; // both photons pass the acceptance cut
                            if((!passA && passB) || (passA && !passB)){ // only one photon pass the acceptance cut
                                N1++;
                                if(passA){
                                    Double_t Mx2_A = (beam + p0 - kpvec - photonA).Mag2();
                                    phXc.push_back(xcA);
                                    phYc.push_back(ycA);
                                    phPx.push_back(photonA.Px());
                                    phPy.push_back(photonA.Py());
                                    phPz.push_back(photonA.Pz());
                                    phE.push_back(eneA);
                                    phMx2.push_back(Mx2_A);
                                }
                                if(passB){
                                    Double_t Mx2_B = (beam + p0 - kpvec - photonB).Mag2();
                                    phXc.push_back(xcB);
                                    phYc.push_back(ycB);
                                    phPx.push_back(photonB.Px());
                                    phPy.push_back(photonB.Py());
                                    phPz.push_back(photonB.Pz());
                                    phE.push_back(eneB);
                                    phMx2.push_back(Mx2_B);
                                }
                            } // if only one photon pass
                        } // loop over ndecay
                        weight = 1./N2;
                        t_pi0sub->Fill();
                    }
                }
            } // clus==2

            if(it == 0) t_prod->Fill();
            else if(it == 1) t_accdt1->Fill();
            else if(it == 2) t_accdt2->Fill();

            // reinitialization for next event
            caloev->Reset();
        } // Event loop: production and accidental tree
    }// Loop over different mean time(-8,0,8)ns

    cout<<endl;
    cout<<"Number of event summary ===================="<<endl;
    cout<<"T tree: "<<nevt_T<<endl;
    cout<<"waveform tree: "<<nevt_wf<<endl;
    cout<<"production tree: "<<t_prod->GetEntries()<<endl;
    cout<<"accidental#1 tree: "<<t_accdt1->GetEntries()<<endl;
    cout<<"accidental#2 tree: "<<t_accdt2->GetEntries()<<endl;
    cout<<"pi0 contamination tree: "<<t_pi0sub->GetEntries()<<endl;
    cout<<"============================================"<<endl;
    cout<<endl;
    if(nevt_T != t_prod->GetEntries()){
        std::cerr<<"Error: number of events in T tree and production tree do not match!"<<endl;
    }

    // Save the output file__________________________________________________________
    outfile->cd();
    t_prod->Write();
    t_accdt1->Write();
    t_accdt2->Write();
    t_pi0sub->Write();
    
    h_reactX->Write();
    h_reactY->Write();
    h_reactZ->Write();
    h_beta->Write();
    h_etottracknorm->Write();
    h_npeSum->Write();
    h_th->Write();
    hh_dp_ph->Write();
    h_nclust->Write();
    
    outfile->Close();
}
