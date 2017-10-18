//  used initially to create angular correlation histograms and graphs from individually simmulated "data" and the Z distributions


TH1D * hT, *Z0, *Z2, *Z4;
TGraphAsymmErrors *acT, *Z0ac, *Z2ac, *Z4ac;


void collect(double E1a, double E2a){

TAngularCorrelation* ac = new TAngularCorrelation();
ac->GenerateMaps(16,110);

double E1b = E1a+1;
double E2b = E2a+1; 

// template
for (int i=1;i<=10;++i){

	TFile* file_data = new TFile(Form("template/Converted%i.root",i));
	TDirectory* NDdir_data = file_data->GetDirectory("GriffinND");
	THnSparse* my_thn = (THnSparse*) NDdir_data->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
	TH2D* h2 = ac->Create2DSlice(my_thn,E1a,E1b,kFALSE,kFALSE);
    TH1D* h1 = ac->IntegralSlices(h2,E2a,E2b);
	TH1D* h = ac->DivideByWeights(h1,false,false); //angular index
	hT->Add(h);
}


// Z0
for (int i=1;i<=10;++i){

	TFile* file_data = new TFile(Form("Z0/Converted%i.root",i));
	TDirectory* NDdir_data = file_data->GetDirectory("GriffinND");
	THnSparse* my_thn = (THnSparse*) NDdir_data->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
	TH2D* h2 = ac->Create2DSlice(my_thn,E1a,E1b,kFALSE,kFALSE);
    TH1D* h1 = ac->IntegralSlices(h2,E2a,E2b);
	TH1D* h = ac->DivideByWeights(h1,false,false); //angular index
	Z0->Add(h);
}

// Z2
for (int i=1;i<=10;++i){

	TFile* file_data = new TFile(Form("Z2/Converted%i.root",i));
	TDirectory* NDdir_data = file_data->GetDirectory("GriffinND");
	THnSparse* my_thn = (THnSparse*) NDdir_data->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
	TH2D* h2 = ac->Create2DSlice(my_thn,E1a,E1b,kFALSE,kFALSE);
    TH1D* h1 = ac->IntegralSlices(h2,E2a,E2b);
	TH1D* h = ac->DivideByWeights(h1,false,false); //angular index
	Z2->Add(h);
}

// Z4
for (int i=1;i<=10;++i){

	TFile* file_data = new TFile(Form("Z4/Converted%i.root",i));
	TDirectory* NDdir_data = file_data->GetDirectory("GriffinND");
	THnSparse* my_thn = (THnSparse*) NDdir_data->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
	TH2D* h2 = ac->Create2DSlice(my_thn,E1a,E1b,kFALSE,kFALSE);
    TH1D* h1 = ac->IntegralSlices(h2,E2a,E2b);
	TH1D* h = ac->DivideByWeights(h1,false,false); //angular index
	Z4->Add(h);
}

// create angular correlations 
acT = ac->CreateGraphFromHst(temp,false,false);
Z0ac = ac->CreateGraphFromHst(Z0,false,false);
Z2ac = ac->CreateGraphFromHst(Z2,false,false);
Z4ac = ac->CreateGraphFromHst(Z4,false,false);

// write to file
TFile *file = new TFile("simulations.root","NEW");

TCanvas* hT_c = new TCanvas("hT");
hT_c->cd();
hT->Write("hT");

TCanvas* acT_c = new TCanvas("acT");
acT_c->cd();
acT->Write("acT");

TCanvas* Z0_c = new TCanvas("Z0");
Z0_c->cd();
Z0->Write("Z0");
TCanvas* Z0ac_c = new TCanvas("Z0ac");
Z0ac_c->cd();
Z0ac->Write("Z0ac");

TCanvas* Z2_c = new TCanvas("Z2");
Z2_c->cd();
Z2->Write("Z2");
TCanvas* Z2ac_c = new TCanvas("Z2ac");
Z2ac_c->cd();
Z2ac->Write("Z2ac");

TCanvas* Z4_c = new TCanvas("Z4");
Z4_c->cd();
Z4->Write("Z4");
TCanvas* Z4ac_c = new TCanvas("Z4ac");
Z4ac_c->cd();
Z4ac->Write("Z4ac");
file->Close(); 

}

