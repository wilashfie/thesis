TH1D *W_0, *W_2, *W_4;

bool normalize = false;

double a_2 = 0.1020408163265306;
double a_4 = 0.009070294784580489;

//double a2 = 0.3571428571428571;
//double a4 = 1.1428571428571415;

Double_t Lpoly(Double_t *x, Double_t *p)
{
	Double_t val = p[0]*(1 + p[1]*::ROOT::Math::legendre(2,x[0]) + p[2]*::ROOT::Math::legendre(4,x[0]));
	return val;
}

// define linear combination fit 
Double_t fitFcn(Double_t *x, Double_t *p){
	Double_t i = x[0];
	Int_t bin_0 = W_0->GetXaxis()->FindBin(i);
	Double_t w0 = (1-p[1]-p[2])*W_0->GetBinContent(bin_0);
	Int_t bin_2 = W_2->GetXaxis()->FindBin(i);
	Double_t w2 = p[1]*W_2->GetBinContent(bin_2);
	Int_t bin_4 = W_4->GetXaxis()->FindBin(i);
	Double_t w4 = p[2]*W_4->GetBinContent(bin_4);
	Double_t out = p[0]*(w0 + w2 + w4);
	return out;
}

void fitHistoRAWnu(){

gStyle->SetOptFit(0);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

// create angular correlation for data 
TDirectory* NDdir = gFile->GetDirectory("GriffinND");
THnSparse* my_thn = (THnSparse*) NDdir->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");

TAngularCorrelation* ac = new TAngularCorrelation();
ac->GenerateMaps(16,110);

TH2D* h2 = ac->Create2DSlice(my_thn,1173,1174,kFALSE,kFALSE);
TH1D* h1 = ac->IntegralSlices(h2,1332,1333);
TH1D* hsto = ac->DivideByWeights(h1,false,false); //angular index

if(normalize){
Double_t norm = hsto->GetEntries(); 
hsto->Scale(1/norm); // normalize 
cout << "norm = " << norm << endl;
}

hsto->SetTitle(
	";"
	"Angular Index;"
	"Counts;"
	);
hsto->SetMarkerColor(38);
hsto->SetMarkerStyle(8);
hsto->SetMarkerSize(2);
hsto->GetYaxis()->SetNdivisions(505);
hsto->Draw();
//hsto->GetYaxis()->SetRangeUser(3000,4000);

// load template files 
TFile* file_a0 = new TFile("ultra_a0Nu.root");
TFile* file_a2 = new TFile("ultra_a2Nu.root");
TFile* file_a4 = new TFile("ultra_a4Nu.root");

TDirectory* NDdir_a0 = file_a0->GetDirectory("GriffinND");
TDirectory* NDdir_a2 = file_a2->GetDirectory("GriffinND");
TDirectory* NDdir_a4 = file_a4->GetDirectory("GriffinND");

THnSparse* thn_a0 = (THnSparse*) NDdir_a0->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
THnSparse* thn_a2 = (THnSparse*) NDdir_a2->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
THnSparse* thn_a4 = (THnSparse*) NDdir_a4->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");

TAngularCorrelation* acNU = new TAngularCorrelation();
acNU->GenerateMaps(16,110);



// create angular index for each template 
TH2D* h2_a0 = acNU->Create2DSlice(thn_a0,1173,1174,kFALSE,kFALSE);
h2_a0->SetName("h2_a0");
TH2D* h2_a2 = acNU->Create2DSlice(thn_a2,1173,1174,kFALSE,kFALSE);
h2_a2->SetName("h2_a2");
TH2D* h2_a4 = acNU->Create2DSlice(thn_a4,1173,1174,kFALSE,kFALSE);
h2_a4->SetName("h2_a4");

TH1D* h1_a0 = acNU->IntegralSlices(h2_a0,1332,1333);
TH1D* h1_a2 = acNU->IntegralSlices(h2_a2,1332,1333);
TH1D* h1_a4 = acNU->IntegralSlices(h2_a4,1332,1333);

W_0 = acNU->DivideByWeights(h1_a0,false,false); 
W_2 = acNU->DivideByWeights(h1_a2,false,false); 
W_4 = acNU->DivideByWeights(h1_a4,false,false); 



if(normalize){
Double_t norm0 = W_0->GetEntries(); 
W_0->Scale(1/norm0); // normalize 

Double_t norm2 = W_2->GetEntries(); 
W_2->Scale(1/norm2); // normalize 

Double_t norm4 = W_4->GetEntries(); 
W_4->Scale(1/norm4); // normalize
} 

// fit data AC with templates 
TF1 *fit = new TF1("fit",fitFcn,1,52,3);
fit->SetParameters((1-a_2-a_4),a_2,a_4);
fit->SetParNames("x","y","z"); 
hsto->Fit(fit,"0","",1,52); //******!*!*!*!!*!*!*!!*!!*!!*!!*

double A0 = fit->GetParameter(0);
double a2 = fit->GetParameter(1);
double a4 = fit->GetParameter(2);
double A0Err = fit->GetParError(0);
double a2Err = fit->GetParError(1);
double a4Err = fit->GetParError(2);
/*
double a2 = y/(x+y+z);
double a4 = z/(x+y+z);
double error = 10*sqrt(((pow(xErr,2)*pow(y,2)) + (pow(x,2)*pow(yErr,2)) + (2*x*pow(yErr,2)*z) + 
	(pow(yErr,2)*pow(z,2)) + (pow(y,2)*pow(zErr,2)))/pow(x + y + z,4));
*/
double ChiSqr = fit->GetChisquare();
double NDF = fit->GetNDF();
double X2 = ChiSqr/NDF;

W_0->Scale(1-a2-a4);
W_2->Scale(a2);
W_4->Scale(a4);

W_0->Add(W_2);
W_0->Add(W_4);

W_0->SetMarkerColor(kRed);
W_0->SetMarkerStyle(4);
W_0->SetMarkerSize(2);
W_0->Draw("same");

cout << "X2/NDF = " << X2 << endl;
cout << "a2 = "<< a2 << endl;
cout << "a4 = "<< a4 << endl;


}