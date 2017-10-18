
// fits the Z distribution to individually simulated distributions

#include "TRandom3.h"

// define global histograms 
TH1D *W_0, *W_2, *W_4;

// option to normalize all histograms 
bool normalize = false;

// input for physical a2,a4 coefficients 
double a_2 = 0.3571428571428571;
double a_4 = 1.1428571428571415;

Double_t acFnc(Double_t *x, Double_t *p)
{
	Double_t val = p[0]*(1 + p[1]*::ROOT::Math::legendre(2,x[0]) + p[2]*::ROOT::Math::legendre(4,x[0]));
	return val;
}

// define linear combination fit 
Double_t fitFcn(Double_t *x, Double_t *p){
	Double_t i = x[0];
	Int_t bin_0 = W_0->GetXaxis()->FindBin(i);
	Double_t w0 = p[0]*W_0->GetBinContent(bin_0);
	Int_t bin_2 = W_2->GetXaxis()->FindBin(i);
	Double_t w2 = p[1]*W_2->GetBinContent(bin_2);
	Int_t bin_4 = W_4->GetXaxis()->FindBin(i);
	Double_t w4 = p[2]*W_4->GetBinContent(bin_4);
	return w0 + w2 + w4;
}

void convergeHisto(){

// create angular correlation from simulated data 
TDirectory* NDdir = gFile->GetDirectory("GriffinND");
THnSparse* my_thn = (THnSparse*) NDdir->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");

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
TGraphAsymmErrors* angular_correlation = ac->CreateGraphFromHst(hsto,false,false); //angular correlation

TF1 *legendre = new TF1("legendre",acFnc,-1,1,3);
legendre->SetParameters(1,1,1);
angular_correlation->Fit(legendre,"q0","",-1,1); 
double c2 = legendre->GetParameter(1);
double c4 = legendre->GetParameter(2);	



// load template files and create Z distribution
TFile* file_a0 = new TFile("ultra_a0Nu.root");
TFile* file_a2 = new TFile("ultra_a2Nu.root");
TFile* file_a4 = new TFile("ultra_a4Nu.root");

TDirectory* NDdir_a0 = file_a0->GetDirectory("GriffinND");
TDirectory* NDdir_a2 = file_a2->GetDirectory("GriffinND");
TDirectory* NDdir_a4 = file_a4->GetDirectory("GriffinND");

THnSparse* thn_a0 = (THnSparse*) NDdir_a0->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");
THnSparse* thn_a2 = (THnSparse*) NDdir_a2->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");
THnSparse* thn_a4 = (THnSparse*) NDdir_a4->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");

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

// start converge test 

int n = 501;
Double_t N[n], X2[n],a2[n],a2err[n],a4[n],out[n];

for(int i=0;i<n;i++){

// replicate experimental data
gRandom = new TRandom3(0);	
TH1D* sample = new TH1D("sample","sample distribution",52,1,52);
gRandom->SetSeed(0);
sample->FillRandom(hsto,50*i);

// fit data AC with templates 
TF1 *fit = new TF1("fit",fitFcn,1,52,3);
fit->SetParameters((1-a_2-a_4),a_2,a_4);
fit->SetParNames("x","y","z"); 
sample->Fit(fit,"","",1,52); 

double x = fit->GetParameter(0);
double y = fit->GetParameter(1);
double z = fit->GetParameter(2);

double xErr = fit->GetParError(0);
double yErr = fit->GetParError(1);
double zErr = fit->GetParError(2);

a2[i] = y/(x + y + z);

a2err[i] = 10*sqrt(((pow(xErr,2)*pow(y,2)) + (pow(x,2)*pow(yErr,2)) + (2*x*pow(yErr,2)*z) + 
      (pow(yErr,2)*pow(z,2)) + (pow(y,2)*pow(zErr,2)))/pow(x + y + z,4));

double ChiSqr = fit->GetChisquare();
double NDF = fit->GetNDF();

N[i] = 50*i;
X2[i] = ChiSqr/NDF;
out[i] = a2[i]/a_2;

if (sample) delete sample;

}

// plot
TGraphErrors* plot = new TGraphErrors(n,N,out,0,a2err);
TF1 *f1 = new TF1("f1","[0]",N[0],N[n]);
plot->SetTitle(
	";"
	"Number of Coincidences;"
	"#frac{a_{2,fit}}{a_{2,phys}};"
	);
plot->SetMarkerStyle(1);
plot->SetLineColorAlpha(38,0.55);
plot->SetFillColor(38);
plot->SetMarkerColorAlpha(38,0.55);
plot->Draw("ap");
plot->Fit("f1","0");
double ymax4 = f1->GetParameter(0);
TLine *line4 = new TLine(N[0],ymax4,N[100],ymax4);
line4->SetLineColor(kRed);
line4->SetLineStyle(2);
line4->Draw("same");


}