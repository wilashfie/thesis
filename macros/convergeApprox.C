// convergece test of calulated a2 vs. fitted a2 


#include "TRandom3.h"


// global variables 
Double_t x0,x2,x4,y00,y2,y4,z0,z2,z4;



// define Legendre Polynomial 
Double_t acFnc(Double_t *x, Double_t *p)
{
	Double_t val = p[0]*(1 + p[1]*::ROOT::Math::legendre(2,x[0]) + p[2]*::ROOT::Math::legendre(4,x[0]));
	return val;
}


// define functions for linear approx 

Double_t A2(double c2, double c4, double x2, double x4, double y2, double y4, double z2, double z4)
{
	Double_t out = (-c4 + x4 + ((c2 - x2)*(x4 - z4))/(x2 - z2))/(x4 - y4 - ((x2 - y2)*(x4 - z4))/(x2 - z2));
	return out;
}

Double_t A4(double c2, double c4, double x2, double x4, double y2, double y4, double z2, double z4)
{
	Double_t out = (-c4 - x4 + ((c2 - x2)*(x4 - y4))/(x2 - y2))/(x4 - ((x4 - y4)*(x2 - z2))/(x2 - y2) - z4);
	return out;
}

Double_t A2error(double c2, double c4, double x2, double x4, double y2, double y4, double z2, double z4, double x2error, double x4error, double y2error, double y4error, double z2error, double z4error)
{
  Double_t out = sqrt((pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*pow(z2error,2)*
       pow(x4 - z4,2))/pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(x2error,2)*pow(x4 - z4,2)*
       pow(-(y4*z2) + c4*(-y2 + z2) + c2*(y4 - z4) + y2*z4,2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(y4error,2)*pow(x2 - z2,2)*
       pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(y2error,2)*pow(x4 - z4,2)*
       pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(x4error,2)*pow(x2 - z2,2)*
       pow(c4*(y2 - z2) + y4*z2 - y2*z4 + c2*(-y4 + z4),2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*pow(x2 - z2,2)*
       pow(z4error,2))/pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4));
  return out;
}

Double_t A4error(double c2, double c4, double x2, double x4, double y2, double y4, double z2, double z4, double x2error, double x4error, double y2error, double y4error, double z2error, double z4error)
{
  Double_t out = sqrt((pow(x4 - y4,2)*pow(2*x2*x4 + c4*(x2 - y2) - x4*y2 - x2*y4 + c2*(-x4 + y4),2)*
       pow(z2error,2))/pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(x2error,2)*pow(x4 - y4,2)*
       pow(-(c2*y4) + c4*(y2 - z2) + 2*x4*(y2 - z2) + y4*z2 + c2*z4 - y2*z4,2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(y2error,2)*pow(x4 - y4,2)*
       pow(2*x2*x4 + c4*(x2 - z2) - x4*z2 - x2*z4 + c2*(-x4 + z4),2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(x2 - y2,2)*pow(y4error,2)*
       pow(2*x2*x4 + c4*(x2 - z2) - x4*z2 - x2*z4 + c2*(-x4 + z4),2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(x4error,2)*pow(x2 - y2,2)*
       pow(2*x2*y4 + c4*(y2 - z2) - y4*z2 - 2*x2*z4 + y2*z4 + c2*(-y4 + z4),2))/
     pow(x4*(y2 - z2) + y4*z2 - y2*z4 + x2*(-y4 + z4),4) + 
    (pow(c4 + x4 - ((c2 - x2)*(x4 - y4))/(x2 - y2),2)*pow(z4error,2))/
     pow(-x4 + ((x4 - y4)*(x2 - z2))/(x2 - y2) + z4,4));

  return out;
}



void convergeNu(){


// create angular correlation for data 
TDirectory* NDdir = gFile->GetDirectory("GriffinND");
THnSparse* my_thn = (THnSparse*) NDdir->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");

TAngularCorrelation* ac = new TAngularCorrelation();
ac->GenerateMaps(16,110);

TH2D* h2 = ac->Create2DSlice(my_thn,1173,1174,kFALSE,kFALSE);
TH1D* h1 = ac->IntegralSlices(h2,1332,1333);
TH1D* hsto = ac->DivideByWeights(h1,false,false); //angular index




// load template files 
TFile* file_a0 = new TFile("ultra_a0Nu.root");
TFile* file_a2 = new TFile("ultra_a2Nu.root");
TFile* file_a4 = new TFile("ultra_a4Nu.root");

TDirectory* NDdir_a0 = file_a0->GetDirectory("GriffinND");
TDirectory* NDdir_a2 = file_a2->GetDirectory("GriffinND");
TDirectory* NDdir_a4 = file_a4->GetDirectory("GriffinND");

THnSparse* thn_a0 = (THnSparse*) NDdir_a0->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");
THnSparse* thn_a2 = (THnSparse*) NDdir_a2->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");
THnSparse* thn_a4 = (THnSparse*) NDdir_a4->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");


// create angular index for each template 
TH2D* h2_a0 = ac->Create2DSlice(thn_a0,1173,1174,kFALSE,kFALSE);
h2_a0->SetName("h2_a0");
TH2D* h2_a2 = ac->Create2DSlice(thn_a2,1173,1174,kFALSE,kFALSE);
h2_a2->SetName("h2_a2");
TH2D* h2_a4 = ac->Create2DSlice(thn_a4,1173,1174,kFALSE,kFALSE);
h2_a4->SetName("h2_a4");

TH1D* h1_a0 = ac->IntegralSlices(h2_a0,1332,1333);
TH1D* h1_a2 = ac->IntegralSlices(h2_a2,1332,1333);
TH1D* h1_a4 = ac->IntegralSlices(h2_a4,1332,1333);

Z_0 = ac->DivideByWeights(h1_a0,false,false); 
Z_2 = ac->DivideByWeights(h1_a2,false,false); 
Z_4 = ac->DivideByWeights(h1_a4,false,false); 


//create angular correlation for each template 
TGraphAsymmErrors* W_0 = ac->CreateGraphFromHst(Z_0,false,false); 
TGraphAsymmErrors* W_2 = ac->CreateGraphFromHst(Z_2,false,false); 
TGraphAsymmErrors* W_4 = ac->CreateGraphFromHst(Z_4,false,false); 

// define initial fit functions
TF1 *l_fit_0 = new TF1("l_fit_0",acFnc,-1,1,3);
TF1 *l_fit_2 = new TF1("l_fit_2",acFnc,-1,1,3);
TF1 *l_fit_4 = new TF1("l_fit_4",acFnc,-1,1,3);

//fit W_0
l_fit_0->SetParameters(1,0,0);
W_0->Fit(l_fit_0,"","",-1,1);
double x0 = l_fit_0->GetParameter(0);
double x2 = l_fit_0->GetParameter(1);
double x4 = l_fit_0->GetParameter(2);
cout << "x0 = " << x0 << " x2 = " << x2 << " x4 = " << x4 << endl;
double x0error = l_fit_0->GetParError(0);
double x2error = l_fit_0->GetParError(1);
double x4error = l_fit_0->GetParError(2);
cout << "x0error = " << x0error << " x2error = " << x2error << " x4error = " << x4error << endl;

//fit W_2 
l_fit_2->SetParameters(1,1,0); 
W_2->Fit(l_fit_2,"","",-1,1);
double y00 = l_fit_2->GetParameter(0);
double y2 = l_fit_2->GetParameter(1);
double y4 = l_fit_2->GetParameter(2);
cout << "y00 = " << y00 << " y2 = " << y2 << " y4 = " << y4 << endl;
double y00error = l_fit_2->GetParError(0);
double y2error = l_fit_2->GetParError(1);
double y4error = l_fit_2->GetParError(2);
cout << "y0error = " << y00error << " y2error = " << y2error << " y4error = " << y4error << endl;

//fit W_4
l_fit_4->SetParameters(1,0,1);
W_4->Fit(l_fit_4,"","",-1,1);
double z0 = l_fit_4->GetParameter(0);
double z2 = l_fit_4->GetParameter(1);
double z4 = l_fit_4->GetParameter(2);
cout << "z0 = " << z0 << " z2 = " << z2 << " z4 = " << z4 << endl;
double z0error = l_fit_4->GetParError(0);
double z2error = l_fit_4->GetParError(1);
double z4error = l_fit_4->GetParError(2);
cout << "z0error = " << z0error << " z2error = " << z2error << " z4error = " << z4error << endl;

double a_2 = 0.1044025523625067;
double a_4 = 0.011377043820133971;

double a2 = A2(a_2,a_4,x2,x4,y2,y4,z2,z4); // !!!!!!
double a4 = A4(a_2,a_4,x2,x4,y2,y4,z2,z4);
double a2error = A2error(a_2,a_4,x2,x4,y2,y4,z2,z4,x2error,x4error,y2error,y4error,z2error,z4error);
double a4error = A4error(a_2,a_4,x2,x4,y2,y4,z2,z4,x2error,x4error,y2error,y4error,z2error,z4error);



// loop
int n = 20;
double c2[n], c4[n],c2error[n],c4error[n];
double outA2[n], outA4[n], outA2error[n], outA4error[n];
double N[n];

for(int i=0;i<n;i++){

if(gRandom) delete gRandom;
gRandom = new TRandom3(0);

// replicate experimental data
gRandom = new TRandom3(0);  
TH1D* sample = new TH1D("sample","sample distribution",52,1,52);
gRandom->SetSeed(0);
sample->FillRandom(hsto,400000*i);

TGraphAsymmErrors* angular_correlation = ac->CreateGraphFromHst(sample,false,false); //angular correlation

TF1 *legendre = new TF1("legendre",acFnc,-1,1,3);
legendre->SetParameters(1,1,1);
angular_correlation->Fit(legendre,"q0","",-1,1); 
c2[i] = legendre->GetParameter(1);
c4[i] = legendre->GetParameter(2);
c2error[i] = legendre->GetParError(1);
c4error[i] = legendre->GetParError(2);

N[i] = 400000*i;

outA2[i] = a2/c2[i];
outA2error[i] = outA2[i]*sqrt(pow(a2error/a2,2)+pow(c2error[i]/c2[i],2));

cout << "i = " << i << endl;
cout << "a2/c2 = " << outA2[i] << endl;
cout << "_______________" << endl;

if(sample) delete sample;

}

// plot
TGraphErrors* plot = new TGraphErrors(n,N,outA2,0,outA2error);
TF1 *f1 = new TF1("f1","[0]",N[0],N[n]);
plot->SetTitle(
  ";"
  "N Events;"
  "#frac{a_{2,calc}}{a_{2,fit}};"
  );
plot->SetMarkerStyle(1);
plot->SetLineColorAlpha(38,0.75);
plot->SetFillColor(38);
plot->SetMarkerColorAlpha(38,0.75);
plot->Draw("ap");
plot->Fit("f1","0");
double ymax4 = f1->GetParameter(0);
TLine *line4 = new TLine(N[0],ymax4,N[100],ymax4);
line4->SetLineColor(kRed);
line4->SetLineStyle(2);
line4->Draw("same");



}