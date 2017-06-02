// used to find calculated a2,a4 and c2,c4 coefficients 




// global variables 
Double_t x0,x2,x4,y00,y2,y4,z0,z2,z4;

// define angular correlation function
Double_t acFnc(Double_t *x, Double_t *p)
{
	Double_t val = p[0]*(1 + p[1]*::ROOT::Math::legendre(2,x[0]) + p[2]*::ROOT::Math::legendre(4,x[0]));
	return val;
}


// define functions for Z approximations 

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


Double_t C2(double c2, double c4, double x2, double y2, double z2)
{
	Double_t out = x2 - c2*(x2 - y2) - c4*(x2 - z2);
	return out;
}

Double_t C2error(double c2, double c4, double x2, double y2, double z2, double x2error, double y2error, double z2error)
{
	Double_t out = sqrt(pow(1 - c2 - c4,2)*pow(x2error,2) + pow(c2,2)*pow(y2error,2) + 
    pow(c4,2)*pow(z2error,2));
    return out;
}

Double_t C4(double c2, double c4, double x4, double y4, double z4)
{
  Double_t out = x4 - c2*(x4 - y4) - c4*(x4 - z4);
  return out;
}

Double_t C4error(double c2, double c4, double x4, double y4, double z4, double x4error, double y4error, double z4error)
{
  Double_t out = sqrt(pow(1 - c2 - c4,2)*pow(x4error,2) + pow(c2,2)*pow(y4error,2) + 
    pow(c4,2)*pow(z4error,2));
  return out;
}



void approxDis(){

// define initial fit functions
TF1 *l_fit_0 = new TF1("l_fit_0",acFnc,-1,1,3);
TF1 *l_fit_2 = new TF1("l_fit_2",acFnc,-1,1,3);
TF1 *l_fit_4 = new TF1("l_fit_4",acFnc,-1,1,3);

// load template files 
// ultra_aX.root -> Z distributions for 0th, 2nd, and 4th order Legendre polynomials
TFile* file_a0 = new TFile("ultra_a0Nu.root");
TFile* file_a2 = new TFile("ultra_a2Nu.root");
TFile* file_a4 = new TFile("ultra_a4Nu.root");

// extract gamma-gamma matrix
TDirectory* NDdir_a0 = file_a0->GetDirectory("GriffinND");
TDirectory* NDdir_a2 = file_a2->GetDirectory("GriffinND");
TDirectory* NDdir_a4 = file_a4->GetDirectory("GriffinND");

THnSparse* thn_a0 = (THnSparse*) NDdir_a0->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");
THnSparse* thn_a2 = (THnSparse*) NDdir_a2->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");
THnSparse* thn_a4 = (THnSparse*) NDdir_a4->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse_0res");


// extract and create angular index histogram for each template 
TAngularCorrelation* ac = new TAngularCorrelation();
ac->GenerateMaps(16,110);

TH2D* h2_a0 = ac->Create2DSlice(thn_a0,1173,1174,kFALSE,kFALSE);
h2_a0->SetName("h2_a0");
TH2D* h2_a2 = ac->Create2DSlice(thn_a2,1173,1174,kFALSE,kFALSE);
h2_a2->SetName("h2_a2");
TH2D* h2_a4 = ac->Create2DSlice(thn_a4,1173,1174,kFALSE,kFALSE);
h2_a4->SetName("h2_a4");

TH1D* h1_a0 = ac->IntegralSlices(h2_a0,1332,1333);
TH1D* h1_a2 = ac->IntegralSlices(h2_a2,1332,1333);
TH1D* h1_a4 = ac->IntegralSlices(h2_a4,1332,1333);

TH1D* Z_0 = ac->DivideByWeights(h1_a0,false,false); 
Double_t norm0 = Z_0->GetEntries(); 
Z_0->Scale(1/norm0); // normalize 

TH1D* Z_2 = ac->DivideByWeights(h1_a2,false,false); 
Double_t norm2 = Z_2->GetEntries(); 
Z_2->Scale(1/norm2); // normalize 

TH1D* Z_4 = ac->DivideByWeights(h1_a4,false,false); 
Double_t norm4 = Z_4->GetEntries(); 
Z_4->Scale(1/norm4); // normalize


//create angular correlation for each template 
TGraphAsymmErrors* W_0 = ac->CreateGraphFromHst(Z_0,false,false); 
TGraphAsymmErrors* W_2 = ac->CreateGraphFromHst(Z_2,false,false); 
TGraphAsymmErrors* W_4 = ac->CreateGraphFromHst(Z_4,false,false); 

// pull out approximated coefficients 𝛼,𝛽,𝛾
// 𝛼 = x
// 𝛽 = y
// 𝛾 = z

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





// load data + create angular correlations + fit + extract coefficients 
int n = 16;
TGraphAsymmErrors* acgraphs[n];

// pre-determined angular correlation coefficients for each constructed AC
double a2[16] = {0.1020408163265306,0.1044025523625067-0.25,0.24725274725274723,0.3863636363636361,0.41601255886970173,
0.17857142857142844,0.4336734693877551,0.3571428571428571,0.9383116883116882,0.7499999999999996,0.9272727272727279,1.0227272727272727,
0.7386363636363631,0.4220779220779217,0.42857142857142855,0.31249999999999967};
double a4[16] = {0.009070294784580489,0.011377043820133971,-0.02930402930402929,0.008582326764144938,-0.17725132010846295,
-0.004329004329004331,-0.1836734693877548,1.1428571428571415,0.364135864135864,0.045454545454545504,0.377622377622378,0.5967365967365967,
0.08815426997245186,0.1847041847041847,-0.09523809523809529,0.002066115702479341};

Double_t a2out[n], a4out[n], a2error_out[n], a4error_out[n];
Double_t c2out[n], c4out[n], c2error_out[n], c4error_out[n];
Double_t c0[n], c2[n], c4[n], c0error[n], c2error[n], c4error[n];


for (int j=0;j<n;j++){

TH1D* Z0 = (TH1D*)Z_0->Clone("Z0");
TH1D* Z2 = (TH1D*)Z_2->Clone("Z2");
TH1D* Z4 = (TH1D*)Z_4->Clone("Z4");

Z0->Scale(1 - a_2[j] - a_4[j]);
Z2->Scale(a_2[j]);
Z4->Scale(a_4[j]);

Z0->Add(Z2);
Z0->Add(Z4);

acgraphs[j] = ac->CreateGraphFromHst(Z0,false,false); //angular correlation 
acgraphs[j]->SetName(Form("angular_correlation%i",j));

TF1* l_fitData = new TF1("l_fitData",acFnc,-1,1,3);
l_fitData->SetParameters(1,a_2[j],a_4[j]);
acgraphs[j]->Fit(l_fitData,"","",-1,1);
c0[j] = l_fitData->GetParameter(0);
c2[j] = abs(l_fitData->GetParameter(1));
c4[j] = l_fitData->GetParameter(2);
c0error[j] = l_fitData->GetParError(0);
c2error[j] = l_fitData->GetParError(1);	
c4error[j] = l_fitData->GetParError(2);



c2out[j] = C2(a2[j],a4[j],x2,y2,z2); 
c4out[j] = C4(a2[j],a4[j],x4,y4,z4);
c2error_out[j] = C2error(a2[j],a4[j],x2,y2,z2,x2error,y2error,z2error);
c4error_out[j] = C4error(a2[j],a4[j],x4,y4,z4,x4error,y4error,z4error);


a2out[j] = A2(c2[j],c4[j],x2,x4,y2,y4,z2,z4);
a4out[j] = A4(c2[j],c4[j],x2,x4,y2,y4,z2,z4);
a2error_out[j] = A2error(c2[j],c4[j],x2,x4,y2,y4,z2,z4,x2error,x4error,y2error,y4error,z2error,z4error);
a4error_out[j] = A4error(c2[j],c4[j],x2,x4,y2,y4,z2,z4,x2error,x4error,y2error,y4error,z2error,z4error);


if(Z0) delete Z0;
if(Z2) delete Z2;
if(Z4) delete Z4;

}


bool drawPlota2 = false;
bool drawPlotc4 = true;



// plot a2 v c2 (or a4 v c4)
gStyle->SetOptFit(000);
gStyle->SetOptStat(0);

if(drawPlota2){
TGraphErrors* points = new TGraphErrors(n,c2,c2out,0,c2error);
TF1 *fita2 = new TF1("fita2","[0]*x",0,1);
points->SetTitle(
	";"
	"c_{2};"
	"#tilde{c}_{2};"
	);
points->GetXaxis()->SetRangeUser(0,1.5);
points->Fit("fita2","0");
//points->Fit("pol1","");
points->SetMarkerStyle(5);
points->SetMarkerColor(1);
points->Draw("ap");
TLine *line4 = new TLine(0.05,0.05,0.45,0.45);
line4->SetLineColor(kRed);
line4->SetLineStyle(2);
//line4->Draw("same");
double slopea2 = fita2->GetParameter(0);
TLine *fitres1 = new TLine(0.05,(slopea2*0.05),1.05,(slopea2*1.05));
fitres1->SetLineColor(38);
fitres1->SetLineStyle(2);
fitres1->Draw("same");

}




// plot or a4 v c4
if(drawPlotc4){
gStyle->SetOptFit(0);
gStyle->SetOptStat(0);
TGraphErrors* pointsa4 = new TGraphErrors(n,a4,a4out,0,a4error);
TF1 *fita4 = new TF1("fita4","[0]*x",0,1);
pointsa4->SetTitle(
	";"
	"c_{4};"
	"#tilde{c}_{4};"
	);
//pointsa4->Fit("pol1","");
//pointsa4->GetYaxis()->SetRangeUser(-0.4,0.8);
//pointsa4->GetXaxis()->SetDomainUser(-0.95,0.9);
pointsa4->Fit("fita4","0");
pointsa4->SetMarkerStyle(5);
pointsa4->Draw("ap");
TLine *line4 = new TLine(-0.16,-0.16,0.02,0.02);
line4->SetLineColor(kRed);
line4->SetLineStyle(2);
double slope = fita4->GetParameter(0);
TLine *fitres2 = new TLine(-0.21,(slope*-0.21),01.05,(slope*1.05));
fitres2->SetLineColor(46);
fitres2->SetLineStyle(2);
fitres2->Draw("same");
//line4->Draw("same");
}


}





