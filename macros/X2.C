// runs X2 analysis for best fit of a2/a4 for a given cascade 


// functions for a2,a4
double Factorial(double val)
{
	double fac;
	if(val>1)
	{
		fac = val*Factorial(val-1);
	}
	else
	{
		fac=1;
	}
	return fac;
}


double ClebschGordan(double j1, double m1, double j2, double m2, double j, double m)
{
  // Conditions check
  if( 2*j1 != floor(2*j1) || 2*j2 != floor(2*j2) || 2*j != floor(2*j) || 2*m1 != floor(2*m1) || 2*m2 != floor(2*m2) || 2*m != floor(2*m) )
  {
  	cout << "All arguments must be integers or half-integers." << endl;
    return 0;
  }
  if(m1 + m2 != m){return 0;}
  if( j1 - m1 != floor ( j1 - m1 ) ){return 0;}
  if( j2 - m2 != floor ( j2 - m2 ) ){return 0;}
  if( j - m != floor ( j - m ) ){return 0;}
  if(j > j1 + j2 || j < fabs(j1 - j2)){return 0;}
  if(fabs(m1) > j1){return 0;}
  if(fabs(m2) > j2){return 0;}
  if(fabs(m) > j){return 0;}

  double term, cg;
  double term1 = pow((((2*j+1)/Factorial(j1+j2+j+1))*Factorial(j2+j-j1)*Factorial(j+j1-j2)*Factorial(j1+j2-j)*Factorial(j1+m1)*Factorial(j1-m1)*Factorial(j2+m2)*Factorial(j2-m2)*Factorial(j+m)*Factorial(j-m)),(0.5));
  double sum = 0;

  for(int k=0;k <= 99;k++)
  {
    if( (j1+j2-j-k < 0) || (j-j1-m2+k < 0) || (j-j2+m1+k < 0) || (j1-m1-k < 0) || (j2+m2-k < 0) )
    {
    }
    else
    {
      term = Factorial(j1+j2-j-k)*Factorial(j-j1-m2+k)*Factorial(j-j2+m1+k)*Factorial(j1-m1-k)*Factorial(j2+m2-k)*Factorial(k);
      if((k%2) == 1)
      {
        term = -1*term;
      }
      sum = sum + 1.0/term;
    }
  }
  cg = term1*sum;
  return cg;
}


double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3)
{
  // Conditions check
  if( 2*j1 != floor(2*j1) || 2*j2 != floor(2*j2) || 2*j3 != floor(2*j3) || 2*m1 != floor(2*m1) || 2*m2 != floor(2*m2) || 2*m3 != floor(2*m3) )
  {
    cout << "All arguments must be integers or half-integers." << endl;
    return 0;
  }
  if(m1 + m2 + m3 != 0)  {return 0;}
  if( j1 + j2 + j3 != floor(j1 + j2 + j3) )  {return 0;}
  if(j3 > j1 + j2 || j3 < fabs(j1 - j2))  {return 0;}
  if(fabs(m1) > j1)  {return 0;}
  if(fabs(m2) > j2)  {return 0;}
  if(fabs(m3) > j3)  {return 0;}

  double out = (pow((-1),(j1-j2-m3)))/(pow((2*j3+1),(1.0/2.0)))*ClebschGordan(j1,m1,j2,m2,j3,-1*m3);
  return out;
}


double Wigner6j(double J1, double J2, double J3, double J4, double J5, double J6)
{
  // Conditions check
  if(J3 > J1 + J2 || J3 < fabs(J1 - J2)){return 0;}
  if(J3 > J4 + J5 || J3 < fabs(J4 - J5)){return 0;}
  if(J6 > J2 + J4 || J6 < fabs(J2 - J4)){return 0;}
  if(J6 > J1 + J5 || J6 < fabs(J1 - J5)){return 0;}

  double j1 = J1;
  double j2 = J2;
  double j12 = J3;
  double j3 = J4;
  double j = J5;
  double j23 = J6;
  double sum = 0;

  for(double m1 = -j1 ; m1 <= j1 ; m1++ )
  {
    for(double m2 = -j2 ; m2 <= j2 ; m2++ )
    {
      for(double m3 = -j3 ; m3 <= j3 ; m3++ )
      {
        for(double m12 = -j12 ; m12 <= j12 ; m12++ )
        {
          for(double m23 = -j23 ; m23 <= j23 ; m23++ )
          {
            for(double m = -j ; m <= j ; m++ )
            {
              sum = sum + pow((-1),(j3+j+j23-m3-m-m23))*Wigner3j(j1,j2,j12,m1,m2,m12)*Wigner3j(j1,j,j23,m1,-m,m23)*Wigner3j(j3,j2,j23,m3,m2,-m23)*Wigner3j(j3,j,j12,-m3,m,m12);
            }
          }
        }
      }
    }
  }
  return sum;
}


double RacahW(double a, double b, double c, double d, double e, double f)
{
  double out = pow((-1),(a+b+d+c))*Wigner6j(a,b,e,d,c,f);
  return out;
}

double F(int k, double jf, int L1, int L2, double ji)
{
  double out;
  double CG = ClebschGordan(L1,1,L2,-1,k,0);
  if(CG == 0)
  {
    return 0;
  }
  double W = RacahW(ji,ji,L1,L2,k,jf);
  if(W == 0)
  {
    return 0;
  }
  out = pow((-1),(jf-ji-1))*(pow((2*L1+1)*(2*L2+1)*(2*ji+1),(1.0/2.0)))*CG*W;
  return out;
}

double B(int k, double ji, double jf, int L1, int L2, double delta)
{
  double out;
  out = (1/(1+pow(delta,2)))*(F(k,jf,L1,L1,ji)+(pow((-1),((double)(L1+L2))))*2*delta*F(k,jf,L1,L2,ji)+delta*delta*F(k,jf,L2,L2,ji) );
  return out;
}

double A(int k, double ji, double jf, int L1, int L2, double delta)
{
  double out;
  out = (1/(1+pow(delta,2)))*(F(k,ji,L1,L1,jf)+2*delta*F(k,ji,L1,L2,jf)+delta*delta*F(k,ji,L2,L2,jf));
  return out;
}

//--------------------------------------------------------------

// define global histograms 
TH1D *Z0, *Z2, *Z4;

Double_t acFnc(Double_t *x, Double_t *p)
{
  Double_t val = p[0]*(1 + p[1]*::ROOT::Math::legendre(2,x[0]) + p[2]*::ROOT::Math::legendre(4,x[0]));
  return val;
}

// define linear combination fit 
Double_t fitFcn(Double_t *x, Double_t *p){
  Double_t i = x[0];
  Int_t bin_0 = Z0->GetXaxis()->FindBin(i);
  Double_t z0 = p[0]*(Z0->GetBinContent(bin_0));
  Int_t bin_2 = Z2->GetXaxis()->FindBin(i);
  Double_t z2 = p[1]*(Z2->GetBinContent(bin_2));
  Int_t bin_4 = Z4->GetXaxis()->FindBin(i);
  Double_t z4 = p[2]*(Z4->GetBinContent(bin_4));
  return z0 + z2 + z4;
}



void X2(){

// option to normalize all histograms created
bool normalize = false;

// load simulated data for wanted cascade 
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
}
TGraphAsymmErrors* angular_correlation = ac->CreateGraphFromHst(hsto,false,false); //angular correlation 	

// get scaling factor a0
TF1 *scaleFit = new TF1("scaleFit",acFnc,-1,1,3);
scaleFit->SetParameters(1,0.1,0.1);
angular_correlation->Fit(scaleFit,"0","",-1,1);
double a0 = scaleFit->GetParameter(0);


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

TAngularCorrelation* acNu = new TAngularCorrelation();
acNu->GenerateMaps(16,110);

// create angular index for each template 
TH2D* h2_a0 = acNu->Create2DSlice(thn_a0,1173,1174,kFALSE,kFALSE);
h2_a0->SetName("h2_a0");
TH2D* h2_a2 = acNu->Create2DSlice(thn_a2,1173,1174,kFALSE,kFALSE);
h2_a2->SetName("h2_a2");
TH2D* h2_a4 = acNu->Create2DSlice(thn_a4,1173,1174,kFALSE,kFALSE);
h2_a4->SetName("h2_a4");

TH1D* h1_a0 = acNu->IntegralSlices(h2_a0,1332,1333);
TH1D* h1_a2 = acNu->IntegralSlices(h2_a2,1332,1333);
TH1D* h1_a4 = acNu->IntegralSlices(h2_a4,1332,1333);

Z0 = acNu->DivideByWeights(h1_a0,false,false); 
Z2 = acNu->DivideByWeights(h1_a2,false,false); 
Z4 = acNu->DivideByWeights(h1_a4,false,false); 


if(normalize){
Double_t norm0 = Z0->GetEntries(); 
Z0->Scale(1/norm0); // normalize 

Double_t norm2 = Z2->GetEntries(); 
Z2->Scale(1/norm2); // normalize 

Double_t norm4 = Z4->GetEntries(); 
Z4->Scale(1/norm4); // normalize 
}

//highest level in cascade is 2505.77 keV
//first gamma in cascade is 1173.24 keV
//second gamma in cascade is 1332.5 keV
// input information for selected cascade below:	
double ji = 1; 
double jo = 2;
double jf = 0;
double L1 = 2;
double L1p = 2;
double L2 = 2;
double L2p = 2;
//double delta1 = 0;
double delta2 = 0;

double Ji = fabs(ji); 
double Jo = fabs(jo);
double Jf = fabs(jf);

// # of iterations (or points)
int n = 10001;

double a2[n], a4[n], delta[n], X2[n], chi[n], x[n], y[n], z[n], a2_fit[n], a4_fit[n];

// vary mixing ratio
for(int i=0;i<n;i++){

double pi = 3.1415926535;
double t = atan(pi);
double mix = tan(t);
double value = tan(-pi/2 + i*pi/10000);

a2[i] = B(2,Jo,Ji,L1,L1p,value)*A(2,Jf,Jo,L2,L2p,delta2); 
a4[i] = B(4,Jo,Ji,L1,L1p,value)*A(4,Jf,Jo,L2,L2p,delta2);

TF1 *fit = new TF1("fit",fitFcn,1,52,3);
fit->SetParameter(0,(1-a2[i]-a4[i]));
fit->FixParameter(1,a2[i]);
fit->FixParameter(2,a4[i]);
fit->SetParNames("x","y","z"); 
hsto->Fit(fit,"q0","",1,52);

double ChiSqr = fit->GetChisquare();
double NDF = fit->GetNDF();

delta[i] = atan(value);
X2[i] = ChiSqr / NDF;
chi[i] = ChiSqr;

}

// fix for plot
double xmin=delta[0];
double xmax=delta[0];
for(int i=0;i<n;i++){
if (xmax<delta[(i+1)]){
  xmax=delta[(i+1)];
}
if (xmin>delta[i+1]){
  xmin=delta[(i+1)];
}}
xmin = xmin - 0.0005;
xmax = xmax + 0.0005;

// plot 
TGraph* plot = new TGraph(n,delta,X2);
plot->SetTitle(
	";"
	"arctan(#delta);"
	"#chi^{2}/NDF;"
	);
plot->SetMarkerStyle(6);
plot->Draw("ap");
TLine *line2 = new TLine(delta[0],2.7,delta[n-1],2.7);
line2->SetLineColor(kRed);
line2->SetLineStyle(2);
line2->Draw("same");



}