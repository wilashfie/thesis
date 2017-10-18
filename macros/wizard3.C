// master script for angular correation data analysis


// define global objects, histograms, etc. 
TH1D *Z0, *Z2, *Z4;
TH1D *W_0, *W_2, *W_4;
TH1D *hD, *hT;

TGraphAsymmErrors *acD, *acT;
TGraphAsymmErrors *Z0ac, *Z2ac, *Z4ac;

TF1* L;
TF1* L_temp;

double *covarmatrixZ, *covarmatrixAL;
double x0,x2,x4,y00,y2,y4,z0,z2,z4;
double x0error,x2error,x4error,y00error,y2error,y4error,z0error,z2error,z4error;
Int_t npfits;

/*
double a_2 = 0.1020408163265306;
double a_4 = 0.009070294784580489;

double a_2 = 0.3571428571428571;
double a_4 = 1.1428571428571415;

double a_2 = 0.17857142857142;
double a_4 = -0.76190476190;

double a_2 = 0.1044025523625067;
double a_4 = 0.011377043820133971;
*/

double a_2 = 0.1020408163265306;
double a_4 = 0.009070294784580489;



//--------------------------------------------------------------------------------------------------------------
// functions for calculating physical a2,a4 coeffients 
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

//--------------------------------------------------------------------------------------------------------------
// fit functions

Double_t ACfcn(Double_t *x, Double_t *p)
{
  Double_t val = p[0]*(1 + p[1]*::ROOT::Math::legendre(2,x[0]) + p[2]*::ROOT::Math::legendre(4,x[0]));
  return val;
}

Double_t P2(Double_t *x)
{
  Double_t val = ::ROOT::Math::legendre(2,x[0]);
  return val;
}

Double_t P4(Double_t *x)
{
  Double_t val = ::ROOT::Math::legendre(4,x[0]);
  return val;
}

// define linear combination fit 
Double_t Zfit(Double_t *x, Double_t *p){
  Double_t i = x[0];
  Int_t bin_0 = Z0->GetXaxis()->FindBin(i);
  Double_t z_0 = (1-p[1]-p[2])*Z0->GetBinContent(bin_0);
  Int_t bin_2 = Z2->GetXaxis()->FindBin(i);
  Double_t z_2 = p[1]*Z2->GetBinContent(bin_2);
  Int_t bin_4 = Z4->GetXaxis()->FindBin(i);
  Double_t z_4 = p[2]*Z4->GetBinContent(bin_4);
  Double_t out = p[0]*(z_0 + z_2 + z_4);
  return out;
}

// define albegraic fit 
double Lfit(Double_t *x, Double_t *p){
  Double_t l0 = (1-p[1]-p[2])*(x0*(1 + x2*P2(x) + x4*P4(x)));
  Double_t l2 = p[1]*(y00*(1 + y2*P2(x) + y4*P4(x)));
  Double_t l4 = p[2]*(z0*(1 + z2*P2(x) + z4*P4(x)));
  Double_t out = p[0]*(l0 + l2 + l4);
  return out;
}

void ZFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  TAxis *xaxis1  = hD->GetXaxis();
  TAxis *yaxis1  = hD->GetYaxis();

  int nbinX1 = hD->GetNbinsX(); 
  int nbinY1 = hD->GetNbinsY(); 

  double chi2 = 0; 
  double npfits = 0; // number of points fit (used for NDF)
  double content,data_error,func_val,error; 
  double x[1];
  double tmp;
  int bin;
  double z_0,z_2,z_4;
  double z_0error,z_2error,z_4error,Zerror;
  double dA0, dZ0, dZ2, dZ4;

  // iterating through first histogram, x-axis
  for (int i = 1; i <= nbinX1; ++i) { 

    bin = xaxis1->FindBin(i);

    // calculate error in Z distribution
    z_0 = Z0->GetBinContent(bin);
    z_2 = Z2->GetBinContent(bin);
    z_4 = Z4->GetBinContent(bin);
    z_0error = Z0->GetBinError(bin);
    z_2error = Z2->GetBinError(bin);
    z_4error = Z4->GetBinError(bin);
    dA0 = (1-p[1]-p[2])*z_0 + p[1]*z_2 + p[2]*z_4;
    dZ0 = p[0]*(1-p[1]-p[2]);
    dZ2 = p[0]*p[1];
    dZ4 = p[0]*p[2];
    Zerror = sqrt(pow(dZ0*z_0error,2)+pow(dZ2*z_2error,2)+pow(dZ4*z_4error,2));

    // chi-square analysis 
    x[0] = xaxis1->GetBinCenter(bin);
    content = hD->GetBinContent(bin);
    if (content<=0) continue;
    data_error = hD->GetBinError(bin);
    if (data_error<=0) continue;
    func_val = Zfit(x,p);
    error = sqrt(pow(data_error,2) + pow(Zerror,2));
    tmp = (content-func_val)/error;
    chi2 += tmp*tmp;
  }
  fval = chi2; // final value
}


void LFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
 
  double Nbins = acD->GetN();
  double x[1],y[1];
  double chi2 = 0; 
  double npfits = 0; // number of points fit (used for NDF)
  double content,data_error,func_val,error; 
  double tmp;
  int bin;


  // iterating through first histogram, x-axis
  for (int i = 0; i < Nbins; ++i) { 

    acD->GetPoint(i,x[0],y[0]);

    // calculate error in L distribution 
    double Lerror = sqrt(pow(P2(x),2)*pow(x0,2)*pow(x2error,2) + 
      pow(x0error,2)*pow(1 + P2(x)*x2 + P4(x)*x4,2) + 
      pow(P4(x),2)*pow(x0,2)*pow(x4error,2) + 
      pow(P2(x),2)*pow(y00,2)*pow(y2error,2) + 
      pow(y00error,2)*pow(1 + P2(x)*y2 + P4(x)*y4,2) + 
      pow(P4(x),2)*pow(y00,2)*pow(y4error,2) + 
      pow(P2(x),2)*pow(z0,2)*pow(z2error,2) + 
      pow(z0error,2)*pow(1 + P2(x)*z2 + P4(x)*z4,2) + 
      pow(P4(x),2)*pow(z0,2)*pow(z4error,2));

    // chi-square analysis 

    // get data bin error
    data_error = acD->GetErrorY(i);
    if (data_error<=0) continue;
    func_val = Lfit(x,p);

    error = sqrt(pow(data_error,2) + pow(Lerror,2));

    tmp = (y[0]-func_val)/data_error;
    chi2 += tmp*tmp;
    npfits++;
  }
  fval = chi2; // final value
}


Double_t fcnA2(double c2, double c4, double x0, double x2, double x4, double y0, double y2, double y4, double z0, double z2, double z4)
{
  Double_t out = (x0*z0*(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4)))/
   ((x0*x2 + c2*(-x0 + y0) - y0*y2)*(x0*x4 + c4*(-x0 + z0) - z0*z4)*
     (1 - ((x0*x4 + c4*(-x0 + y0) - y0*y4)*(x0*x2 + c2*(-x0 + z0) - z0*z2))/
        ((x0*x2 + c2*(-x0 + y0) - y0*y2)*(x0*x4 + c4*(-x0 + z0) - z0*z4))));
  return out;
}

Double_t fcnA4(double c2, double c4, double x0, double x2, double x4, double y0, double y2, double y4, double z0, double z2, double z4)
{
  Double_t out = (x0*y0*(-(x4*y2) + c4*(-x2 + y2) + c2*(x4 - y4) + x2*y4))/
   ((x0*x2 + c2*(-x0 + y0) - y0*y2)*(x0*x4 + c4*(-x0 + z0) - z0*z4)*
     (1 - ((x0*x4 + c4*(-x0 + y0) - y0*y4)*(x0*x2 + c2*(-x0 + z0) - z0*z2))/
        ((x0*x2 + c2*(-x0 + y0) - y0*y2)*(x0*x4 + c4*(-x0 + z0) - z0*z4))));
  return out;
}

Double_t fcnA2error(double c2, double c4, double x0, double x2, double x4, double y0, double y2, double y4, double z0, double z2, double z4, double c2error, double c4error, double x0error, double x2error, double x4error, double y0error, double y2error, double y4error, double z0error, double z2error, double z4error)
{
  Double_t out = sqrt((pow(x0,2)*pow(x4error,2)*pow(y0,2)*pow(z0,2)*
       pow(x0*x2 + c2*(-x0 + z0) - z0*z2,2)*
       pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2) + 
      pow(c2error,2)*pow(x0,2)*pow(y0,2)*pow(z0,2)*
       pow(x0*x2 + c2*(-x0 + z0) - z0*z2,2)*
       pow(x4*y2 - x2*y4 - x4*z2 + y4*z2 + x2*z4 - y2*z4,2) + 
      pow(x0,2)*pow(y0,2)*pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*
       pow(z0,2)*pow(z2error,2)*pow(x0*x4 + c4*(-x0 + z0) - z0*z4,2) + 
      pow(x0,2)*pow(x2error,2)*pow(y0,2)*pow(z0,2)*
       pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2)*
       pow(x0*x4 + c4*(-x0 + z0) - z0*z4,2) + 
      pow(c2error,2)*pow(x0,2)*pow(y0,2)*pow(z0,2)*
       pow(x4*y2 - x2*y4 - x4*z2 + y4*z2 + x2*z4 - y2*z4,2)*
       pow(x0*x4 + c4*(-x0 + z0) - z0*z4,2) + 
      pow(x0,4)*pow(y0,2)*pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*
       pow(z0error,2)*pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2) + 
      pow(x0,2)*pow(y0,2)*pow(y4error,2)*pow(z0,2)*
       pow(x0*x2 + c2*(-x0 + z0) - z0*z2,2)*
       pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2) + 
      pow(x0error,2)*pow(y0,2)*pow(z0,4)*
       pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2)*
       pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2) + 
      pow(x0,2)*pow(y0,2)*pow(y2error,2)*pow(z0,2)*
       pow(x0*x4 + c4*(-x0 + z0) - z0*z4,2)*
       pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2) + 
      pow(x0,2)*pow(y0error,2)*pow(z0,2)*
       pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2)*
       pow(x0*x4*y2 - x0*x2*y4 + c4*(x0*(x2 - y2) + z0*(y2 - z2)) + y4*z0*z2 - 
         y2*z0*z4 + c2*(-(x0*x4) + x0*y4 - y4*z0 + z0*z4),2) + 
      pow(x0,2)*pow(y0,2)*pow(z0,2)*pow(x0*x2 + c2*(-x0 + z0) - z0*z2,2)*
       pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2)*pow(z4error,2))/
    pow(x0*x4*y0*y2 - x0*x2*y0*y4 - x0*x4*z0*z2 + y0*y4*z0*z2 + 
      c4*(y0*z0*(y2 - z2) + x0*(x2*y0 - y0*y2 - x2*z0 + z0*z2)) + x0*x2*z0*z4 - 
      y0*y2*z0*z4 + c2*(y0*z0*(-y4 + z4) + x0*(-(x4*y0) + y0*y4 + x4*z0 - z0*z4)),4));
  return out;
}

Double_t fcnA4error(double c2, double c4, double x0, double x2, double x4, double y0, double y2, double y4, double z0, double z2, double z4, double c2error, double c4error, double x0error, double x2error, double x4error, double y0error, double y2error, double y4error, double z0error, double z2error, double z4error)
{
  Double_t out = sqrt((pow(x0,2)*pow(y0,2)*pow(x0*x4 + c4*(-x0 + y0) - y0*y4,2)*
       pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*pow(z0,2)*
       pow(z2error,2) + pow(x0,2)*pow(x4error,2)*pow(y0,2)*
       pow(x0*x2 + c2*(-x0 + y0) - y0*y2,2)*pow(z0,2)*
       pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2) + 
      pow(x0,2)*pow(x2error,2)*pow(y0,2)*pow(x0*x4 + c4*(-x0 + y0) - y0*y4,2)*
       pow(z0,2)*pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2) + 
      pow(x0error,2)*pow(y0,4)*
       pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*pow(z0,2)*
       pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2) + 
      pow(c2error,2)*pow(x0,2)*pow(y0,2)*pow(x0*x2 + c2*(-x0 + y0) - y0*y2,2)*
       pow(z0,2)*pow(x4*y2 - x2*y4 - x4*z2 + y4*z2 + x2*z4 - y2*z4,2) + 
      pow(c2error,2)*pow(x0,2)*pow(y0,2)*pow(x0*x4 + c4*(-x0 + y0) - y0*y4,2)*
       pow(z0,2)*pow(x4*y2 - x2*y4 - x4*z2 + y4*z2 + x2*z4 - y2*z4,2) + 
      pow(x0,2)*pow(y0,2)*pow(y2error,2)*pow(x0*x4 + c4*(-x0 + y0) - y0*y4,2)*
       pow(z0,2)*pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2) + 
      pow(x0,4)*pow(y0error,2)*
       pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*pow(z0,2)*
       pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2) + 
      pow(x0,2)*pow(y0,2)*pow(x0*x2 + c2*(-x0 + y0) - y0*y2,2)*pow(y4error,2)*
       pow(z0,2)*pow(c4*(x2 - z2) + x4*z2 - x2*z4 + c2*(-x4 + z4),2) + 
      pow(x0,2)*pow(y0,2)*pow(c4*(x2 - y2) + x4*y2 - x2*y4 + c2*(-x4 + y4),2)*
       pow(z0error,2)*pow(x0*x4*z2 - y0*y4*z2 + 
         c4*(x0*(x2 - z2) + y0*(-y2 + z2)) - x0*x2*z4 + y0*y2*z4 + 
         c2*(-(x0*x4) + y0*y4 + x0*z4 - y0*z4),2) + 
      pow(x0,2)*pow(y0,2)*pow(x0*x2 + c2*(-x0 + y0) - y0*y2,2)*pow(z0,2)*
       pow(c4*y2 - c2*y4 - c4*z2 + y4*z2 + c2*z4 - y2*z4,2)*pow(z4error,2))/
    pow(x0*x4*y0*y2 - x0*x2*y0*y4 - x0*x4*z0*z2 + y0*y4*z0*z2 + 
      c4*(y0*z0*(y2 - z2) + x0*(x2*y0 - y0*y2 - x2*z0 + z0*z2)) + x0*x2*z0*z4 - 
      y0*y2*z0*z4 + c2*(y0*z0*(-y4 + z4) + x0*(-(x4*y0) + y0*y4 + x4*z0 - z0*z4)),4));
  return out;
}


void wizard3(){

TAngularCorrelation* ac = new TAngularCorrelation();
ac->GenerateMaps(16,110);

// load and create AC from data 
TFile* data_file = new TFile("Converted.root");
TDirectory* NDdir_data = data_file->GetDirectory("Griffin3D");
TH3I* hst = (TH3I*) NDdir_data->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry");
hst->GetXaxis()->SetRangeUser(1173, 1174);
TH2D* h2_data = (TH2D*) hst->Project3D("zyeo");
TH1D* h1_data = ac->IntegralSlices(h2_data,1332,1333);
hD = ac->DivideByWeights(h1_data,false,false); //angular index
acD = ac->CreateGraphFromHst(hD,false,false); //angular correlation  

// get simulations made from collect.C
TFile* sim_file = new TFile("simulations.root");
hT = (TH1D*) sim_file->Get("hT");
Z0 = (TH1D*) sim_file->Get("Z0");
Z2 = (TH1D*) sim_file->Get("Z2");
Z4 = (TH1D*) sim_file->Get("Z4");
acT = (TGraphAsymmErrors*) sim_file->Get("acT");
Z0ac = (TGraphAsymmErrors*) sim_file->Get("Z0ac");
Z2ac = (TGraphAsymmErrors*) sim_file->Get("Z2ac");
Z4ac = (TGraphAsymmErrors*) sim_file->Get("Z4ac");

// fit funcitons for Z distributions 
TF1 *l_fit_0 = new TF1("l_fit_0",ACfcn,-1,1,3);
TF1 *l_fit_2 = new TF1("l_fit_2",ACfcn,-1,1,3);
TF1 *l_fit_4 = new TF1("l_fit_4",ACfcn,-1,1,3);

//fit Z0
l_fit_0->SetParameters(1,0,0);
Z0ac->Fit(l_fit_0,"","",-1,1);
x0 = l_fit_0->GetParameter(0);
x2 = l_fit_0->GetParameter(1);
x4 = l_fit_0->GetParameter(2);
x0error = l_fit_0->GetParError(0);
x2error = l_fit_0->GetParError(1);
x4error = l_fit_0->GetParError(2);

//fit Z2 
l_fit_2->SetParameters(1,1,0); 
Z2ac->Fit(l_fit_2,"","",-1,1);
y00 = l_fit_2->GetParameter(0);
y2 = l_fit_2->GetParameter(1);
y4 = l_fit_2->GetParameter(2);
y00error = l_fit_2->GetParError(0);
y2error = l_fit_2->GetParError(1);
y4error = l_fit_2->GetParError(2);

//fit Z4
l_fit_4->SetParameters(1,0,1);
Z4ac->Fit(l_fit_4,"","",-1,1);
z0 = l_fit_4->GetParameter(0);
z2 = l_fit_4->GetParameter(1);
z4 = l_fit_4->GetParameter(2);
z0error = l_fit_4->GetParError(0);
z2error = l_fit_4->GetParError(1);
z4error = l_fit_4->GetParError(2);


// open file 
TFile *file = new TFile("plots.root","RECREATE");


// ---------------------------------------------------------------------
// data graph (x = cos(theta))

TF1 *bare_fit = new TF1("bare_fit",ACfcn,-1,1,3);

acD->SetTitle(
  "Bare Fit to Data;"
  "Cos(#theta);"
  "Counts;"
  );

bare_fit->SetParameters(1,a_2,a_4);
bare_fit->SetParNames("c0","c2","c4"); 

TFitResultPtr fit_result = acD->Fit(bare_fit,"0S","",-1,1);

double c0 = bare_fit->GetParameter(0);
double c2 = bare_fit->GetParameter(1);
double c4 = bare_fit->GetParameter(2);
double c0error = bare_fit->GetParError(0);
double c2error = bare_fit->GetParError(1);
double c4error = bare_fit->GetParError(2);

TF1 *bare_fit_draw = new TF1("bare_fit_draw",ACfcn,-1,1,3);
bare_fit_draw->FixParameter(0,c0);
bare_fit_draw->FixParameter(1,c2);
bare_fit_draw->FixParameter(2,c4);

TMatrixDSym covB(3);
covB = fit_result->GetCovarianceMatrix();

TCanvas* data_ac_canvas = new TCanvas("data_ac");
data_ac_canvas->cd();

acD->Draw("AP");

bare_fit_draw->SetLineColor(kRed);
bare_fit_draw->Draw("same");

TPaveText* pt_dataAC = new TPaveText(0.05,175,0.95,180);
pt_dataAC->SetTextFont(133);
pt_dataAC->SetTextSize(20);
pt_dataAC->AddText(Form("c_{2} = %f #pm %f",c2,c2error));
pt_dataAC->AddText(Form("c_{4} = %f #pm %f",c4,c4error));
pt_dataAC->Draw();

data_ac_canvas->Write();


// ---------------------------------------------------------------------
// data-to-template comparison 


TH1D* hDnorm = (TH1D*) hD->Clone();
hDnorm->Sumw2();
hDnorm->Scale(100000/hD->Integral());

TH1D *hTnorm = (TH1D*) hT->Clone();
hTnorm->Sumw2();
hTnorm->Scale(100000/hT->Integral());


double Chi=0;
  for(int i=1;i<52;i++){
    int tempbin = hTnorm->GetXaxis()->FindBin(i);
    int databin = hDnorm->GetXaxis()->FindBin(i);

    double temp = hTnorm->GetBinContent(tempbin);
    double tempError = hTnorm->GetBinError(tempbin);
    double data = hDnorm->GetBinContent(tempbin);
    double dataError = hDnorm->GetBinError(tempbin);
    double error = sqrt(pow(dataError,2) + pow(tempError,2));
    double chi2 = pow((data-temp)/error,2);
    Chi = Chi + chi2;
  }

cout << "Chi = " << Chi <<endl;

TCanvas* data_temp_canvas = new TCanvas("data_temp_comparison");
data_temp_canvas->cd();

TGraphAsymmErrors* acDnorm = ac->CreateGraphFromHst(hDnorm,false,false);
TGraphAsymmErrors* acTnorm = ac->CreateGraphFromHst(hTnorm,false,false);

TMultiGraph* multi_datatemp = new TMultiGraph();

acDnorm->SetMarkerStyle(8);
multi_datatemp->Add(acDnorm);

acTnorm->SetMarkerColor(kRed);
acTnorm->SetLineColor(kRed);
acTnorm->SetMarkerStyle(4);
multi_datatemp->Add(acTnorm);

multi_datatemp->SetTitle(
  "data-to-template comparison;"
  "Cos(#theta);"
  "Counts;"
  );
multi_datatemp->Draw("ap");

TPaveText* pt_data_temp = new TPaveText(0.25,2110,0.75,2160);
pt_data_temp->AddText(Form("#chi^{2}/#nu = %f",Chi/48));
pt_data_temp->Draw();

data_temp_canvas->Write();





// ---------------------------------------------------------------------
// best fit - Z distribution

// how many parameters does the fit have?
static const int nfitpars = 3;

// create fit function 
TF1 *RSfit = new TF1("RSfit",Zfit,1,52,nfitpars);
double iniParams[nfitpars] = {1,1,1};
RSfit->SetParameters(iniParams);
RSfit->SetParNames("A_{0}","a_{2}","a_{4}"); 

// initialize the fitter
TVirtualFitter *minuitZ = TVirtualFitter::Fitter(0,nfitpars); // arguments here are TObject *obj and maximum number of parameters
// initialize the parameters
for (int i = 0; i < nfitpars; ++i) {  
  minuitZ->SetParameter(i, RSfit->GetParName(i), RSfit->GetParameter(i), 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
}

// set function to minimize
minuitZ->SetFCN(ZFcn);

double arglist[100];
// set print level
arglist[0] = 0;
minuitZ->ExecuteCommand("SET PRINT",arglist,2);

// set minimization parameters
arglist[0] = 5000; // number of function calls
arglist[1] = 0.001; // tolerance

// execute minimization of function
// list of minimization algorithms in https://root.cern.ch/fitting
// here's a good website for minuit: http://iminuit.readthedocs.org/en/latest/api.html
minuitZ->ExecuteCommand("MIGRAD",arglist,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments
minuitZ->ExecuteCommand("MINOS",arglist,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments

// get result and errors
double minParamsZ[nfitpars];
double parErrorsZ[nfitpars];
for (int i = 0; i < nfitpars; ++i) {  
  minParamsZ[i] = minuitZ->GetParameter(i);
  parErrorsZ[i] = minuitZ->GetParError(i);
}

// get statistics about fit
double Zchi2, edm, errdef; 
int nvpar, nparx;
minuitZ->GetStats(Zchi2,edm,errdef,nvpar,nparx);
covarmatrixZ = minuitZ->GetCovarianceMatrix();

// assign values to function
RSfit->SetParameters(minParamsZ);
RSfit->SetParErrors(parErrorsZ);
RSfit->SetChisquare(Zchi2);
int ndf = npfits-nvpar;
RSfit->SetNDF(ndf);

double A0 = minParamsZ[0];
double a2 = minParamsZ[1];
double a4 = minParamsZ[2];


// find reduced chi-square for Z distribution
double combined_error[53] = {0};

for(int i=1;i<52;i++){
  int bin = Z0->FindBin(i);

  double z0 = Z0->GetBinContent(bin);
  double z2 = Z2->GetBinContent(bin);
  double z4 = Z4->GetBinContent(bin);

  double z0error = Z0->GetBinError(bin);
  double z2error = Z2->GetBinError(bin);
  double z4error = Z4->GetBinError(bin);

  double dA0 = (1-a_2-a_4)*z0 + a_2*z2 + a_4*z4;
  double dZ0 = A0*(1-a_2-a_4);
  double dZ2 = A0*a_2;
  double dZ4 = A0*a_4;

  combined_error[i] = sqrt(pow(dZ0*z0error,2)+pow(dZ2*z2error,2)+pow(dZ4*z4error,2));
}


TH1D* Z = (TH1D*) Z0->Clone();
TH1D* Z_2 = (TH1D*) Z2->Clone();
TH1D* Z_4 = (TH1D*) Z4->Clone();

Z->Scale(1-a2-a4);
Z->Add(Z_2,a2);
Z->Add(Z_4,a4);
Z->Scale(A0);

double X2 = 0;
double diff[52];
double index[52];
double error[52];

for(int i=1;i<52;i++){
  int data_bin = hD->GetXaxis()->FindBin(i);
  int Z_bin = Z->GetXaxis()->FindBin(i);

  double data = hD->GetBinContent(data_bin);
  double data_error = hD->GetBinError(data_bin);
  double z = Z->GetBinContent(Z_bin);
  double z_error = combined_error[i]; 

  error[i] = sqrt((data_error*data_error) + (z_error*z_error));
  index[i] = i;
  diff[i] = data - z;

  double chi2 = pow(diff[i]/error[i],2);
  X2 = X2 + chi2;
}

double reducedZ = X2/48;

TCanvas* Z_fit_canvas = new TCanvas("Z_fit");
Z_fit_canvas->cd();

TGraphAsymmErrors *Zac = ac->CreateGraphFromHst(Z,false,false);

TMultiGraph* multi_Zplot = new TMultiGraph();
acD->SetMarkerStyle(8);
multi_Zplot->Add(acD);
Zac->SetMarkerStyle(4);
Zac->SetMarkerColor(kRed);
Zac->SetLineColor(kRed);
multi_Zplot->Add(Zac);
multi_Zplot->SetTitle(
  "Best Z fit;"
  "Cos(#theta);"
  "Counts;"
  );
multi_Zplot->Draw("ap");

TPaveText* pt_Z = new TPaveText(0,172,0.75,177);
pt_Z->SetTextFont(133);
pt_Z->SetTextSize(20);
pt_Z->AddText(Form("a_{2} = %f #pm %f",a2,parErrorsZ[1]));
pt_Z->AddText(Form("a_{4} = %f #pm %f",a4,parErrorsZ[2]));
pt_Z->AddText(Form("#chi^{2}/#nu = %f",reducedZ));
pt_Z->Draw();

Z_fit_canvas->Write();


// residual plot
TCanvas* Z_res_canvas = new TCanvas("Z_res");
Z_res_canvas->cd();

double cosine[52];
for (int i=1;i<=52;i++)
{
  double angle = ac->GetAngleFromIndex(i-1);
  cosine[i] = TMath::Cos(angle); // maybe deg->rad needed
}

TGraphErrors* Zres = new TGraphErrors(51,cosine,diff,0,error);
Zres->SetMarkerStyle(8);
Zres->SetTitle(
  "Z residual;"
  "Cos(#theta);"
  "Counts;"
  );
Zres->Draw("ap");

TLine *zero = new TLine(-1,0,1,0);
zero->SetLineColor(kBlue);
zero->SetLineStyle(2);
zero->Draw("same");

Z_res_canvas->Write();


// ---------------------------------------------------------------------
// best fit - L distribution

// create fit function 
TF1 *ALfit = new TF1("ALfit",Lfit,-1,1,nfitpars);
ALfit->SetParameters(iniParams);
ALfit->SetParNames("A_{0}","a_{2}","a_{4}"); 

// initialize the fitter
TVirtualFitter *minuitAL = TVirtualFitter::Fitter(0,nfitpars); // arguments here are TObject *obj and maximum number of parameters
// initialize the parameters
for (int i = 0; i < nfitpars; ++i) {  
  minuitAL->SetParameter(i, ALfit->GetParName(i), ALfit->GetParameter(i), 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
}

// set function to minimize
minuitAL->SetFCN(LFcn);

double arglistAL[100];
// set print level
arglistAL[0] = 0;
minuitAL->ExecuteCommand("SET PRINT",arglistAL,2);

// set minimization parameters
arglistAL[0] = 5000; // number of function calls
arglistAL[1] = 0.001; // tolerance

// execute minimization of function
// list of minimization algorithms in https://root.cern.ch/fitting
// here's a good website for minuit: http://iminuit.readthedocs.org/en/latest/api.html
minuitAL->ExecuteCommand("MIGRAD",arglistAL,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments
minuitAL->ExecuteCommand("MINOS",arglistAL,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments

// get result and errors
double minParamsAL[nfitpars];
double parErrorsAL[nfitpars];
for (int i = 0; i < nfitpars; ++i) {  
  minParamsAL[i] = minuitAL->GetParameter(i);
  parErrorsAL[i] = minuitAL->GetParError(i);
}

// get statistics about fit
double chi2AL, edmAL, errdefAL; 
int nvparAL, nparxAL;
minuitAL->GetStats(chi2AL,edmAL,errdefAL,nvparAL,nparxAL);
covarmatrixAL = minuitAL->GetCovarianceMatrix();

// assign values to function
ALfit->SetParameters(minParamsAL);
ALfit->SetParErrors(parErrorsAL);
ALfit->SetChisquare(chi2AL);
int ndfAL = npfits-nvparAL;
ALfit->SetNDF(ndfAL);

double A2 = minParamsAL[1];
double A4 = minParamsAL[2];

// find a2,a4 via algebraic expressions
double A2test = fcnA2(c2,c4,x0,x2,x4,y00,y2,y4,z0,z2,z4);
double A4test = fcnA4(c2,c4,x0,x2,x4,y00,y2,y4,z0,z2,z4);
double A2error = fcnA2error(c2,c4,x0,x2,x4,y00,y2,y4,z0,z2,z4,c2error,c4error,x0error,x2error,x4error,y00error,y2error,y4error,z0error,z2error,z4error);
double A4error = fcnA4error(c2,c4,x0,x2,x4,y00,y2,y4,z0,z2,z4,c2error,c4error,x0error,x2error,x4error,y00error,y2error,y4error,z0error,z2error,z4error);


L = new TF1("L",Lfit,-1,1,3);
L->FixParameter(0,minParamsAL[0]);
L->FixParameter(1,A2);
L->FixParameter(2,A4);

double Nbins = acD->GetN();
double x[1],y[1];
double chi2 = 0; 
double npfits = 0; // number of points fit (used for NDF)
double content,func_val; 
double tmp;
int bin;

double c_error[52], L_diff[52], data_error[52];


// iterating through first histogram, x-axis
for (int i = 0; i < Nbins; ++i) { 

  acD->GetPoint(i,x[0],y[0]);

  double X = x[0];

  // calculate error in L distribution 
  double Lerror = sqrt(pow(P2(x),2)*pow(x0,2)*pow(x2error,2) + 
  pow(x0error,2)*pow(1 + P2(x)*x2 + P4(x)*x4,2) + 
  pow(P4(x),2)*pow(x0,2)*pow(x4error,2) + 
  pow(P2(x),2)*pow(y00,2)*pow(y2error,2) + 
  pow(y00error,2)*pow(1 + P2(x)*y2 + P4(x)*y4,2) + 
  pow(P4(x),2)*pow(y00,2)*pow(y4error,2) + 
  pow(P2(x),2)*pow(z0,2)*pow(z2error,2) + 
  pow(z0error,2)*pow(1 + P2(x)*z2 + P4(x)*z4,2) + 
  pow(P4(x),2)*pow(z0,2)*pow(z4error,2));

  // chi-square analysis 

  // get data bin error
  data_error[i] = acD->GetErrorY(i);
  if (data_error<=0) continue;
  func_val = L->Eval(X);

  c_error[i] = sqrt(pow(data_error[i],2) + pow(Lerror,2));
  L_diff[i] = y[0] - func_val;

  tmp = (L_diff[i])/c_error[i];
  chi2 += tmp*tmp;
  npfits++;
}


double reducedL = chi2/48;

TCanvas *L_fit_canvas = new TCanvas("L_fit");
L_fit_canvas->cd();

acD->SetTitle(
  "Best L fit;"
  "Cos(#theta);"
  "Counts;"
  );
acD->Draw("AP");

L->SetLineColor(kBlue);
L->SetLineColor(kBlue);
L->Draw("same");

TPaveText* pt_L = new TPaveText(0,174,0.75,181);
pt_L->SetTextFont(133);
pt_L->SetTextSize(20);
pt_L->AddText(Form("a_{2} = %f #pm %f",A2,parErrorsAL[1]));
pt_L->AddText(Form("a_{4} = %f #pm %f",A4,parErrorsAL[2]));
pt_L->AddText(Form("#chi^{2}/#nu = %f",reducedL));
pt_L->Draw();

L_fit_canvas->Write();


// residual plot
TCanvas* L_res_canvas = new TCanvas("L_res");
L_res_canvas->cd();

TGraphErrors* Lres = new TGraphErrors(51,cosine,L_diff,0,c_error);
Lres->SetMarkerStyle(8);
Lres->SetTitle(
  "L residual;"
  "Cos(#theta);"
  "Residual;"
  );
Lres->Draw("ap");

zero->Draw("same");

L_res_canvas->Write();



// ---------------------------------------------------------------------
// error ellipse plot

// bare fit
TMatrixD matrixB(2,2);
matrixB[0][0] = covB[1][1];
matrixB[0][1] = covB[1][2];
matrixB[1][0] = covB[2][1];
matrixB[1][1] = covB[2][2];
TMatrixDEigen eigenmatrixB(matrixB);
TVectorD eigenvalB = eigenmatrixB.GetEigenValuesRe();
TMatrixD eigenvecB = eigenmatrixB.GetEigenVectors();


double r1_B = sqrt(eigenvalB[0]);
double r2_B = sqrt(eigenvalB[1]);
double maxrB = TMath::Max(r1_B,r2_B);
double thetaB = TMath::RadToDeg()*atan(eigenvecB[0][1]/eigenvecB[0][0]);

TEllipse *eB = new TEllipse(c2,c4,r1_B*sqrt(6.18),r2_B*sqrt(6.18),0,360,thetaB);

// Z fit 
TMatrixD matrixZ(2,2);
matrixZ[0][0] = covarmatrixZ[4]; //why?
matrixZ[0][1] = covarmatrixZ[5];
matrixZ[1][0] = covarmatrixZ[7];
matrixZ[1][1] = covarmatrixZ[8];

TMatrixDEigen eigenmatrixZ(matrixZ);
TVectorD eigenvalZ = eigenmatrixZ.GetEigenValuesRe();
TMatrixD eigenvecZ = eigenmatrixZ.GetEigenVectors();

double r1_Z = sqrt(eigenvalZ[0]);
double r2_Z = sqrt(eigenvalZ[1]);
double maxrZ = TMath::Max(r1_Z,r2_Z);
double thetaZ = TMath::RadToDeg()*atan(eigenvecZ[0][1]/eigenvecZ[0][0]);

TEllipse *eZ = new TEllipse(a2,a4,r1_Z*sqrt(6.18),r2_Z*sqrt(6.18),0,360,thetaZ);

// L fit
TMatrixD matrixAL(2,2);
matrixAL[0][0] = covarmatrixAL[4];
matrixAL[0][1] = covarmatrixAL[5];
matrixAL[1][0] = covarmatrixAL[7];
matrixAL[1][1] = covarmatrixAL[8];
TMatrixDEigen eigenmatrixAL(matrixAL);
TVectorD eigenvalAL = eigenmatrixAL.GetEigenValuesRe();
TMatrixD eigenvecAL = eigenmatrixAL.GetEigenVectors();

double r1_AL = sqrt(eigenvalAL[0]);
double r2_AL = sqrt(eigenvalAL[1]);
double maxrAL = TMath::Max(r1_AL,r2_AL);
double thetaAL = TMath::RadToDeg()*atan(eigenvecAL[0][1]/eigenvecAL[0][0]);

TEllipse *eAL = new TEllipse(A2,A4,r1_AL*sqrt(6.18),r2_AL*sqrt(6.18),0,360,thetaAL);


cout << "testZ = " << covarmatrixAL[4] << endl;
cout << "testL = " << covarmatrixZ[4] << endl;


TCanvas *error_ellipse_canvas = new TCanvas("error_ellipse");
error_ellipse_canvas->cd();

TGraph* g = new TGraph(2);
g->SetTitle(";a_{2};a_{4}");
g->SetPoint(0,a2-maxrB*sqrt(35),a4-maxrB*sqrt(35));
g->SetPoint(1,a2+maxrB*sqrt(35),a4+maxrB*sqrt(35));
g->Draw("pa");

g->GetXaxis()->SetRangeUser(a2-0.015,a2+0.008);
g->GetYaxis()->SetRangeUser(a4-0.08,a4+0.08);

TMarker* physical = new TMarker(a_2,a_4,0);
physical->SetMarkerStyle(29);
physical->Draw();

eB->SetLineColor(kRed);
eB->Draw("same");
TMarker* bare = new TMarker(c2,c4,0);
bare->SetMarkerStyle(21);
bare->Draw("same");

eZ->SetLineColor(kBlue);
eZ->SetFillColorAlpha(kWhite,0);
eZ->Draw("");
TMarker* RS = new TMarker(a2,a4,0);
RS->SetMarkerStyle(20);
RS->Draw("");

eAL->SetLineColor(kGreen);
eAL->SetFillColorAlpha(kWhite,0);
eAL->Draw("");
TMarker* AL = new TMarker(A2,A4,0);
AL->SetMarkerStyle(22);
AL->Draw("");

error_ellipse_canvas->Write();


// ---------------------------------------------------------------------
// reduced chi-square v delta for Z and L distributions 

int n = 2000;
double a2_calc[n], a4_calc[n], delta[n], X2pearsons[n], X2combined[n], X2algebra[n], chi[n], a2_fit[n], a4_fit[n];

int ngraphs = 5;
TMultiGraph* mg[ngraphs];
TCanvas* c_chi[ngraphs];


for (int i=0;i<=4;++i){

// input information for selected cascade below:	
double ji = i; //i 
double jo = 2;
double jf = 0;
double diff = abs(jo-ji);
double L1 = 1;
double L1p = 2;
double L2 = 2;
double L2p = 2;
//double delta1 = 0;
double delta2 = 0;

double Ji = fabs(ji); 
double Jo = fabs(jo);
double Jf = fabs(jf);

// get scaling factor from template angular index 
TF1 *fit_scale = new TF1("fit_scale",Zfit,1,52,3);
fit_scale->SetParameters((1-a_2-a_4),a_2,a_4);
hT->Fit(fit_scale,"0q","",1,52);

double A0 = fit_scale->GetParameter(0);
double A0error = fit_scale->GetParError(0);

// vary mixing ratio
for(int i=0;i<n;i++){

  if(W_0) delete W_0;
  if(W_2) delete W_2;
  if(W_4) delete W_4;
  if(L_temp) delete L_temp;


  W_0 = (TH1D*) Z0->Clone();
  W_2 = (TH1D*) Z2->Clone();
  W_4 = (TH1D*) Z4->Clone();

  double pi = 3.1415926535;
  double t = atan(pi);
  double mix = tan(t);

  double value = tan(-pi/2 + i*pi/(n+1));

  a2_calc[i] = B(2,Jo,Ji,L1,L1p,value)*A(2,Jf,Jo,L2,L2p,delta2); //(-10+i)*0.1
  a4_calc[i] = B(4,Jo,Ji,L1,L1p,value)*A(4,Jf,Jo,L2,L2p,delta2);

  delta[i] = atan(value);


  // root chi-square test (pearsons)
  TF1 *fit_chi = new TF1("fit_chi",Zfit,1,52,3);
  fit_chi->FixParameter(0,A0);
  fit_chi->FixParameter(1,a2_calc[i]);
  fit_chi->FixParameter(2,a4_calc[i]);
  fit_chi->SetParNames("a0","a2","a4"); 
  hT->Fit(fit_chi,"q0","",1,52);

  double ChiSqr = fit_chi->GetChisquare();
  double NDF = fit_chi->GetNDF();
  X2pearsons[i] = ChiSqr/NDF;


  // chi-square using Z-dist.
  // calculate combined error
  double zError[53] = {0};
  for(int i=1;i<52;i++){
    int bin = W_0->FindBin(i);

    double z0 = W_0->GetBinContent(bin);
    double z2 = W_2->GetBinContent(bin);
    double z4 = W_4->GetBinContent(bin);
    double z0error = W_0->GetBinError(bin);
    double z2error = W_2->GetBinError(bin);
    double z4error = W_4->GetBinError(bin);

    double dA0 = (1-a2_calc[i]-a4_calc[i])*z0 + a2_calc[i]*z2 + a4_calc[i]*z4;
    double dZ0 = A0*(1-a2_calc[i]-a4_calc[i]);
    double dZ2 = A0*a2_calc[i];
    double dZ4 = A0*a4_calc[i];

    zError[i] = sqrt(pow(dA0*A0error,2)+pow(dZ0*z0error,2)+pow(dZ2*z2error,2)+pow(dZ4*z4error,2));
  }

  // construct AC using Z distribution 
  W_0->Scale(1-a2_calc[i]-a4_calc[i]);
  W_0->Add(W_2,a2_calc[i]);
  W_0->Add(W_4,a4_calc[i]);
  W_0->Scale(A0);

  double chiSquare=0;
  for(int i=1;i<52;i++){
    int tempbin = hT->GetXaxis()->FindBin(i);
    int zbin = W_0->GetXaxis()->FindBin(i);

    double temp = hT->GetBinContent(tempbin);
    double tempError = hT->GetBinError(tempbin);
    double z = W_0->GetBinContent(zbin);
    double zerror = zError[i];
    double error = pow(zerror,2) + pow(tempError,2);
    double chi = pow((z-temp),2)/error;

    chiSquare = chiSquare + chi;
  }

  X2combined[i] = chiSquare/48;




  L_temp = new TF1("L",Lfit,-1,1,3);
  L_temp->FixParameter(0,minParamsAL[0]);
  L_temp->FixParameter(1,a2_calc[i]);
  L_temp->FixParameter(2,a4_calc[i]);

  double Nbins = acD->GetN();
  double x[1],y[1];
  double chi2 = 0; 
  double npfits = 0; // number of points fit (used for NDF)
  double content,func_val; 
  double tmp;
  int bin;
  double c_error[52], L_diff[52], data_error[52];

  // iterating through first histogram, x-axis
  for (int i = 0; i < Nbins; ++i) { 

    acD->GetPoint(i,x[0],y[0]);
    double X = x[0];

    // calculate error in L distribution 
    double Lerror = sqrt(pow(P2(x),2)*pow(x0,2)*pow(x2error,2) + 
    pow(x0error,2)*pow(1 + P2(x)*x2 + P4(x)*x4,2) + 
    pow(P4(x),2)*pow(x0,2)*pow(x4error,2) + 
    pow(P2(x),2)*pow(y00,2)*pow(y2error,2) + 
    pow(y00error,2)*pow(1 + P2(x)*y2 + P4(x)*y4,2) + 
    pow(P4(x),2)*pow(y00,2)*pow(y4error,2) + 
    pow(P2(x),2)*pow(z0,2)*pow(z2error,2) + 
    pow(z0error,2)*pow(1 + P2(x)*z2 + P4(x)*z4,2) + 
    pow(P4(x),2)*pow(z0,2)*pow(z4error,2));

    // chi-square analysis:
    // get data bin error
    data_error[i] = acD->GetErrorY(i);
    if (data_error<=0) continue;
    func_val = L_temp->Eval(X);

    c_error[i] = sqrt(pow(data_error[i],2) + pow(Lerror,2));
    L_diff[i] = y[0] - func_val;

    tmp = (L_diff[i])/c_error[i];
    chi2 += tmp*tmp;
    npfits++;
  }

  X2algebra[i] = chi2/48;


}

c_chi[i] = new TCanvas(Form("#chi^{2}/#nu vs. #delta for J_{i} = %d",i));
c_chi[i]->cd();
mg[i] = new TMultiGraph();
mg[i]->SetTitle(Form("#chi^{2}/#nu vs. #delta for J_{i} = %d",i));

TGraph* plot_pearsons = new TGraph(n,delta,X2pearsons);
plot_pearsons->SetMarkerStyle(6);
mg[i]->Add(plot_pearsons);

TGraph* plot_combined = new TGraph(n,delta,X2combined);
plot_combined->SetMarkerStyle(6);
plot_combined->SetMarkerColor(kBlue);
mg[i]->Add(plot_combined);

TGraph* plot_algebra = new TGraph(n,delta,X2algebra);
plot_algebra->SetMarkerStyle(6);
plot_algebra->SetMarkerColor(kCyan);
mg[i]->Add(plot_algebra);

mg[i]->Draw("ap");

c_chi[i]->Write();

}

file->Close();

}




