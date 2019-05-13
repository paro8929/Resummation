#include <finiteT2d-R3.h>


//Run parameters which you might want to change
const int N0=20;
const int N1=20;
const int mIT=20;
const int ITmax=40; //maximum number of iterations in root finder
const double Mmax=60; //maximal delta m2 in root finder
const double m02=1.0; //fiducial mass squared scale in units of lambda
double Temp=1.0; //temperature scale in units of lambda
const double cut=atan(1000); //maximum p2 where to cut off integrals
//

//don't change anything below unless
//you know what you're doing


int ABORT=0;

const double EulerGamma=0.5772156649015329;


double s2p0vals[N1],c2p0vals[N1],tp0vals[N1];
double tp1vals[N1];
double p0vals[N1];
double p0weights[N1];
double modp0weights[N1];
double wn[N0];
double wnweights[N0];

double **Sigma;
double **dPi,**ndPi;
vertex Gamma(N0,N1,2*N0,2*N1);
vertex nGamma(N0,N1,2*N0,2*N1);


double dx=0.1;
//double piprime=0; //derivative of pi @ Q^2=0

const int ITlen=20000;
double ITstep=0.01;
double ITx[ITlen];
double ITy[ITlen]; //Thermal integral data, discretized, from file
double JTy[ITlen]; //Thermal integral data, discretized, from file

#include <rootfinding.cpp>

void load_Bessel_data()
{
  
  
  char fname[255];
  sprintf(fname,"IT.dat");

  printf("===> Info: Loading Thermal Integral I_T(x) data from File %s\n",fname);

  fstream parameters;
  parameters.open(fname, ios::in);
  if(parameters.is_open())
    {
      for (int i=0;i<ITlen;i++)
	{
	  parameters >> ITx[i];
	  parameters >> ITy[i];	  
	}

      printf("finished\n");
    }

  else
    {
      printf("FAILED\n");
    }
  parameters.close();
}


void load_BesselJ_data()
{
  
  
  char fname[255];
  sprintf(fname,"JT.dat");

  printf("===> Info: Loading Thermal Integral J_T(x) data from File %s\n",fname);

  fstream parameters;
  parameters.open(fname, ios::in);
  if(parameters.is_open())
    {
      for (int i=0;i<ITlen;i++)
	{
	  parameters >> JTy[i];
	  parameters >> JTy[i];
	  //printf("i=%i j=%f\n",i,JTy[i]);
	}
      printf("finished\n");
    }

  else
    {
      printf("FAILED\n");
    }
  parameters.close();
}


void load_quadrature_data()
{
  

  char fname[255];
  sprintf(fname,"LegendreData/LegQuad-N%i.dat",N1-1);
  
  printf("===> Info: Loading Quadrature Data from File %s\n",fname);

  fstream parameters;
  parameters.open(fname, ios::in);
  if(parameters.is_open())
    {
      
      for (int i=0;i<N1;i++)
	{
	  parameters >> p0vals[i];
	  parameters >> p0weights[i];
	  //printf("i=%i (%.16f,%.16f)\n",i,p0vals[i],p0weights[i]);
	}

      printf("finished\n");
    }
  else
    {
      printf("FAILED\n");
    }
  parameters.close();
}


void derived_data()
{

  for (int i=0;i<N1;i++)
    {
      s2p0vals[i]=sin(p0vals[i]*M_PI*0.5)*sin(p0vals[i]*M_PI*0.5);
      c2p0vals[i]=cos(p0vals[i]*M_PI*0.5)*cos(p0vals[i]*M_PI*0.5);
      tp0vals[i]=tan(p0vals[i]*M_PI*0.5);
      modp0weights[i]=p0weights[i];
      tp1vals[i]=tan(p0vals[i]*cut);
    }
  modp0weights[0]=0.5*p0weights[0];

  for (int i=0;i<N0;i++)
    {
      wn[i]=2*M_PI*i*Temp;
      wnweights[i]=1;
    }
  wnweights[0]=0.5;

}

void allocate_memory()
{
  Sigma = new double*[N0];
  for (int i=0;i<N0;i++)
    Sigma[i]=new double[N1];
  
  dPi = new double*[N0];
  for (int i=0;i<N0;i++)
    dPi[i]=new double[N1];

  ndPi = new double*[N0];
  for (int i=0;i<N0;i++)
    ndPi[i]=new double[N1];

}

void free_memory()
{
  for (int i=0;i<N0;i++)
    {
      delete [] Sigma[i];
      delete [] dPi[i];
      delete [] ndPi[i];
    }
  
  delete [] Sigma;
  delete [] dPi;
  delete [] ndPi;
  
}

double max(double a,double b)
{
  double rr=0;
  if (a>=b)
    rr=a;
  else
    rr=b;

  return rr;
}

double newbigS(double p0,double p1,double shiftm)
{

  
  double res=0;

#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double s0=wn[i]-p0;
	double s1=tp0vals[j]-p1;
	double pp=int_field(wn[i],tp0vals[j],dPi,wn,tp1vals);
	double pi2=int_field(s0,s1,dPi,wn,tp1vals);

	double q0=wn[i];
	double q1=tp0vals[j];
	
	double ga1=Gamma.ipt4(q0,q1,s0,s1,wn,tp0vals);
	
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1;
	//printf("i=%i j=%i pp=%g pi2=%g res=%f\n",i,j,pp,pi2,res);
	
	s1=tp0vals[j]+p1;
	pi2=int_field(s0,s1,dPi,wn,tp1vals);
	ga1=Gamma.ipt4(q0,-q1,s0,-s1,wn,tp0vals);
	
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1;
	
	s0=wn[i]+p0;
	s1=tp0vals[j]-p1;
	pi2=int_field(s0,s1,dPi,wn,tp1vals);
	ga1=Gamma.ipt4(-q0,q1,-s0,s1,wn,tp0vals);
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1;
	
	s1=tp0vals[j]+p1;
	pi2=int_field(s0,s1,dPi,wn,tp1vals);
	ga1=Gamma.ipt4(-q0,-q1,-s0,-s1,wn,tp0vals);
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1;	
      }
  res/=4.;
  res*=2*Temp;
  return res;
}


double newbigSprime(double p0,double p1,double shiftm)
{

  
  double res=0;

#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double s0=wn[i]-p0;
	double s1=tp0vals[j]-p1;
	double pp=int_field(wn[i],tp0vals[j],dPi,wn,tp1vals);
	double pi2=int_field(s0,s1,dPi,wn,tp1vals);

	double q0=wn[i];
	double q1=tp0vals[j];
	
	double ga1=Gamma.ipt4(q0,q1,s0,s1,wn,tp0vals);	
	
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1*(ga1-1);
	//printf("i=%i j=%i pp=%g pi2=%g res=%f\n",i,j,pp,pi2,res);
	
	s1=tp0vals[j]+p1;
	pi2=int_field(s0,s1,dPi,wn,tp1vals);
	ga1=Gamma.ipt4(q0,-q1,s0,-s1,wn,tp0vals);
	
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1*(ga1-1);
	
	s0=wn[i]+p0;
	s1=tp0vals[j]-p1;
	pi2=int_field(s0,s1,dPi,wn,tp1vals);
	ga1=Gamma.ipt4(-q0,q1,-s0,s1,wn,tp0vals);
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1*(ga1-1);
	
	s1=tp0vals[j]+p1;
	pi2=int_field(s0,s1,dPi,wn,tp1vals);
	ga1=Gamma.ipt4(-q0,-q1,-s0,-s1,wn,tp0vals);
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2)*ga1*(ga1-1);	
      }
  res/=4.;
  res*=2*Temp;
  return res;
}



double highS(double p0,double p1,double shiftm)
{

  
  double res=0;

#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double s0=wn[i]-p0;
	double s1=tp0vals[j]-p1;
	double pp=0;
	double pi2=0;
	
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2);
	
	s1=tp0vals[j]+p1;

	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2);
	
	s0=wn[i]+p0;
	s1=tp0vals[j]-p1;
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2);
	
	s1=tp0vals[j]+p1;
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(s0*s0+s1*s1+m02+shiftm-pi2);	
      }
  res/=4.;
  res*=2*Temp;
  //res=0;
  return res;
}

void getall_sigma(double shiftm)
{
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
    {
      Sigma[i][j]=newbigS(wn[i],tp1vals[j],shiftm);
    }
}

double newdeltaPi(double p0,double p1,double shiftm)
{
  double res=0;

#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double s0=wn[i]-p0;
	double s1=tp0vals[j]-p1;
	
	double sigmaval=int_field(s0,s1,Sigma,wn,tp1vals);
	double tprop1=2*sigmaval/(1+2*sigmaval);
	double ga1=Gamma.ipt4(p0,p1,-wn[i],-tp0vals[j],wn,tp0vals);
	tprop1*=ga1;
	tprop1+=(1-ga1);
	
	s1=tp0vals[j]+p1;
	sigmaval=int_field(s0,s1,Sigma,wn,tp1vals);
	double tprop2=2*sigmaval/(1+2*sigmaval);
	ga1=Gamma.ipt4(p0,p1,-wn[i],tp0vals[j],wn,tp0vals);
	tprop2*=ga1;
	tprop2+=(1-ga1);

	
	s0=wn[i]+p0;
	s1=tp0vals[j]-p1;
	sigmaval=int_field(s0,s1,Sigma,wn,tp1vals);
	double tprop3=2*sigmaval/(1+2*sigmaval);
	ga1=Gamma.ipt4(p0,p1,wn[i],-tp0vals[j],wn,tp0vals);
	tprop3*=ga1;
	tprop3+=(1-ga1);
	
	s1=tp0vals[j]+p1;
	sigmaval=int_field(s0,s1,Sigma,wn,tp1vals);
	double tprop4=2*sigmaval/(1+2*sigmaval);
	ga1=Gamma.ipt4(p0,p1,wn[i],tp0vals[j],wn,tp0vals);
	tprop4*=ga1;
	tprop4+=(1-ga1);

	
	double pp=int_field(wn[i],tp0vals[j],dPi,wn,tp1vals);
	
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)*(tprop1+tprop2+tprop3+tprop4);
      }
  res*=2*Temp;
  return res;
}

void getall_dPi(double shiftm)
{
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
    {
      ndPi[i][j]=newdeltaPi(wn[i],tp1vals[j],shiftm);
    }
}

double newGamma(double k0,double k1,double p0,double p1,double shiftm)
{
  double res=0;
  

  //#pragma omp parallel for reduction( + : res)
    for (int i=0;i<N0;i++)
      for (int j=0;j<N1;j++)
	{
	  double q0=wn[i];
	  double q1=tp0vals[j];

	  
	  double s1,s2,t1,t2,pi1,pi2,g1,g2;
	  double ga1=1,ga2=1,ga3=1;
	  double reg1a=(q0+k0)*(q0+k0)+(q1+k1)*(q1+k1);
	  double arg1a=(q0+p0+k0)*(q0+p0+k0)+(q1+p1+k1)*(q1+p1+k1);
	  double arg2a=q0*q0+q1*q1;
	  
	  
	  s1=int_field(q0+k0,q1+k1,Sigma,wn,tp1vals);
	  t1=2./(1+2*s1);
	  
	  pi1=int_field(q0+p0+k0,q1+p1+k1,dPi,wn,tp1vals);
	  pi2=int_field(q0,q1,dPi,wn,tp1vals);	
	  g1=1./(arg1a+1+shiftm-pi1);
	  g2=1./(arg2a+1+shiftm-pi2);

	
	  ga1=Gamma.ipt4(q0,q1,k0,k1,wn,tp0vals);
	  ga2=Gamma.ipt4(p0,p1,-p0-k0-q0,-p1-k1-q1,wn,tp0vals);
	  ga3=Gamma.ipt4(p0+k0+q0,p1+k1+q1,-q0,-q1,wn,tp0vals);
	  
	  res+=wnweights[i]*modp0weights[j]/c2p0vals[j]*t1*g1*g2*ga1*ga2*ga3;

	  
	  q0=-wn[i];
	  q1=tp0vals[j];

	  reg1a=(q0+k0)*(q0+k0)+(q1+k1)*(q1+k1);
	  arg1a=(q0+p0+k0)*(q0+p0+k0)+(q1+p1+k1)*(q1+p1+k1);
	  arg2a=q0*q0+q1*q1;
	  	    
	    s1=int_field(q0+k0,q1+k1,Sigma,wn,tp1vals);
	    
	  
	  t1=2./(1+2*s1);
	  
	  pi1=int_field(q0+p0+k0,q1+p1+k1,dPi,wn,tp1vals);
	  pi2=int_field(q0,q1,dPi,wn,tp1vals);	
	  g1=1./(arg1a+1+shiftm-pi1);
	  g2=1./(arg2a+1+shiftm-pi2);
	  
	  ga1=Gamma.ipt4(q0,q1,k0,k1,wn,tp0vals);
	  ga2=Gamma.ipt4(p0,p1,-p0-k0-q0,-p1-k1-q1,wn,tp0vals);
	  ga3=Gamma.ipt4(p0+k0+q0,p1+k1+q1,-q0,-q1,wn,tp0vals);	  
	  res+=wnweights[i]*modp0weights[j]/c2p0vals[j]*t1*g1*g2*ga1*ga2*ga3;

	  q0=wn[i];
	  q1=-tp0vals[j];

	  reg1a=(q0+k0)*(q0+k0)+(q1+k1)*(q1+k1);
	  arg1a=(q0+p0+k0)*(q0+p0+k0)+(q1+p1+k1)*(q1+p1+k1);
	  arg2a=q0*q0+q1*q1;

	  s1=int_field(q0+k0,q1+k1,Sigma,wn,tp1vals);

	  
	  t1=2./(1+2*s1);
	  
	  pi1=int_field(q0+p0+k0,q1+p1+k1,dPi,wn,tp1vals);
	  pi2=int_field(q0,q1,dPi,wn,tp1vals);	
	  g1=1./(arg1a+1+shiftm-pi1);
	  g2=1./(arg2a+1+shiftm-pi2);

	  ga1=Gamma.ipt4(q0,q1,k0,k1,wn,tp0vals);
	  ga2=Gamma.ipt4(p0,p1,-p0-k0-q0,-p1-k1-q1,wn,tp0vals);
	  ga3=Gamma.ipt4(p0+k0+q0,p1+k1+q1,-q0,-q1,wn,tp0vals);
	  res+=wnweights[i]*modp0weights[j]/c2p0vals[j]*t1*g1*g2*ga1*ga2*ga3;

	  q0=-wn[i];
	  q1=-tp0vals[j];

	  reg1a=(q0+k0)*(q0+k0)+(q1+k1)*(q1+k1);
	  arg1a=(q0+p0+k0)*(q0+p0+k0)+(q1+p1+k1)*(q1+p1+k1);
	  arg2a=q0*q0+q1*q1;
	  
	  s1=int_field(q0+k0,q1+k1,Sigma,wn,tp1vals);

	  t1=2/(1+2*s1);
	  
	  pi1=int_field(q0+p0+k0,q1+p1+k1,dPi,wn,tp1vals);
	  pi2=int_field(q0,q1,dPi,wn,tp1vals);	
	  g1=1/(arg1a+1+shiftm-pi1);
	  g2=1/(arg2a+1+shiftm-pi2);

	  ga1=Gamma.ipt4(q0,q1,k0,k1,wn,tp0vals);
	  ga2=Gamma.ipt4(p0,p1,-p0-k0-q0,-p1-k1-q1,wn,tp0vals);
	  ga3=Gamma.ipt4(p0+k0+q0,p1+k1+q1,-q0,-q1,wn,tp0vals);
	  res+=wnweights[i]*modp0weights[j]/c2p0vals[j]*t1*g1*g2*ga1*ga2*ga3;
	  
      }
  
  res/=4.;
  res*=4*Temp;
  
  
  return res;
  
}

void Pgetall_gamma(double shiftm)
{

  int steppi=N0/10;

 
  #pragma omp parallel for 
  for (int k0=0;k0<N0;k0++)
    for (int k1=0;k1<N1;k1++)
      for (int p0=-N0+1;p0<N0;p0++)
	for (int p1=-N1+1;p1<N1;p1++)
	  {
	    double x0,x1;
	    
	    if (p0>=0)
	      x0=wn[p0];
	    else
	      x0=-wn[-p0];
	    
	    if (p1>=0)
	      x1=tp0vals[p1];
	    else
	      x1=-tp0vals[-p1];
	    
	    double val=newGamma(wn[k0],tp0vals[k1],x0,x1,shiftm);
	    nGamma.set(k0,k1,N0+p0,N1+p1,1-val);
	    

	  }

  
  for (int jj=0;jj<24;jj++)
    {
      double x0=0.401972;
      double x1=0.401972*cos(2*M_PI*jj/24.);
      double x2=0.401972*sin(2*M_PI*jj/24.);
      
    }
  
   
 
}


double int_field(double p0,double p1,double **field,double *xvals,double *yvals)
{
  double out=0;
  //assume symmetric arguments
  if (p0<0)
    p0*=-1;
  if (p1<0)
    p1*=-1;
  int a=0;
  int b=0;
  if (p0==xvals[0])
    {
      if (p1==yvals[0])
	{
	  out=field[0][0];
	}
      else if (p1<yvals[N1-1])
	{
	  while (p1>yvals[b])
	    b++;
	  out=linear_int(yvals[b-1],yvals[b],field[0][b-1],field[0][b],p1);
	}
      else
	{
	  out=0;
	}
    }
  else if (p0<xvals[N0-1])
    {
      while (p0>xvals[a])
	    a++;
      
      if (p1==yvals[0])
	{
	  out=linear_int(xvals[a-1],xvals[a],field[a-1][0],field[a][0],p0);
	}
      else if (p1<yvals[N1-1])
	{
	  while (p1>yvals[b])
	    b++;
	  double f0=field[a-1][b-1];
	  double f1=field[a][b-1];
	  double f2=field[a][b];
	  double f3=field[a-1][b];
	  
	  out=bilinear_int(xvals[a-1],xvals[a],yvals[b-1],yvals[b],f0,f1,f2,f3,p0,p1);			   
	}
      else
	{
	  out=0;
	}
    }
  else
    {
      out=0;
    }
  //printf("found case a=%i b=%i out=%f\n",a,b,out);
  
  return out;
}

double intIT(double x)
{
  double out=0;
  int a=1;
  while (x>ITx[a])
    a++;
  out=linear_int(ITx[a-1],ITx[a],ITy[a-1],ITy[a],x);
  return out;
}

double intJT(double x)
{
  double out=0;
  int a=0;
  while (x>ITx[a])
    a++;
  out=linear_int(ITx[a-1],ITx[a],JTy[a-1],JTy[a],x);
  return out;
}



double get_shiftm(double shiftm)
{
  double num=0;
  //double den=0;

  
#pragma omp parallel for reduction( + : num)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double pp=int_field(wn[i],tp0vals[j],dPi,wn,tp1vals);
	num+=pp*wnweights[i]*modp0weights[j]/c2p0vals[j]/(wn[i]*wn[i]+tp0vals[j]*tp0vals[j]+m02+shiftm)/(wn[i]*wn[i]+tp0vals[j]*tp0vals[j]+m02+shiftm-pp);	
      }
  num*=12*Temp;
  
  num+=12*intIT(Temp/sqrt(m02+shiftm));
  num+=3./M_PI*log(m02/(m02+shiftm));
  
  double res=num;
  return res;
}

double shiftfind(double *args)
{
  double sm=args[0];
  return get_shiftm(sm)-sm;
}


double altget_allshiftm(int verbose)
{
  if (verbose)
    printf("looking for root > %f\n",dPi[0][0]-m02+1.e-6);
  double Rag[6]={1,40,max(dPi[0][0]-m02+1.e-6,0),Mmax,double(verbose),0};
  double newshiftm=convergeroot2(&shiftfind,Rag);
  return newshiftm;
}

void exportfield(double **field,const char *name)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)      
      out << atan(wn[i]/(2*M_PI*Temp))*2/M_PI << "\t " << atan(tp1vals[j])*2/M_PI << "\t" << field[i][j] << "\n";
  
  out.close();
}

void exportfieldw0(double **field,const char *name)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int j=0;j<N1;j++)      
    out << atan(tp1vals[j]*tp1vals[j])*2/M_PI << "\t" << field[0][j] << "\n";
  
  out.close();
}

void exportSprime(const char *name,double shiftm)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int j=0;j<N1;j++)      
    out << atan(tp0vals[j]*tp0vals[j])*2/M_PI << "\t" << newbigSprime(0,tp0vals[j],shiftm) << "\n";
  
  out.close();
}


double altmeasureCP(double shiftm)
{
  double res=0;
  double ss;
#pragma omp parallel for reduction( + : res)
    for (int j=0;j<N1;j++)
      for (int i=0;i<N0;i++)
	{
	  
	  double pp;
	  pp=int_field(wn[i],tp0vals[j],dPi,wn,tp1vals);
	  
	  res+=wnweights[i]*modp0weights[j]/c2p0vals[j]*(log((tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp)/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm))+0.5*pp/(tp0vals[j]*tp0vals[j]+wn[i]*wn[i]+m02+shiftm-pp));
	  
      }
  
  res*=Temp;
  return res;
}

double altmeasureCS(double shiftm)
{
  double res=0;
#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double ss;
	
	ss=int_field(wn[i],tp0vals[j],Sigma,wn,tp1vals);
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]*(log(1+2*ss)-2*ss/(1+2*ss));
	
    }
  res*=Temp;
  return res;
}

double altmeasureextra(double shiftm)
{
  double res=0;
 
  for (int i=0;i<N0;i++)
    for (int j=0;j<N1;j++)
      {
	double ss;
	
	ss=int_field(wn[i],tp0vals[j],Sigma,wn,tp1vals);
	double pp=newbigSprime(wn[i],tp0vals[j],shiftm);
	res+=wnweights[i]*modp0weights[j]/c2p0vals[j]*2./(1+2*ss)*pp;
	
    }
  res*=Temp;
  return res;
}

double altmeasureVL(double shiftm)
{
  
  double vr=(m02+shiftm)*(m02+shiftm)/(48.);

  double logv=-(m02+shiftm)/24.*(m02+3/M_PI*(1+log(m02/(m02+shiftm)))); 
  return (vr+logv-(m02+shiftm)*intJT(Temp/sqrt(m02+shiftm)));


}

int main(int argc, char* argv[])
{
  load_quadrature_data();
  load_Bessel_data(); 
  load_BesselJ_data(); 
 
  
  derived_data();
  allocate_memory();

  
  char buffer[255];
  double shiftm=0.0;

  fstream hist;
  sprintf(buffer,"store/history-R3-NT%i-N%i-cut%i-m%f-T%f-IT%i.dat",N0,N1,int(tan(cut)),m02,Temp,mIT);
  hist.open(buffer,ios::out);


   for (int i=0;i<N0;i++)
     for (int j=0;j<N1;j++)
       {
	 Sigma[i][j]=0;
	 dPi[i][j]=0;
	 ndPi[i][j]=0;
       }
   Gamma.set_unity();

  

   double cp=0.5*altmeasureCP(shiftm);
   double cs=0.5*altmeasureCS(shiftm);
   double vl=altmeasureVL(shiftm);
   double extra=altmeasureextra(shiftm);
   
   printf("cs=%f\t cp=%f vl=%f extra=%f sum=%f\n",cs,cp,vl,extra,-cp-cs+vl-extra);
   
  for(int IT=1;IT<mIT;IT++)
    {
      printf("IT=%i\n",IT);
      getall_sigma(shiftm);
      printf("\t sigma done %f\n",int_field(wn[0],1.0,Sigma,wn,tp1vals));
      sprintf(buffer,"finTR3-Sdump%i",IT);        
      exportfieldw0(Sigma,buffer);

      
      
      Pgetall_gamma(shiftm);
      
      Gamma.addtx(&nGamma,0.5);
      

      

      
      getall_dPi(shiftm);
      printf("\t delta Pi done %f\n",int_field(wn[0],1.0,ndPi,wn,tp1vals));
      for (int i=0;i<N0;i++)
	for (int j=0;j<N1;j++)
	  dPi[i][j]=ndPi[i][j];

      
      sprintf(buffer,"finTR3-Pdump%i",IT);        
      exportfieldw0(dPi,buffer);
      
      shiftm=altget_allshiftm(0);
      printf("\t shiftm=%f\n",shiftm);
     

      sprintf(buffer,"finTR3-Sprimedump%i",IT);        
      exportSprime(buffer,shiftm);
      cp=0.5*altmeasureCP(shiftm);
      cs=0.5*altmeasureCS(shiftm);
      vl=altmeasureVL(shiftm);
      extra=-0.25*altmeasureextra(shiftm);
   
      printf("cs=%f\t cp=%f vl=%f extra=%f sum=%f\n",cs,cp,vl,extra,-cp-cs+vl-extra);
          
      sprintf(buffer,"%i\t %.3f\t %.3f\t %.10f\t 0\n",IT,shiftm,dPi[0][0],-cp-cs+vl-extra);
      hist << buffer;
     
      
      hist.flush();
      
    }
  
  hist.close();
  free_memory();
}
