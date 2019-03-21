#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <complex>
#include <cmath>
#include <omp.h>
#include <linear_int.h>
#include <structs.h>
#include <2d-R3.h>

using namespace std;

//Run parameters which you might want to change
const int N0=31;
const int N1=12;
const int mIT=10;
const double Mmax=6; //maximal delta m2 in root finder
const int ITmax=40;
const double interaction=1.0;
const double cut=atan(1000); //maximum p2 where to cut off integrals
//

//don't change anything below unless
//you know what you're doing


int ADD=2;
char T;


int ABORT=0;


const double EulerGamma=0.5772156649015329;


double s2p0vals[N0],c2p0vals[N0],tp0vals[N0];
double c2p1vals[N0],tp1vals[N0];
double p0vals[N0];
double p0weights[N0];
double modp0weights[N0];
double Smodp0weights[N0];

double phivals[2*N1],phiweights[2*N1];

vertex Gamma(N0,N0,2*N1); //vertex
vertex nGamma(N0,N0,2*N1); //new vertex
vertex SGamma(N0,N0,2*N1); //subtracted vertex
vertex SnGamma(N0,N0,2*N1); //new subtracted vertex
double Gamma1[N0];

double deltaPi[N0],ndeltaPi[N0];
double Sigma[N0];


double help[N0];

double oldSigma[N0],oldPi[N0];




#include <rootfinding.cpp>


void load_quadrature_data()
{
  

  char fname[255];
  sprintf(fname,"LegendreData/LegQuad-N%i.dat",N0-1);
  
  printf("===> Info: Loading Quadrature Data from File %s\n",fname);

  fstream parameters;
  parameters.open(fname, ios::in);
  if(parameters.is_open())
    {
      
      for (int i=0;i<N0;i++)
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
  //double cut=atan(500);
  for (int i=0;i<N0;i++)
    {
      s2p0vals[i]=sin(p0vals[i]*M_PI*0.5)*sin(p0vals[i]*M_PI*0.5);
      c2p0vals[i]=cos(p0vals[i]*M_PI*0.5)*cos(p0vals[i]*M_PI*0.5);
      tp0vals[i]=tan(p0vals[i]*M_PI*0.5);

      c2p1vals[i]=cos(p0vals[i]*cut)*cos(p0vals[i]*cut);
      tp1vals[i]=tan(p0vals[i]*cut);
      modp0weights[i]=p0weights[i];
      
    }
  modp0weights[0]=0.5*p0weights[0];

  for (int i=0;i<N0;i++)
    Smodp0weights[i]=modp0weights[i]/c2p0vals[i];

  printf("===> Info: Using %i points for angular quadrature integration\n",2*N1);
  
  for (int i=0;i<2*N1;i++)
    {
      phivals[i]=i*M_PI/N1;
      phiweights[i]=M_PI/N1;
    }
}

void exportfield(double *field,const char *name,double *xvals)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=0;i<N0;i++)
    out << atan(xvals[i])*2/M_PI << "\t" << field[i] << "\n";
  
  out.close();
}



void measureexp(double *field,const char *name,double *xvals,double shiftm)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=1;i<N0-1;i++)
    {
      out << xvals[i] << "\t" << log((xvals[i]+1+shiftm-field[i])/(xvals[i-1]+1+shiftm-field[i-1]))/(log(xvals[i]/xvals[i-1]));
      out << "\t" << log((xvals[i]+1+shiftm-field[i])/(xvals[i-1]+1+shiftm-field[i-1]))/(xvals[i]-xvals[i-1])*(xvals[i]+xvals[i-1])*0.5;
      out << "\t" << xvals[i]+1+shiftm-field[i] << "\n";
    }
  out.close();
}

double calcc1(double *field,const char *name,double *xvals,double shiftm,double p2)
{

  int a=0;
  while (p2>xvals[a])
    a++;
  double yval=log((xvals[a]+1+shiftm-field[a])/(xvals[a-1]+1+shiftm-field[a-1]));
  double xv=log(xvals[a]/xvals[a-1]);

  return yval/xv;
}

void Nexportfield(double *field,const char *name,double *xvals,int len)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=0;i<len;i++)
    out << xvals[i] << "\t" << field[i] << "\n";
  
  out.close();
}

double int_field(double req,double *field,double *xvals)
{
  //printf("if %f\n",req);
  int a=0;
  double out=0;
  if (req==xvals[0])
    {
      //printf("hello %f\n",field[0]);
      out=field[0];
    }
  else if (req<xvals[N0-1])
    {
      while(req>xvals[a])
	a++;

      out=linear_int(xvals[a-1],xvals[a],field[a-1],field[a],req);
    }
  else
    {
      //printf("req=%f out of bounds returning 0 max %f\n",req,xvals[N0-1]);
      out=0;
    }
    //out=field[N0-1];
  //printf("fi tp0[%i]=%f<%f<tp0[%i]=%f\n",a-1,tp0vals[a-1],req,a,tp0vals[a]);
  
  return out;
}



double highS(double p2,double shiftm)
{
  double res=0;
  double m2=1+shiftm;
  double ag=sqrt(p2/(4*m2+p2));
  res=4*atanh(ag)*ag/p2;
  res*=2/(4*M_PI);
  return res;
}

double newbigS(double p2,double shiftm)
{
  double res=0;
  double s1,s2;
  double max=tp1vals[N0-1];
  for (int i=0;i<N0;i++)
    for (int j=0;j<N0;j++)
      {
        double sval1=(tp0vals[i]-sqrt(p2))*(tp0vals[i]-sqrt(p2))+tp0vals[j]*tp0vals[j];
        double sval2=(tp0vals[i]+sqrt(p2))*(tp0vals[i]+sqrt(p2))+tp0vals[j]*tp0vals[j];
	s1=int_field(sval1,deltaPi,tp1vals);
	s2=int_field(sval2,deltaPi,tp1vals);
	double reg=tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j];
	double pp=int_field(reg,deltaPi,tp1vals);
       
        double tprop1=1/(sval1+1+shiftm-s1);
	double tprop2=1/(sval2+1+shiftm-s2);

	double ang1=acos((sqrt(p2)*tp0vals[i]-reg)/sqrt(reg*(p2+reg-2*sqrt(p2)*tp0vals[i])));
	double ga1=Gamma.fullinterpolate(reg,sval1,ang1,tp0vals,phivals,N0,2*N1,0,1);

	double ang2=acos((-sqrt(p2)*tp0vals[i]-reg)/sqrt(reg*(p2+reg+2*sqrt(p2)*tp0vals[i])));
	double ga2=Gamma.fullinterpolate(reg,sval2,ang2,tp0vals,phivals,N0,2*N1,0,1);
	
	res+=modp0weights[i]*modp0weights[j]/c2p0vals[i]/c2p0vals[j]/(reg+1+shiftm-pp)*(tprop1*ga1+tprop2*ga2);
	//res+=modp0weights[i]*modp0weights[j]/c2p0vals[i]/c2p0vals[j]/(reg+1+shiftm-pp)*(tprop1+tprop2);
      }
  res/=8.;
  res*=2;
  return res;
}


double bigSprime(double p2,double shiftm)
{
  double res=0;
  double s1,s2;
  double max=tp1vals[N0-1];
  for (int i=0;i<N0;i++)
    for (int j=0;j<N0;j++)
      {
	double sval1=(tp0vals[i]-sqrt(p2))*(tp0vals[i]-sqrt(p2))+tp0vals[j]*tp0vals[j];
	double sval2=(tp0vals[i]+sqrt(p2))*(tp0vals[i]+sqrt(p2))+tp0vals[j]*tp0vals[j];
	s1=int_field(sval1,deltaPi,tp1vals);
	s2=int_field(sval2,deltaPi,tp1vals);
	double reg=tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j];
	double pp=int_field(reg,deltaPi,tp1vals);

	double tprop1=1/(sval1+1+shiftm-s1);
	double tprop2=1/(sval2+1+shiftm-s2);

	double ang1=acos((sqrt(p2)*tp0vals[i]-reg)/sqrt(reg*(p2+reg-2*sqrt(p2)*tp0vals[i])));
	double ga1=Gamma.fullinterpolate(reg,sval1,ang1,tp0vals,phivals,N0,2*N1,0,1);

	double ang2=acos((-sqrt(p2)*tp0vals[i]-reg)/sqrt(reg*(p2+reg+2*sqrt(p2)*tp0vals[i])));
	double ga2=Gamma.fullinterpolate(reg,sval2,ang2,tp0vals,phivals,N0,2*N1,0,1);

	res+=modp0weights[i]*modp0weights[j]/c2p0vals[i]/c2p0vals[j]/(reg+1+shiftm-pp)*(tprop1*ga1*(ga1-1)+tprop2*ga2*(ga2-1));

      }
  res/=8.;
  res*=2;
  return res;
}




 
double newgetdeltaPi(double p2,double shiftm,double lambda)
{
  double res=0;
  double s1,s2;
  double max=tp1vals[N0-1];
  
  for (int i=0;i<N0;i++)
    for (int j=0;j<N0;j++)
      {
        double sval1=(tp0vals[i]-sqrt(p2))*(tp0vals[i]-sqrt(p2))+tp0vals[j]*tp0vals[j];
        double sval2=(tp0vals[i]+sqrt(p2))*(tp0vals[i]+sqrt(p2))+tp0vals[j]*tp0vals[j];
	if (sval1<max)
	  s1=int_field(sval1,Sigma,tp1vals);
	else
	  s1=highS(sval1,shiftm);
	if (sval2<max)
          s2=int_field(sval2,Sigma,tp1vals);
        else
          s2=highS(sval2,shiftm);


	double tprop1=2*lambda*s1/(1+2*lambda*s1);
	double tprop2=2*lambda*s2/(1+2*lambda*s2);

	double reg=tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j];
	double pp=int_field(reg,deltaPi,tp1vals);
        

	double ang1=atan2(-tp0vals[j],-tp0vals[i]);	
	double ang2=atan2(-tp0vals[j],tp0vals[i]);
	double ga1=Gamma.fullinterpolate(p2,reg,ang1,tp0vals,phivals,N0,2*N1,0,1);
	double ga2=Gamma.fullinterpolate(p2,reg,ang2,tp0vals,phivals,N0,2*N1,0,1);
	
	
	res+=modp0weights[i]*modp0weights[j]/c2p0vals[i]/c2p0vals[j]/(reg+1+shiftm-pp)*(tprop1*ga1+tprop2*ga2);

	
	res-=modp0weights[i]*modp0weights[j]/c2p0vals[i]/c2p0vals[j]/(reg+1+shiftm-pp)*(ga1+ga2-2);
       	
	
      }
  res*=lambda;

  //hm.free();
  return res;
}


double get_shiftm(double shiftm, double lambda)
{
  
  double num=0;
  double den=0;
  for (int i=0;i<N0;i++)
    {
      double pp=int_field(tp0vals[i],deltaPi,tp1vals);
      //printf("x=%f pp=%f\n",tp0vals[i],pp);
      num+=modp0weights[i]/c2p0vals[i]*pp/(tp0vals[i]+1+shiftm-pp)/(tp0vals[i]+1);
      den+=modp0weights[i]/c2p0vals[i]/(tp0vals[i]+1+shiftm-pp)/(tp0vals[i]+1);

    }
  num/=8.;
  den/=8.;

  double res=12*lambda*num/(1+12*lambda*den);
  return res;
}

double NdGamma(struct coord *hm,double shiftm,double lambda)
{
  double res=0;
  double max=tp1vals[N0-1];
  double p0=sqrt(tp0vals[hm->x[0]]);
  double k0=cos(phivals[hm->x[2]])*sqrt(tp0vals[hm->x[1]]);
  double k1=sin(phivals[hm->x[2]])*sqrt(tp0vals[hm->x[1]]);

  

  //#pragma omp parallel for reduction( + : res)
    for (int i=0;i<N0;i++)
      for (int j=0;j<N0;j++)
	{
	  double s1,s2,t1,t2,pi1,pi2,g1,g2;
	  double ga1=1,ga2=1,ga3=1;
	double reg1a=(tp0vals[i]+k0)*(tp0vals[i]+k0)+(tp0vals[j]+k1)*(tp0vals[j]+k1);
	double reg1b=(tp0vals[i]+k0)*(tp0vals[i]+k0)+(-tp0vals[j]+k1)*(-tp0vals[j]+k1);
	double reg1c=(tp0vals[i]-k0)*(tp0vals[i]-k0)+(tp0vals[j]+k1)*(tp0vals[j]+k1);
	double reg1d=(tp0vals[i]-k0)*(tp0vals[i]-k0)+(-tp0vals[j]+k1)*(-tp0vals[j]+k1);
        double arg1a=(tp0vals[i]+p0+k0)*(tp0vals[i]+p0+k0)+(tp0vals[j]+k1)*(tp0vals[j]+k1);
	double arg1b=(tp0vals[i]+p0+k0)*(tp0vals[i]+p0+k0)+(tp0vals[j]-k1)*(tp0vals[j]-k1);
	double arg1c=(tp0vals[i]-p0-k0)*(tp0vals[i]-p0-k0)+(tp0vals[j]+k1)*(tp0vals[j]+k1);
	double arg1d=(tp0vals[i]-p0-k0)*(tp0vals[i]-p0-k0)+(tp0vals[j]-k1)*(tp0vals[j]-k1);
	double arg2a=(tp0vals[i])*(tp0vals[i])+(tp0vals[j])*(tp0vals[j]);

	if (reg1a<max)
	  s1=int_field(reg1a,Sigma,tp1vals);
	else
	  s1=highS(reg1a,shiftm);
	
	t1=2*lambda/(1+2*lambda*s1);

	pi1=int_field(arg1a,deltaPi,tp1vals);
	pi2=int_field(arg2a,deltaPi,tp1vals);	
	g1=1/(arg1a+1+shiftm-pi1);
	g2=1/(arg2a+1+shiftm-pi2);

	double ang1=-atan2(tp0vals[j],tp0vals[i])+atan2(k1,k0);
	ga1=Gamma.fullinterpolate(arg2a,k0*k0+k1*k1,ang1,tp0vals,phivals,N0,2*N1,0,1);
	double ang2=atan2(-k1-tp0vals[j],-p0-k0-tp0vals[i]);
	ga2=Gamma.fullinterpolate(p0*p0,arg1a,ang2,tp0vals,phivals,N0,2*N1,0,1);
	double ang3=-atan2(k1+tp0vals[j],p0+k0+tp0vals[i])+atan2(-tp0vals[j],-tp0vals[i]);
	ga3=Gamma.fullinterpolate(arg1a,arg2a,ang3,tp0vals,phivals,N0,2*N1,0,1);
				  
	res+=Smodp0weights[i]*Smodp0weights[j]*t1*g1*g2*ga1*ga2*ga3;

	if (reg1b<max)
	  s1=int_field(reg1b,Sigma,tp1vals);
	else
	  s1=highS(reg1b,shiftm);
	
	t1=2*lambda/(1+2*lambda*s1);

	pi1=int_field(arg1b,deltaPi,tp1vals);
	pi2=int_field(arg2a,deltaPi,tp1vals);
	g1=1/(arg1b+1+shiftm-pi1);
	g2=1/(arg2a+1+shiftm-pi2);

	ang1=-atan2(-tp0vals[j],tp0vals[i])+atan2(k1,k0);
	ga1=Gamma.fullinterpolate(arg2a,k0*k0+k1*k1,ang1,tp0vals,phivals,N0,2*N1,0,1);
	ang2=atan2(-k1+tp0vals[j],-p0-k0-tp0vals[i]);
	ga2=Gamma.fullinterpolate(p0*p0,arg1b,ang2,tp0vals,phivals,N0,2*N1,0,1);
	ang3=-atan2(k1-tp0vals[j],p0+k0+tp0vals[i])+atan2(tp0vals[j],-tp0vals[i]);
	ga3=Gamma.fullinterpolate(arg1b,arg2a,ang3,tp0vals,phivals,N0,2*N1,0,1);

	res+=Smodp0weights[i]*Smodp0weights[j]*t1*g1*g2*ga1*ga2*ga3;

	if (reg1c<max)
	  s1=int_field(reg1c,Sigma,tp1vals);
	else
	  s1=highS(reg1c,shiftm);
	
	t1=2*lambda/(1+2*lambda*s1);
	
	pi1=int_field(arg1c,deltaPi,tp1vals);
	pi2=int_field(arg2a,deltaPi,tp1vals);	
	g1=1/(arg1c+1+shiftm-pi1);
	g2=1/(arg2a+1+shiftm-pi2);

	ang1=-atan2(tp0vals[j],-tp0vals[i])+atan2(k1,k0);
	ga1=Gamma.fullinterpolate(arg2a,k0*k0+k1*k1,ang1,tp0vals,phivals,N0,2*N1,0,1);
	ang2=atan2(-k1-tp0vals[j],-p0-k0+tp0vals[i]);
	ga2=Gamma.fullinterpolate(p0*p0,arg1c,ang2,tp0vals,phivals,N0,2*N1,0,1);
	ang3=-atan2(k1+tp0vals[j],p0+k0-tp0vals[i])+atan2(-tp0vals[j],tp0vals[i]);
	ga3=Gamma.fullinterpolate(arg1c,arg2a,ang3,tp0vals,phivals,N0,2*N1,0,1);
	
	res+=Smodp0weights[i]*Smodp0weights[j]*t1*g1*g2*ga1*ga2*ga3;

	
	if (reg1d<max)
	  s1=int_field(reg1d,Sigma,tp1vals);
	else
	  s1=highS(reg1d,shiftm);
	
	t1=2*lambda/(1+2*lambda*s1);
	
	
	pi1=int_field(arg1d,deltaPi,tp1vals);	
	pi2=int_field(arg2a,deltaPi,tp1vals);
	g1=1/(arg1d+1+shiftm-pi1);
	g2=1/(arg2a+1+shiftm-pi2);

	ang1=-atan2(-tp0vals[j],-tp0vals[i])+atan2(k1,k0);
	ga1=Gamma.fullinterpolate(arg2a,k0*k0+k1*k1,ang1,tp0vals,phivals,N0,2*N1,0,1);
	ang2=atan2(-k1+tp0vals[j],-p0-k0+tp0vals[i]);
	ga2=Gamma.fullinterpolate(p0*p0,arg1d,ang2,tp0vals,phivals,N0,2*N1,0,1);
	ang3=-atan2(k1-tp0vals[j],p0+k0-tp0vals[i])+atan2(tp0vals[j],tp0vals[i]);
	ga3=Gamma.fullinterpolate(arg1d,arg2a,ang3,tp0vals,phivals,N0,2*N1,0,1);
	
	res+=Smodp0weights[i]*Smodp0weights[j]*t1*g1*g2*ga1*ga2*ga3;
      }
  
  res/=16.;
  res*=4;
  
  
  
  return res;
  
}


double NdGamma(int p,int k, int f, double shiftm,double lambda)
{

  struct coord hm;
  hm.d=3;
  hm.allocate();

  hm.x[0]=p;
  hm.x[1]=k;
  hm.x[2]=f;
  double val=NdGamma(&hm,shiftm,lambda);		
  hm.free();
  return val;
}

  


double shiftfind(double *args)
{
  double sm=args[0];
  double lambda=args[5];
  return get_shiftm(sm,lambda)-sm;
}


double altget_allshiftm(double shiftm,double lambda,int verbose)
{
  
  double Rag[6]={1,ITmax,max(deltaPi[0]-1+0.001,0),Mmax,double(verbose),lambda};
  double newshiftm=convergeroot2(&shiftfind,Rag);
  return newshiftm;
}

void getall_sigma(double shiftm)
{
#pragma omp parallel for 
  for (int i=0;i<N0;i++)
    {
      Sigma[i]=newbigS(tp1vals[i],shiftm);
      //printf("%i done\n",i);
    }
}

void getall_Gamma1(double shiftm,double lambda)
{
  for (int i=0;i<N0;i++)
    Gamma1[i]=1/(1+4*lambda*Sigma[i])-1;	
}

void getall_Pi(double shiftm,double lambda)
{
#pragma omp parallel for 
  for (int i=0;i<N0;i++)
    {
      ndeltaPi[i]=newgetdeltaPi(tp1vals[i],shiftm,lambda);
      //printf("%i done %f\n",i,deltaPi[i]);                                                                                   
    }
}



void Pgetall_gamma(double shiftm,double lambda,int verbose)
{

  int steppi=N0/10;

  if (verbose)
    {
       printf("\t(");
      
      
      for (int a=0;a<N0/steppi;a++)
	printf("%i",a);
      
      printf(" )\n\t ");
    }
  
#pragma omp parallel for 
  for (int i=0;i<N0;i++)
    {
      
      struct coord hm;
      hm.d=3;
      hm.allocate();
      hm.x[0]=i;
      for (hm.x[1]=0;hm.x[1]<=i;hm.x[1]++)
	for (hm.x[2]=0;hm.x[2]<2*N1;hm.x[2]++)
	  {

	    double val=NdGamma(&hm,shiftm,lambda);
	    nGamma.set(hm,1-val);

	  }
      if (verbose)
	{
	  if (i%steppi==0)
	    {	  
	      printf("*");
	      fflush(stdout);
	    }
	}
      hm.free();
    }

  

  //completion of data type
  struct coord hm;
  hm.d=3;
  hm.allocate();
  for (hm.x[0]=0;hm.x[0]<N0;hm.x[0]++)
    for (hm.x[1]=hm.x[0]+1;hm.x[1]<N0;hm.x[1]++)
      for (hm.x[2]=0;hm.x[2]<2*N1;hm.x[2]++)
	{	  
	  double val=nGamma.get(hm.x[1],hm.x[0],hm.x[2]);
	  nGamma.set(hm,val);
	}

  hm.free();


}


#include<mf-2d-R3.cpp>

int main(int argc, char* argv[])
{
  load_quadrature_data();
  derived_data();

  char buffer[255];
  double lambda=interaction;
  double shiftm=0.0;
  double oldshiftm=shiftm;

  T='E';
  
  if(ADD==2)
    T='S';
  if(ADD==3)
    T='T';
  if(ADD==4)
    T='Q';

  printf("===> Info: Maximum deltaM^2 expected in Root finder: %f\n",Mmax);

  fstream hist;
  sprintf(buffer,"store/history-R3zeroTfinal%cN%i-A%i-cut%i-lam%f-IT%i.dat",T,N0,N1,int(tan(cut)),lambda,mIT);
  //printf("hi %s\n",buffer);
  hist.open(buffer,ios::out);

  hist << "#Iterations \t";
  hist << "delta m2 \t";
  hist << "delta Pi[0] \t";
  hist << "Lambda \t";
  hist << "M/mR \t";
  hist << "c2\t";
  hist << "delta c2\t";
  hist << endl;
  hist.flush();
  
  Gamma.set_unity();
  //SGamma.set_zero();
  

  for (int i=0;i<N0;i++)
    {
      Sigma[i]=0;
      deltaPi[i]=0;
      ndeltaPi[i]=0;
    }

  double cs,cp,cv,ex,mest;
  double c1,c2;

  
  for (int IT=1;IT<mIT;IT++)
    {

      
      printf("IT=%i\n",IT);
      getall_sigma(shiftm);
      printf("\t sigma done %f\n",int_field(1.0,Sigma,tp1vals));
      
      for (int hi=0;hi<1;hi++)
	{
	  
	  Pgetall_gamma(shiftm,lambda,1);

	  if (ADD==2)
	    Gamma.addhalf(&nGamma);
	  else if (ADD==3)
	    Gamma.addthird(&nGamma);
	  else if(ADD==4)
	    Gamma.addquarter(&nGamma);
	  else
	    Gamma.setequal(&nGamma);

	  
	}
     
      
      getall_Pi(shiftm,lambda);

      
      sprintf(buffer,"R3-Pidump%i",IT);                                                                       
      exportfield(ndeltaPi,buffer,tp1vals);

      
      
      printf("\t pi done %f\n",int_field(1.0,ndeltaPi,tp1vals));
      for (int i=0;i<N0;i++)
	deltaPi[i]=ndeltaPi[i];

      
      shiftm=altget_allshiftm(shiftm,lambda,0);
      printf("\t shiftm=%f\n",shiftm);
      
      
      cs=0.5*R3measureCS(shiftm,lambda);
      cp=0.5*R3measureCP(shiftm,lambda);
      cv=R3measureVL(shiftm,lambda);
      ex=-0.25*R3measureextra(shiftm,lambda);
      mest=calcmest(shiftm);
      c1=1-calcc1(deltaPi,buffer,tp1vals,shiftm,1);
      c2=1-calcc1(deltaPi,buffer,tp1vals,shiftm,1.5);
      
      printf("cs=%f cp=%f vl=%f ex=%f sum=%f mest=%f c1=%f c2=%f\n",cs,cp,cv,ex,cs+cp+cv+ex,sqrt(mest),c1,c2);

      
      hist << IT << "\t";
      hist << shiftm << "\t";
      hist << deltaPi[0] << "\t";
      hist << cs+cp+cv+ex << "\t";
      hist << sqrt(mest) << "\t";
      hist << 0.5*(c1+c2) << "\t";
      hist << fabs(c1-c2) << "\t";
      hist << endl;
      hist.flush();
      
      
    }

  hist.close();
  
  
  fstream final;
  sprintf(buffer,"store/R3zeroTfinalN%i-A%i-cut%i-lam%f-IT%i.dat",N0,N1,int(tan(cut)),lambda,mIT);
  
  final.open(buffer,ios::out);
  final << "#lambda" << "\t";
  final << "shiftm2" << "\t\t";
  final << "VacEnergy" << "\t";
  final << "M(est)" << "\t";
  final << "c1" << "\t";
  final << "dc1" << "\t";
  final << endl;
  
  final << lambda << "\t";
  final << shiftm << "\t";
  final << cs+cp+cv+ex << "\t";
  final << sqrt(mest) << "\t";
  final << 0.5*(c1+c2) << "\t";
  final << fabs(c1-c2) << "\t";
  final << endl;
  final.flush();
  final.close();

}
