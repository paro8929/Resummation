#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <complex>
#include <cmath>
#include <linear_int.h>


using namespace std;


//Run parameters which you might want to change
const int N0=101;
const int mIT=50;
const double Mmax=60; //maximal delta m2 in root finder
const int ITmax=40; //maximum number of iterations in root finder
const double interaction=1.0;
const double cut=atan(2000); //maximum p2 where to cut off integrals
//

//don't change anything below unless
//you know what you're doing


int ABORT=0;

const double EulerGamma=0.5772156649015329;


double s2p0vals[N0],c2p0vals[N0],tp0vals[N0];
double c2p1vals[N0],tp1vals[N0];
double p0vals[N0];
double p0weights[N0];
double modp0weights[N0];

//double phivals[2*N1],phiweights[2*N1];

double deltaPi[N0],ndeltaPi[N0];
double Sigma[N0];




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

  //printf("===> Info: Using %i points for angular quadrature integration\n",2*N1);
  /*
  for (int i=0;i<2*N1;i++)
    {
      phivals[i]=i*M_PI/N1;
      phiweights[i]=M_PI/N1;
      }*/
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

/*
double testme(double nn)
{
  double res=0;
  for (int j=0;j<2*N1;j++)
    {
      res+=phiweights[j]*cos(nn*phivals[j]);
    }
  return 0.5*res;
  }*/


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
  //  double s1,s2;
  double max=tp1vals[N0-1];
#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N0;j++)
      {
        double sval1=(tp0vals[i]-sqrt(p2))*(tp0vals[i]-sqrt(p2))+tp0vals[j]*tp0vals[j];
        double sval2=(tp0vals[i]+sqrt(p2))*(tp0vals[i]+sqrt(p2))+tp0vals[j]*tp0vals[j];
	double s1=int_field(sval1,deltaPi,tp1vals);
	double s2=int_field(sval2,deltaPi,tp1vals);
	double pp=int_field(tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j],deltaPi,tp1vals);
       
        double tprop1=1/(sval1+1+shiftm-s1);
	double tprop2=1/(sval2+1+shiftm-s2);
	res+=modp0weights[i]*modp0weights[j]/c2p0vals[i]/c2p0vals[j]/(tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j]+1+shiftm-pp)*(tprop1+tprop2);
      }
  res/=8.;
  res*=2;
  return res;
}





 
double newgetdeltaPi(double p2,double shiftm,double lambda)
{
  double res=0;
  //  double s1,s2;
  double max=tp1vals[N0-1];
  //printf("max=%f\n",max);
#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    for (int j=0;j<N0;j++)
      {
        double sval1=(tp0vals[i]-sqrt(p2))*(tp0vals[i]-sqrt(p2))+tp0vals[j]*tp0vals[j];
        double sval2=(tp0vals[i]+sqrt(p2))*(tp0vals[i]+sqrt(p2))+tp0vals[j]*tp0vals[j];
	double s1,s2;
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
	double pp=int_field(tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j],deltaPi,tp1vals);
        //printf("%f adding %f to %f\n",tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j],1/(tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j]+1+shiftm)*(tprop1+tprop2),res);
        //if (fabs(pp)>1.e-6)
	//  printf("i=%i j%i pp=%f\n",i, j,pp);
	//pp=0;
	res+=modp0weights[i]*modp0weights[j]/c2p0vals[i]/c2p0vals[j]/(tp0vals[i]*tp0vals[i]+tp0vals[j]*tp0vals[j]+1+shiftm-pp)*(tprop1+tprop2);
      }
  res*=lambda;
  return res;
}

double get_shiftm(double shiftm, double lambda)
{
  
  double num=0;
  double den=0;
#pragma omp parallel for reduction( + : num,den)
  for (int i=0;i<N0;i++)
    {
      double pp=int_field(tp0vals[i],deltaPi,tp1vals);
      //printf("x=%f pp=%f\n",tp0vals[i],pp);
      num+=modp0weights[i]/c2p0vals[i]*pp/(tp0vals[i]+1+shiftm-pp)/(tp0vals[i]+1);
      den+=modp0weights[i]/c2p0vals[i]/(tp0vals[i]+1+shiftm-pp)/(tp0vals[i]+1);
      //den+=modp0weights[i]/c2p0vals[i]/(tp0vals[i]+1)/(tp0vals[i]+1);
    }
  num/=8.;
  den/=8.;
  //printf("num=%f den=%f\n",num,den);
  double res=12*lambda*num/(1+12*lambda*den);
  return res;
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
  for (int i=0;i<N0;i++)
    {
      Sigma[i]=newbigS(tp1vals[i],shiftm);
      //printf("%i done\n",i);
    }
}

void getall_Pi(double shiftm,double lambda)
{
  for (int i=0;i<N0;i++)
    {
      ndeltaPi[i]=newgetdeltaPi(tp1vals[i],shiftm,lambda);
      //printf("%i done %f\n",i,deltaPi[i]);                                                                                   
    }
}



#include<mf-2d-R2.cpp>

int main(int argc, char* argv[])
{
  load_quadrature_data();
  derived_data();


  char buffer[255];
  double lambda=interaction;
  double shiftm=0.0;

  fstream hist;
  sprintf(buffer,"store/history-R2zeroT-N%i-cut%i-lam%f-IT%i.dat",N0,int(tan(cut)),lambda,mIT);
  hist.open(buffer,ios::out);
  
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
      getall_Pi(shiftm,lambda);
      //printf("test %f\n",newgetdeltaPi(1,shiftm,lambda));
      sprintf(buffer,"R2-Pidump%i",IT);                                                                       
      exportfield(ndeltaPi,buffer,tp1vals);  
      printf("\t pi done %f\n",int_field(1.0,ndeltaPi,tp1vals));
      for (int i=0;i<N0;i++)
	deltaPi[i]=ndeltaPi[i];
      shiftm=altget_allshiftm(shiftm,lambda,0);
      printf("\t shiftm=%f\n",shiftm);

      cs=0.5*altmeasureCS(shiftm,lambda);
      cp=0.5*altmeasureCP(shiftm,lambda);
      cv=altmeasureVL(shiftm,lambda);
      mest=calcmest(shiftm);
      c1=1-calcc1(deltaPi,buffer,tp1vals,shiftm,1);
      c2=1-calcc1(deltaPi,buffer,tp1vals,shiftm,1.5);


      
      hist << IT << "\t";
      hist << shiftm << "\t";
      hist << deltaPi[0] << "\t";
      hist << cs+cp+cv << "\t";
      hist << sqrt(mest) << "\t";
      hist << 0.5*(c1+c2) << "\t";
      hist << fabs(c1-c2) << "\t";
      hist << endl;
      hist.flush();
      
      printf("cs=%f cp=%f vl=%f  sum=%f mest=%f c1=%f c2=%f\n",cs,cp,cv,cs+cp+cv,sqrt(mest),c1,c2);
    }
  hist.close();
  

  printf("cs=%f cp=%f vl=%f\n",cs,cp,cv);
  printf("sum=%g\n",cs+cp+cv);

  fstream final;
  sprintf(buffer,"store/R2zeroTfinalN%i-cut%i-lam%f-IT%i.dat",N0,int(tan(cut)),lambda,mIT);
  
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
  final << cs+cp+cv << "\t";
  final << sqrt(mest) << "\t";
  final << 0.5*(c1+c2) << "\t";
  final << fabs(c1-c2) << "\t";
  final << endl;
  final.flush();
  final.close();

  //printf("Sigma(p2=1) %f vs high S=%f newS=%f\n",bigS(1,0),highS(1,0),newbigS(1,0));
  //getall_sigma(0);
  //printf("sigma done %f\n",intfield(1.0,Sigma,tp1vals));
  //printf("testing %f\n",int_field(1.0,Sigma,tp1vals));
  //printf("Tprop %f\n",propT(1,0,lambda));
  //printf("Tprop %f\n",propT(1000,0,lambda));
  //printf("deltaPi(p2=0) %f new=%f \n",getdeltaPi(0,0,lambda),newgetdeltaPi(0,0,lambda));
  //printf("deltaPi(p2=500) %f new=%f\n",getdeltaPi(500,0,lambda),newgetdeltaPi(500,0,lambda));
  //printf("deltaPi(p2=10000) %f new=%f\n",getdeltaPi(10000,0,lambda),newgetdeltaPi(10000,0,lambda));
  //getall_Pi(0,lambda);
  //sprintf(buffer,"Pidump%i",1);
  //exportfield(deltaPi,buffer,tp1vals);
  //sprintf(buffer,"sigmadump%i",1);
  //exportfield(Sigma,buffer,tp1vals);


  //double Gm=altget_allshiftm(0,lambda,1);
  //printf("shiftm0=%f vs find=%f\n",get_shiftm(shiftm,lambda),Gm);
 
  //double Gm=altget_allshiftm(0,lambda,0);
  //printf("shiftm0=%f vs %f\n",get_shiftm(shiftm,lambda),Gm);

  //getall_sigma(Gm);
  //printf("sigma done %f\n",Sigma[0]);
  //getall_Pi(Gm,lambda);
  //sprintf(buffer,"Pidump%i",2);
  //exportfield(deltaPi,buffer,tp1vals);
  //sprintf(buffer,"sigmadump%i",2);
  //exportfield(Sigma,buffer,tp1vals);
  /*
  fstream final;
  sprintf(buffer,"store/R2zeroTfinalN%i-A%i-cut%i-lam%f-IT%i.dat",N0,N1,int(tan(cut)),lambda,mIT);
  
  final.open(buffer,ios::out);
  final << "#lambda" << "\t";
  final << "shiftm2" << "\t\t";
  final << "VacEnergy" << "\t";
  final << endl;

  final << lambda << "\t";
  final << shiftm << "\t";
  final << cs+cp+cv << "\t";
  final << endl;
  final.flush();
  final.close();*/
}
