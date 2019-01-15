double R3measureCS(double shiftm, double lambda)
{
  double res=0;
  double ss;
  double max=tp1vals[N0-1];
//#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    {
      if (tp0vals[i]<=max)
	ss=newbigS(tp0vals[i],shiftm);
      else
	ss=highS(tp0vals[i],shiftm);
      res=res + modp0weights[i]/c2p0vals[i]*(log(1+2*lambda*ss)-2*lambda*ss/(1+2*lambda*ss));
    }
  res/=8.;
  return res;
}


double R3measureextra(double shiftm, double lambda)
{
  double res=0;
  double ss,pp;
  double max=tp1vals[N0-1];
  //#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    {
      if (tp0vals[i]<=max)
	{
	  pp=bigSprime(tp0vals[i],shiftm);
	  ss=newbigS(tp0vals[i],shiftm);
	}
      else
	{
	  ss=highS(tp0vals[i],shiftm);
	  pp=bigSprime(tp0vals[i],shiftm);
	}
      res=res + modp0weights[i]/c2p0vals[i]*2*lambda/(1+2*lambda*ss)*pp;
    }
  res/=8.;
  return res;
}



double R3measureCP(double shiftm,double lambda)
{
  double res=0;
  double ss;
  double max=tp1vals[N0-1];
  //#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    {
      if (tp0vals[i]<=max)
	ss=int_field(tp0vals[i],deltaPi,tp1vals);
      //newgetdeltaPi(tp0vals[i]+0.000001,shiftm,lambda);
      else
	ss=0;     
      res=res+modp0weights[i]/c2p0vals[i]*(log((tp0vals[i]+1+shiftm-ss)/(tp0vals[i]+1+shiftm))+0.5*ss/(tp0vals[i]+1+shiftm-ss));
      //printf("i=%i res=%f 1+shiftm=%f ss=%f \n",i,res,1+shiftm,ss);
    }
  res/=8.;
  return res;
}


double R3measureVL(double shiftm,double lambda)
{
  double vr=-shiftm*shiftm/(48*lambda);

  double logv=(shiftm+1)*(1+log(1/(shiftm+1)))-1;
  logv/=(8*M_PI);
  //printf("vr=%f logv=%f\n",vr,logv);
  return (vr+logv);


}

double calcmest(double shiftm)
{
  double num=1+shiftm-deltaPi[0];
  double den=1-(deltaPi[1]-deltaPi[0])/(tp1vals[1]-tp1vals[0]);

  //printf("num=%f den=%f\n",num,den);
  return num/den;
}
