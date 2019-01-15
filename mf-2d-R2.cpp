double altmeasureCS(double shiftm, double lambda)
{
  double res=0;
  //double ss;
  double max=tp1vals[N0-1];
#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    {
      double ss;
      if (tp0vals[i]<=max)
	ss=newbigS(tp0vals[i],shiftm);
      else
	ss=highS(tp0vals[i],shiftm);
      res=res + modp0weights[i]/c2p0vals[i]*(log(1+2*lambda*ss)-2*lambda*ss/(1+2*lambda*ss));
    }
  res/=8.;
  return res;
}


double altmeasureCP(double shiftm,double lambda)
{
  double res=0;
  double ss;
  double max=tp1vals[N0-1];
#pragma omp parallel for reduction( + : res)
  for (int i=0;i<N0;i++)
    {
      double ss;
      if (tp0vals[i]<=max)
	ss=int_field(tp0vals[i],deltaPi,tp1vals);
      else
	ss=0;     
      res=res+modp0weights[i]/c2p0vals[i]*(log((tp0vals[i]+1+shiftm-ss)/(tp0vals[i]+1+shiftm))+0.5*ss/(tp0vals[i]+1+shiftm-ss));
    }
  res/=8.;
  return res;
}


double altmeasureVL(double shiftm,double lambda)
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

double calcc1(double *field,const char *name,double *xvals,double shiftm,double p2)
{

  int a=0;
  while (p2>xvals[a])
    a++;
  double yval=log((xvals[a]+1+shiftm-field[a])/(xvals[a-1]+1+shiftm-field[a-1]));
  double xv=log(xvals[a]/xvals[a-1]);

  return yval/xv;
}
