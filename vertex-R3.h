class vertex : public data_array
{
 private:
  struct coord help;
 public:
  vertex(int k0,int p0,int k1,int p1)
    {
      help.d=4;
      help.allocate();
      help.x[0]=k0;
      help.x[1]=p0;
      help.x[2]=k1;
      help.x[3]=p1;
      init(help);
      alloc_data();

    }
  void set_zero();
  void set_unity();
  void set(struct coord,double);
  double get(struct coord);
  void set(int,int,int,int,double);
  double get(int,int,int,int);
  double simple_interpolate4(double,double,double,double,double *,double*);
  //double interpolate(int,double,double,double *,double *,int,int,int,double);
  //double fullinterpolate(double,double,double,double *,double *,int,int,int,double);
  double ipt1(double,int,int,int,double *,double*);
  double ipt2(double,double,int,int,double *,double*);
  double ipt3(double,double,double,int,double *,double*);
  double ipt4(double,double,double,double,double *,double*);
  void setequal(vertex *);
  void addtx(vertex *,double alpha);
  };


void vertex::setequal(vertex *in)
{
  for (int i=0;i<SIZE;i++)
    data[i]=in->data[i];
}

void vertex::addtx(vertex *in,double alpha)
{
  for (int i=0;i<SIZE;i++)
    {
      data[i]=(1-alpha)*data[i]+alpha*in->data[i];
    }
}


void vertex::set_zero()
{
  if (dat_alloc)
    {
      for (int i=0;i<SIZE;i++)
	data[i]=0;
    }
  else
    printf("vertex ERROR: data not initialized\n");
}

void vertex::set_unity()
{
  if (dat_alloc)
    {
      for (int i=0;i<SIZE;i++)
	data[i]=1.0;
    }
  else
    printf("vertex ERROR: data not initialized\n");
}


void vertex::set(struct coord hm,double val)
{
  if (dat_alloc)
    {
      data[getind(&hm)]=val;
      //do somethgin
    }
  else
    printf("vertex ERROR: data not initialized\n");
}


void vertex::set(int a0,int a1,int b0, int b1,double val)
{
  struct coord hm;
  hm.d=4;
  hm.allocate();
  hm.x[0]=a0;
  hm.x[1]=a1;
  hm.x[2]=b0;
  hm.x[3]=b1;
  if (dat_alloc)
    {
      data[getind(&hm)]=val;
      //do somethgin
    }
  else
    printf("vertex ERROR: data not initialized\n");

  hm.free();
}

double vertex::get(struct coord hm)
{
  double val=0;
  if (dat_alloc)
    {
      val=data[getind(&hm)];
      //do somethgin
    }
  else
    printf("vertex ERROR: data not initialized\n");
  return val;
}


double vertex::get(int a0,int a1,int b0, int b1)
{
  double val;
  struct coord hm;
  hm.d=4;
  hm.allocate();
  hm.x[0]=a0;
  hm.x[1]=a1;
  hm.x[2]=b0;
  hm.x[3]=b1;
  if (dat_alloc)
    {
      val=data[getind(&hm)];
      //do somethgin
    }
  else
    printf("vertex ERROR: data not initialized\n");

  hm.free();
  return val;
}

double vertex::ipt1(double k0,int k1,int p0,int p1,double *xvals,double*yvals)
{
  double out=1;
  //assumming k0>0
  if (k0<0)
    printf("ERROR in intk0: found k0<0\n");
  if (k0==0)
    {
      out=get(0,k1,p0,p1);
    }
  else if (k0>xvals[help.x[0]-1])
    out=1;
  else
    {
      int ek0=1;
      while (k0>xvals[ek0])
	ek0++;
      double f0=get(ek0-1,k1,p0,p1);
      double f1=get(ek0,k1,p0,p1);
      out=linear_int(xvals[ek0-1],xvals[ek0],f0,f1,k0);
    }
  return out;
}

double vertex::ipt2(double k0,double k1,int p0,int p1,double *xvals,double*yvals)
{
  double out=1;
  //assumming k1>0
  if (k1<0)
    printf("ERROR in intk0: found k0<0\n");
  if (k1==0)
    {
      out=ipt1(k0,0,p0,p1,xvals,yvals);
    }
  else if (k1>yvals[help.x[1]-1])
    out=1;
  else
    {
      int ek1=1;
      while (k1>yvals[ek1])
	ek1++;
      double f0=ipt1(k0,ek1-1,p0,p1,xvals,yvals);
      double f1=ipt1(k0,ek1,p0,p1,xvals,yvals);
      out=linear_int(yvals[ek1-1],yvals[ek1],f0,f1,k1);
    }
  return out;
}

double vertex::ipt3(double k0,double k1,double p0,int p1,double *xvals,double*yvals)
{
  double out=1;
  if (p0==0)
    {
      out=ipt2(k0,k1,help.x[0],p1,xvals,yvals);
    }
  else if (fabs(p0)>xvals[help.x[0]-1])
    out=1;
  else
    {
      if (p0>0)
	{
	  int ep0=1;
	  while (p0>xvals[ep0])
	    ep0++;
	  double f0=ipt2(k0,k1,help.x[0]+ep0-1,p1,xvals,yvals);
	  double f1=ipt2(k0,k1,help.x[0]+ep0,p1,xvals,yvals);
	  out=linear_int(xvals[ep0-1],xvals[ep0],f0,f1,p0);
	}
      else
	{
	  int ep0=1;
	  while (p0<-xvals[ep0])
	    ep0++;
	  double f0=ipt2(k0,k1,help.x[0]-ep0,p1,xvals,yvals);
	  double f1=ipt2(k0,k1,help.x[0]-ep0+1,p1,xvals,yvals);
	  out=linear_int(-xvals[ep0],-xvals[ep0-1],f0,f1,p0);
	  //printf("where f0=%f f1=%f\n",f0,f1);
	}
    }
  return out;
}

double vertex::ipt4(double k0,double k1,double p0,double p1,double *xvals,double*yvals)
{
  
  /*
  if (k0<0)
    {
      k0*=-1;
      p0*=-1;
    }
  if (k1<0)
    {
      k1*=-1;
      p1*=-1;
    }
  
  double out=1;
  if (p1==0)
    {
      out=ipt3(k0,k1,p0,help.x[1],xvals,yvals);
    }
  else if (fabs(p1)>yvals[help.x[1]-1])
    out=1;
  else
    {
      if (p1>0)
	{
	  int ep1=1;
	  while (p1>yvals[ep1])
	    ep1++;
	  double f0=ipt3(k0,k1,p0,help.x[1]+ep1-1,xvals,yvals);
	  double f1=ipt3(k0,k1,p0,help.x[1]+ep1,xvals,yvals);
	  out=linear_int(yvals[ep1-1],yvals[ep1],f0,f1,p1);
	}
      else
	{
	  int ep1=1;
	  while (p1<-yvals[ep1])
	    ep1++;
	  double f0=ipt3(k0,k1,p0,help.x[1]-ep1,xvals,yvals);
	  double f1=ipt3(k0,k1,p0,help.x[1]-ep1+1,xvals,yvals);
	  out=linear_int(-yvals[ep1],-yvals[ep1-1],f0,f1,p1);
	  //printf("where f0=%f f1=%f\n",f0,f1);
	}
    }
  */
  double out=simple_interpolate4(k0,k1,p0,p1,xvals,yvals);

  //if (fabs(out-out2)>1.e-1)
  // printf("difference sim=%f new=%f at %f %f %f %f maxv=%f,%f\n",out2,out,k0,k1,p0,p1,xvals[help.x[0]-1],yvals[help.x[1]-1]);
  return out;
}


double vertex::simple_interpolate4(double k0,double k1,double p0,double p1,double* xvals,double* yvals)
{
  //printf("V %f %f %f %f \n",k0,k1,p0,p1);
  double res=1;
  if (k0<0)
    {
      k0*=-1;
      p0*=-1;
    }
  if (k1<0)
    {
      k1*=-1;
      p1*=-1;
    }
  //int ak0=1;
  //while (k0>xvals[ak0])
  //  ak0++;

  if ((k0>xvals[help.x[0]-1])||(k1>yvals[help.x[1]-1])||(fabs(p0)>xvals[help.x[0]-1])||(fabs(p1)>yvals[help.x[1]-1]))
    res=1;
  else
    {
      int ek0=0;
      while (k0>xvals[ek0])
	ek0++;
      
      int ak1=1;
      while (k1>yvals[ak1])
	ak1++;
      
      
      int ap1=1;
      if (p1>0)
	{
	  while (p1>yvals[ap1])
	    ap1++;
	}
      else
	{
	  while (p1<-yvals[ap1])
	    ap1++;
	}
      
      //int ap0=1;
      int ep0=0;
      if (p0>0)
	{
	  //while (p0>xvals[ap0])
	  //	ap0++;
	  
	  while (p0>xvals[ep0])
	    ep0++;
	  
	}
      else
	{
	  //while (p0<-xvals[ap0])
	  //	ap0++;
	  while (p0<-xvals[-ep0])
	    ep0--;
	}
      
      //printf("interpolating %i(%i) %i %i(%i) %i\n",ak0,ek0,ak1,ap0,ep0,ap1);
      
      if (p1==0)
	{
	  res=linear_int(yvals[ak1-1],yvals[ak1],get(ek0,ak1-1,help.x[0]+ep0,help.x[1]),get(ek0,ak1,help.x[0]+ep0,help.x[1]),k1);
	}
      else if (p1>0)
	{
	  double f0=get(ek0,ak1-1,help.x[0]+ep0,help.x[1]+ap1-1);
	  double f1=get(ek0,ak1,help.x[0]+ep0,help.x[1]+ap1-1);
	  double f2=get(ek0,ak1,help.x[0]+ep0,help.x[1]+ap1);
	  double f3=get(ek0,ak1-1,help.x[0]+ep0,help.x[1]+ap1);
	  res=bilinear_int(yvals[ak1-1],yvals[ak1],yvals[ap1-1],yvals[ap1],f0,f1,f2,f3,k1,p1);
	}
      else if (p1<0)    
	{
	  double f0=get(ek0,ak1-1,help.x[0]+ep0,help.x[1]-ap1);
	  double f1=get(ek0,ak1,help.x[0]+ep0,help.x[1]-ap1);
	  double f2=get(ek0,ak1,help.x[0]+ep0,help.x[1]-ap1+1);
	  double f3=get(ek0,ak1-1,help.x[0]+ep0,help.x[1]-ap1+1);
	  res=bilinear_int(yvals[ak1-1],yvals[ak1],-yvals[ap1],-yvals[ap1-1],f0,f1,f2,f3,k1,p1);
	}
      
      //printf("\t V %f %f %f %f OK \n",k0,k1,p0,p1);
    }
  return res;
}
