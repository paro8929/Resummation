

double max(double a,double b)
{
  double rr=0;
  if (a>=b)
    rr=a;
  else
    rr=b;

  return rr;
}


class vertex : public data_array
{
 private:
  struct coord help;
 public:
  vertex(int x,int y,int z)
    {
      help.d=3;
      help.allocate();
      help.x[0]=x;
      help.x[1]=y;
      help.x[2]=z;
      init(help);
      alloc_data();

    }
  void set_zero();
  void set_unity();
  void set(struct coord,double);
  double get(struct coord);
  void set(int,int,int,double);
  double get(int,int,int);
  double interpolate(int,double,double,double *,double *,int,int,int,double);
  double fullinterpolate(double,double,double,double *,double *,int,int,int,double);
  void setequal(vertex *);
  void addhalf(vertex *);
  void addthird(vertex *);
  void addquarter(vertex *);
};


void vertex::setequal(vertex *in)
{
  for (int i=0;i<SIZE;i++)
    data[i]=in->data[i];
}

void vertex::addhalf(vertex *in)
{
  for (int i=0;i<SIZE;i++)
    {
      data[i]=0.5*data[i]+0.5*in->data[i];
    }
}

void vertex::addthird(vertex *in)
{
  for (int i=0;i<SIZE;i++)
    {
      data[i]=(2*data[i]+in->data[i])/3.;
    }
}

void vertex::addquarter(vertex *in)
{
  for (int i=0;i<SIZE;i++)
    {
      data[i]=(3*data[i]+in->data[i])/4.;
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

void vertex::set(int x,int y,int z,double val)
{
  struct coord hm;
  hm.d=3;
  hm.allocate();
  hm.x[0]=x;
  hm.x[1]=y;
  hm.x[2]=z;
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

double vertex::get(int x,int y,int z)
{
  struct coord hm;
  hm.d=3;
  hm.allocate();
  hm.x[0]=x;
  hm.x[1]=y;
  hm.x[2]=z;
  double val=0;
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

double vertex::interpolate(int p2,double xx,double ff,double *xvals,double *fvals,int m1,int m2,int verbose,double def)
{
  int a=0,b=0;
  double out=def;
  //make sure ff is positive
  while(ff<0)
    ff+=2*M_PI;
  //and between 0 and 2 Pi
  while(ff>2*M_PI)
    ff-=2*M_PI;
  
  if (xx<xvals[m1-1])
    {
      
      while(xx>xvals[a])
	a++;
      if (verbose)
	printf("\t 2d interpolate x=%i between %f %f\n",a,xvals[a-1],xvals[a]);
      if (ff<fvals[m2-1])
	{
	  while(ff>fvals[b])
	    b++;	  
	   if (verbose)
	     printf("found ang %i between %f %f\n",b,fvals[b-1],fvals[b]);
	   if (b>0)
	     {
	       if (a>0)
		 {
		   double f0=get(p2,a-1,b-1);
		   double f1=get(p2,a,b-1);
		   double f2=get(p2,a,b);
		   double f3=get(p2,a-1,b);
		   out=bilinear_int(xvals[a-1],xvals[a],fvals[b-1],fvals[b],f0,f1,f2,f3,xx,ff);
		   if (verbose)
		     {
		       printf("f[%i,%i,%i]=%f\n",p2,a-1,b-1,f0);
		       printf("f[%i,%i,%i]=%f\n",p2,a,b-1,f1);
		       printf("f[%i,%i,%i]=%f\n",p2,a,b,f2);
		       printf("f[%i,%i,%i]=%f\n",p2,a-1,b,f3);
		       printf("out=%f\n",out);
		     }
		 }
	       else
		 {
		   double f0=get(p2,0,b-1);
		   double f1=get(p2,0,b);
		   out=linear_int(fvals[b-1],fvals[b],f0,f1,ff);		   
		 }
	     }
	   else
	     {
	       if (a>0)
		 {
		   double f0=get(p2,a-1,0);
		   double f1=get(p2,a,0);
		   out=linear_int(xvals[a-1],xvals[a],f0,f1,xx);
		   if (verbose)
		     {
		       printf("%f %f out=%f\n",f0,f1,out);
		     }
		 }
	       else
		 out=get(p2,0,0);
	     }
	}
      else
	{
	  double f0=get(p2,a-1,m2-1);
	  double f1=get(p2,a,m2-1);
	  double f2=get(p2,a,0);
	  double f3=get(p2,a-1,0);
	  out=bilinear_int(xvals[a-1],xvals[a],fvals[m2-1],2*M_PI,f0,f1,f2,f3,xx,ff);
	  //printf("%f %f %f %f out=%f\n",f0,f1,f2,f3,out);
	  //result between last and first entry w/ periodic boundary conditions
	}
      //printf("y=%i between %f %f\n",b,fvals[b-1],fvals[b]);
    }
  else
    {
      out=def;
    }
  
  

  return out;
  
}


double vertex::fullinterpolate(double yy,double xx,double ff,double *xvals,double *fvals,int m1,int m2,int verbose,double def)
{

  
  if (verbose)
    printf("full interpolation request %f %f %f \n",yy,xx,ff);
  
  int a=0,b=0,c=0;
  double out=def;
  //make sure ff is positive
  while(ff<0)
    ff+=2*M_PI;
  //and between 0 and 2 Pi
  while(ff>2*M_PI)
    ff-=2*M_PI;

  if(verbose)
    printf("converted ff=%f\n",ff);	   

  if (yy<xvals[m1-1])
    {
      if (yy==xvals[0])
	{
	  if (verbose)
	    printf("sending to 2d interpolation\n");
	  out=interpolate(0,xx,ff,xvals,fvals,m1,m2,verbose,def);
	}
      else
	{
	  
	  while(yy>xvals[c])
	    c++;
	 
	  
	  if (xx<xvals[m1-1])
	    {
	       
	      if (xx==xvals[0])
		{
		  
		  out=interpolate(0,yy,ff,xvals,fvals,m1,m2,verbose,def);
		  if (verbose)
		    printf("sending to 2d interpolation, found %f\n",out);
		}
	      else
		{
		  while(xx>xvals[a])
		    a++;
		  //printf("x=%i between %f %f\n",a,xvals[a-1],xvals[a]);
		  if (ff<fvals[m2-1])
		    {
		      
		      while(ff>fvals[b])
			b++;
		      if (b>0)
			{
			  double f000=get(c-1,a-1,b-1);
			  double f001=get(c-1,a-1,b);
			  double f010=get(c-1,a,b-1);
			  double f011=get(c-1,a,b);
			  double f100=get(c,a-1,b-1);
			  double f101=get(c,a-1,b);
			  double f110=get(c,a,b-1);
			  double f111=get(c,a,b);
			  
			  if(verbose)
			    {
			      printf("branch a\n");
			      printf("f[%i,%i,%i]=%f\n",c-1,a-1,b-1,f000);
			      printf("f[%i,%i,%i]=%f\n",c-1,a-1,b,f001);
			      printf("f[%i,%i,%i]=%f\n",c-1,a,b-1,f010);
			      printf("f[%i,%i,%i]=%f\n",c-1,a,b,f011);
			      printf("f[%i,%i,%i]=%f\n",c,a-1,b-1,f100);
			      printf("f[%i,%i,%i]=%f\n",c,a-1,b,f101);
			      printf("f[%i,%i,%i]=%f\n",c,a,b-1,f110);
			      printf("f[%i,%i,%i]=%f\n",c,a,b,f111);
			      /*printf("f001=%f\n",f001);
			      printf("f010=%f\n",f010);
			      printf("f011=%f\n",f011);
			      printf("f100=%f\n",f100);
			      printf("f101=%f\n",f101);
			      printf("f110=%f\n",f110);
			      printf("f111=%f\n",f111);*/
			    }
		      
		      
			  out=trilinear_int(xvals[c-1],xvals[c],xvals[a-1],xvals[a],fvals[b-1],fvals[b],f000,f001,f010,f011,f100,f101,f110,f111,yy,xx,ff);
			  if (verbose)
			    {
			      printf("branch a out %f\n",out);
			      double x0=xvals[c-1];
			      double x1=xvals[c];
			      double y0=xvals[a-1];
			      double y1=xvals[a];
			      double z0=fvals[b-1];
			      double z1=fvals[b];
			      double t=(yy-x0)/(x1-x0);
			      double u=(xx-y0)/(y1-y0);
			      double v=(ff-y0)/(z1-z0);
			      //printf("inner workings x =%f %f \n",x0,x1);
			      //printf("inner workings y =%f %f \n",y0,y1);
			      //printf("inner workings z =%f %f \n",z0,z1);
			      //printf("inner workings t =%f u=%f v=%f \n",t,u,v);
			    }
			  //printf("%f %f %f %f out=%f\n",f0,f1,f2,f3,out);
			}
		      else
			{
			  double f00=get(c-1,a-1,0);
			  double f01=get(c-1,a,0);
			  double f10=get(c,a-1,0);
			  double f11=get(c,a,0);
			  out=bilinear_int(xvals[c-1],xvals[c],xvals[a-1],xvals[a],f00,f01,f10,f11,yy,xx);
			}
		    }
		  else
		    {
		      double f000=get(c-1,a-1,m2-1);
		      double f001=get(c-1,a-1,0);
		      double f010=get(c-1,a,m2-1);
		      double f011=get(c-1,a,0);
		      double f100=get(c,a-1,m2-1);
		      double f101=get(c,a-1,0);
		      double f110=get(c,a,m2-1);
		      double f111=get(c,a,0);
		      out=trilinear_int(xvals[c-1],xvals[c],xvals[a-1],xvals[a],fvals[m2-1],2*M_PI,f000,f001,f010,f011,f100,f101,f110,f111,yy,xx,ff);
		  
		      //printf("%f %f %f %f out=%f\n",f0,f1,f2,f3,out);
		      //result between last and first entry w/ periodic boundary conditions
		    }
		  //printf("y=%i between %f %f\n",b,fvals[b-1],fvals[b]);
		}
	    }
	  else
	    {
	      if (verbose)
		printf("xx=%f bigger thn max\n",yy);
	      out=def;
	    }
	}
    }
  else
    {
      if (verbose)
	printf("yy=%f bigger thn max\n",yy);
      out=def;
    }
  if(verbose)
    printf("close to %i %i %i, returning %f \n",a,b,c,out);

  return out;
  
}
