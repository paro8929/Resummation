struct coord
{
  int d;//dimensions
  int* x;//coordinates
  bool dat_alloc;//flag to check if memory has been allocated
 
  void allocate()
  {
    x = new int[d];
    dat_alloc=1;
  }

  void free()
  {
    delete [] x;
    dat_alloc=0;
  }

  void print()
  {
    if (dat_alloc)
      {
	for (int i=0;i<d;i++)
	  printf("x[%i]=%i\t",i,x[i]);
	printf("\n");
      }
    else
      printf("Error, memory not allocated\n");
  }
  
  void set(int* a)
  {
    if (dat_alloc)
      {
	for (int i=0;i<d;i++)
	  x[i]=a[i];
      }    
    else
      printf("Error, memory not allocated\n");
  }
  
  double veclen()
  {
    double sum=0;
     if (dat_alloc)
       {
	 for (int i=0;i<d;i++)
	   sum+=x[i]*x[i];	 
       }
     else
      printf("Error, memory not allocated\n");
     return sum;
  }

};


class data_array
{
 public:
  struct coord gridsize;  //contains info on x,y,z,... extent of grid
  struct coord Middle; //coordinates of middle point


  long int SIZE; //total number of data points
  long int intMiddle; //middle point

  long int *helparr;


  //generic data pointer
   //should contain double array of size SIZE 
   double * data;

   bool dat_alloc;

   data_array()
     {
       dat_alloc=0;
     };

   ~data_array()
     {
       //printf("===> Info :  Deleting data_array\n");
       //if we allocated data, free memory
       delete [] helparr;
       if (dat_alloc)
	 delete [] data;
       dat_alloc=0;
     };


   void init(struct coord tdim);
   void alloc_data()
   {
     data = new double[SIZE];
     dat_alloc=1;
   }

   long int getind(struct coord *co);
   void getcoord(struct coord *co, long int ind);
   long int getplus(long int ind,int dir);
   long int getminus(long int ind,int dir);
};


void data_array::init(struct coord tdim)
{
  SIZE=1;
  gridsize.d=tdim.d;
  gridsize.allocate();

  Middle.d=tdim.d;
  Middle.allocate();

  helparr=new long int[tdim.d];


  for (int i=0;i<tdim.d;i++)
    {
      gridsize.x[i]=tdim.x[i];
      Middle.x[i]=(gridsize.x[i]-1)/2;
      SIZE*=tdim.x[i];
    }
  
  //Middle.print();

  for (int i=0;i<gridsize.d;i++)
    {
      helparr[i]=1;
      for (int j=i+1;j<gridsize.d;j++)
	helparr[i]*=gridsize.x[j];
      //printf("helparr=%li\n",helparr[i]);
    }

  
  
  printf("===> Info: Generating data_array\n");
  for (int i=0;i<tdim.d;i++)
    printf("Dimension %i: %i points\n",i,tdim.x[i]);  
 
  dat_alloc=0;

  intMiddle=getind(&Middle);
}


long int data_array::getind(struct coord *co)
{
  long int temp=0;
  //long int help;
  for (int i=0;i<gridsize.d;i++)
    {
      //help=1;
      // for (int j=i+1;j<gridsize.d;j++)
      //	help*=gridsize.x[j];
      
      //printf("help =%li helparr=%li\n",help,helparr[i]);

      temp+=co->x[i]*helparr[i];
    }
  return temp;
}


//get neighbor in positive direction dir
long int data_array::getplus(long int ind, int dir)
{


  long int temp;
  
  if (dir==0)
    {
      //temp=ind/helparr[0];
      temp=ind;
    }
  else
    {
      temp=ind%helparr[dir-1];
      //temp/=helparr[dir];
    }

  int neighbor=temp;

  neighbor+=helparr[dir]; //neighbor in positive direction
  //neighbor%=gridsize.x[dir]; //periodic boundary conditions
  if (neighbor>=gridsize.x[dir]*helparr[dir])
    {
      //neighbor-=gridsize.x[dir];
      //printf("check one %li vs %li\n",(ind+(neighbor-temp)*helparr[dir]),ind+(1-gridsize.x[dir])*helparr[dir]);
      return (ind+(1-gridsize.x[dir])*helparr[dir]);
    }
  else
    {
      return (ind+helparr[dir]);
    }
  
  
  //return (ind+(neighbor-temp)*helparr[dir]);
}

//get neighbor in negative direction dir
long int data_array::getminus(long int ind, int dir)
{


  long int temp;
  if (dir==0)
    {
      temp=ind/helparr[0];
    }
  else
    {
      temp=ind%helparr[dir-1];
      temp/=helparr[dir];
    }

  int neighbor=temp;

  neighbor--; //neighbor in negative direction
  neighbor%=gridsize.x[dir]; //periodic boundary conditions

  if(neighbor<0)
    neighbor+=gridsize.x[dir];

  return (ind+(neighbor-temp)*helparr[dir]);
}


void data_array::getcoord(struct coord *co, long int ind)
{
  long int temp;
  co->x[0]=ind/helparr[0];
  temp=ind%helparr[0];
  for (int i=1;i<gridsize.d;i++)
    {
      //printf("helparr=%li\n",helparr[i]);
      co->x[i]=temp/helparr[i];
      temp%=helparr[i];
    }
  
}


