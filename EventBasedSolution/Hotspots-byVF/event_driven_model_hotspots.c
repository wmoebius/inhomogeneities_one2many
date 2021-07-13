#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <sys/stat.h>

#define pi 3.14159265358979
#ifdef MPI
#include<mpi.h> 
#endif

int NY;
#define NY (800)

#ifdef MPI
int me,NP;
#else
int me=0;
#endif

double edump,etam,etas,vmean,tend,dt,lenx,leny,robs,volumef;
int nsim,unum,tdump,tnum,nread,seed; 

typedef struct{
  double x,y,r,tstart;
  int flag;
}rad;

typedef struct{
  double x,y;
}front;

typedef struct{
  double max,min,mean,std;
}vel;

void read_parameters(){

  FILE *fin;
  fin=fopen("param.in","r");
  fscanf(fin,"%d \n",&nread);
  fscanf(fin,"%lf %lf %lf \n",&etas,&edump,&etam);
  fscanf(fin,"%lf \n",&vmean);
  fscanf(fin,"%d \n",&nsim);
  fscanf(fin,"%lf %lf \n",&tend,&dt);
  fscanf(fin,"%d \n",&tdump);
  fscanf(fin,"%lf %lf \n",&lenx,&leny);
  fscanf(fin,"%lf %lf \n",&volumef,&robs);

  fclose(fin);
  tnum = floor(tend/(dt));
  unum = floor(-log(1.0-volumef)*lenx*leny/pi/robs/robs);

}

void write_parameters(int me){
  FILE *fid;
  char fname[100];

 if(me==0){
  sprintf(fname,"param.out");
  fid=fopen(fname,"w");
  fprintf(fid,"nread = %d \n",nread);
#ifndef PERFECT
  fprintf(fid,"estart = %lf edump = %lf eta_max = %lf \n",etas,edump,etam);
#else
  fprintf(fid,"eta is infinite");
#endif  
  fprintf(fid,"vmean = %lf \n",vmean);
  fprintf(fid,"nsim = %d \n",nsim);
  fprintf(fid,"tend = %lf dt = %lf \n",tend,dt);
  fprintf(fid,"tdump = %d \n",tdump);
  fprintf(fid,"lenx = %lf leny = %lf \n",lenx,leny);
  fprintf(fid,"volume fraction = %lf, robs = %lf \n",volumef,robs);

  fclose(fid);
 }

}

/* to check if flag = 0 or flag=1 */
int compare_int(const void *a,const void *b)
{
  
  rad *c_a = (rad *)a;
  rad *c_b = (rad *)b;
  
  int temp = c_b[0].flag - c_a[0].flag;
  
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

void read_random_pos(rad *radius,int me,int r, char dirname[100])
{
  int n;
  char fname[100];
  FILE *fin;
  double xr,yr;

  sprintf(fname,"%s/position_highdif_areas.me.%d.isample.%ddat",dirname,me,r);
  fin = fopen(fname,"r");
  
  n = 0; 
  while(fscanf(fin,"%lf %lf ",&xr,&yr) != EOF){
  radius[n].x = xr;
  radius[n].y = yr;
  n ++;
  }
  unum = n; 
  fclose(fin);
}

/* create random positions of robstacles*/
void gen_urandom(rad *radius, int me, int r, char dirname[100])
{

int i;
double rn,rnx,rny;
int min = 0;
char fnameu[100];
FILE *fidu;


sprintf(fnameu,"%s/position_obstacles.me.%d.r.%d.dat",dirname,me,r);
fidu = fopen(fnameu,"w");

 for(i=0; i<unum; i++)
  {
  rnx = ( rand() / ((double)RAND_MAX+1)) * (lenx-min) + min;
  rny = ( rand() / ((double)RAND_MAX+1)) * (leny-min) + min;
  radius[i].x=rnx;
  radius[i].y=rny;

  fprintf(fidu,"%g %g \n",rnx,rny);

  }

fclose(fidu);
}

/* function to compute distances between circlesr*/
/*void dist(rad *radius, double distr[][unum]){ 

int k,l;
double distx,disty;

 for(k=0;k<unum;k++){
  for(l=0;l<unum;l++){
    if(l!=k){
    distx = (radius[k].x-radius[l].x);
    disty = (radius[k].y-radius[l].y);
    if(disty < -0.5*leny)disty = disty+leny; 
    if(disty > 0.5*leny) disty = disty-leny; 
    distr[l][k]=distx*distx+disty*disty;
    }
    else
    { 
    distr[l][k]=-1;
    }
 }
  }
}*/

//initialize y-front
void  init_front(front *f){ 
 
 int j;
 double dy;
  
 dy = leny/NY;

 for(j=0;j<NY;j++){
 f[j].y = dy*j;  
 f[j].x = 0;
 }
}

// initialize circular wave
void init_circles(rad *radius){
int i;

 for(i=0;i<unum;i++){
#ifdef PERFECT 
 radius[i].r = robs;
#else
 radius[i].r = 0;
#endif
 radius[i].flag = 0;
 radius[i].tstart = 0;
  }
}

void dump_front(front *f,int me,int r,double t, char dirname[100])
{
int iy;
char fname[100];
FILE *fid;

  sprintf(fname,"%s/xfront.me.%d.r.%d.t.%.2lf.dat",dirname,me,r,t);
  fid=fopen(fname,"w");
/*  fwrite(f,sizeof(front),NY,fid);
  fclose(fid); */

for(iy=0;iy<NY;iy++){
	fprintf(fid,"%lf %lf \n",f[iy].x,f[iy].y);
 }
 fclose(fid);
}

void dump_radius(rad *radius,int me,int r,double t, char dirname[100])
{
int in;
char fname[100];
FILE *fid;

  sprintf(fname,"%s/radius.me.%d.r.%d.t.%.2lf.dat",dirname,me,r,t);
  fid=fopen(fname,"w");
/*  fwrite(radius,sizeof(radius),unum,fid);
  fclose(fid); */

for(in=0;in<unum;in++){
	fprintf(fid,"%lf %lf %lf %lf %d \n",radius[in].x,radius[in].y,radius[in].r,radius[in].tstart,radius[in].flag);
 }
 fclose(fid);

}

void dump_velocity(vel *effv, int me, int r, double tn, char dirname[100])
{
int  n;
char fname[100];
FILE *fid;

sprintf(fname,"%s/effv.me.%d.r.%d.dat",dirname,me,r);

 fid=fopen(fname,"w");
 for(n = 1;n<tn;n++)
 {
 fprintf(fid,"%lf %lf %lf %lf \n",n*dt,effv[n].max,effv[n].min,effv[n].mean);
 }
 fclose(fid);
}

void dump_meanfront(vel *effpos, int me, int r, double tn, char dirname[100])
{
int  n;
char fname[100];
FILE *fid;

sprintf(fname,"%s/meanfront.me.%d.r.%d.dat",dirname,me,r);

 fid=fopen(fname,"w");
 for(n = 1;n<tn;n++)
 {
 fprintf(fid,"%lf %lf %lf\n",n*dt,effpos[n].mean,effpos[n].std);
 }
 fclose(fid);
}

int main(int argc, char **argv){

int noff,oncheck,check;
int tn,s,r,n,j;
double vr,vf,timelag,bandw,distrad;
double xlin,fmax,xradl,t;
double distx,disty,distn,shift;
double eta,xm, xmp, ypos;
double maxpos,minpos,meanpos,sum,stdpos; 
char dirname[100];

vel *effv,*effpos;
rad *radius;
front *f;

#ifdef MPI
  MPI_Init(&argc,&argv) ;
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
#endif

/* read and write paramters*/
read_parameters();
write_parameters(me);

radius = (rad*) malloc(sizeof(rad)*unum);
f = (front*) malloc(sizeof(front)*NY);

/*allocate memory for effv and rad*/

 radius = (rad*) malloc(sizeof(rad)*unum);
 effv = (vel*) malloc(sizeof(vel)*tnum);
 effpos = (vel*) malloc(sizeof(vel)*tnum);

/* initialize random number */


#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

/* start loop for different perturbations eta*/
#ifndef PERFECT
for(eta=etas;eta<etam;eta=eta+edump){
fprintf(stderr,"eta = %.3lf \n",eta);
/* create directory*/  
  sprintf(dirname,"data.eta.%.3lf",eta); 
  mkdir(dirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

#else
eta=NAN;
fprintf(stderr,"eta = infinite \n");
/* create directory*/  
  sprintf(dirname,"data.eta.inf"); 
  mkdir(dirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

/* loop over nsim realization per processor*/
 for(r=0;r<nsim;r++){

  seed = me*10+time(NULL)+r;
  srand48(seed);   
  srand(seed);

  fprintf(stderr,"r = %d \n",r);

 /* compute (random) positions obstaclesr*/
  if(nread==1)
  {
   if(nsim>3)
   {
   fprintf(stderr,"ERROR: cannot read more then 3 samples \n");
   exit(0);
   }
   char Initdir[100];
   sprintf(Initdir,"Obs_pos"); 
   read_random_pos(radius,me,r,Initdir);
  }
  else
  gen_urandom(radius,me,r,dirname);
/* initialize y-front position*/
  init_front(f); 

  /* set velocition insdie obstacles (vr) and outside obstacle(vf)*/
  vf  = vmean-eta;
#ifndef PERFECT
  vr = vmean+eta;
#else
  vr = NAN;
#endif
  if(nread==1)
  {
  vf = 0.0073;
  vr = 0.0282;
  }

#ifndef PERFECT
  timelag = robs/vr; //compute timelag caused by acceleration inside obstacle
  bandw = vr*2; //range of obstacles t check
#else
  timelag=0.0; //if infinite filling no timelag
  bandw = vf*5; //range of obstacles t check
#endif  
  distrad = timelag*vf;
 
  /* initialize circles and distance between circles*/ 
  init_circles(radius);
  //double distr[unum][unum];
  //dist(radius,distr);

  /*initialize linear front*/
  xlin = 0;
  noff = 0; //start of unactivated points
  
  /*loop over time as long as maximum front position smaller then lenx*/
  maxpos = 0;
  minpos = 0;
  meanpos = 0;
  tn = 1;
  while(tn<tnum && maxpos<lenx){ 
  t = tn*dt; 
 
  /*update linear front, starting from x=0*/
  xlin = xlin+vf*dt;
  
 /* check which points are going to be activated*/
#ifndef PERFECT
  fmax = vr*t+bandw;
#else
  fmax = vf*10*t*bandw; //not sure about which multiplicative factor to use of vf
#endif      

  for(n=noff;n<unum;n++){    
   if(radius[n].x<fmax){
     
     check = 0;
     s = 0;
     oncheck = radius[s].flag;
  

     xradl = radius[n].x + distrad - robs;
     if(radius[n].x-robs<0){
#ifndef PERFECT
     xradl = radius[n].x/vr*vf;}
#else
     xradl = 0.0;} //what??
#endif

     /* check activation by linear front*/ 
     if(xlin>xradl)
     {
     check = 1;
     noff++;
     radius[n].tstart = t-(xlin-xradl)/vf;
     radius[n].flag = 1;
     radius[n].r+=vr/vf*(xlin-xradl);
     }

     while(check<1 && oncheck==1 && s<unum){
     
     //distn=sqrt(distr[s][n]);
    /*compute distance between points*/
    distx = (radius[s].x-radius[n].x);
    disty = (radius[s].y-radius[n].y);
    if(disty < -0.5*leny)disty = disty+leny; 
    if(disty > 0.5*leny) disty = disty-leny; 
    distn = sqrt(distx*distx+disty*disty);

     /*two situations: circles are overlapping, or not */
#ifndef PERFECT
     if(distn<=2*robs && t > (distn/vr+radius[s].tstart)) /*overlapping and activation needed*/ 
#else     
     if(distn<=2*robs && t > (radius[s].tstart)) 
#endif     
     {
     check = 1;   
     noff++;
#ifndef PERFECT
     radius[n].tstart = distn/vr+radius[s].tstart; /*write time of activation*/
#else
     radius[n].tstart = radius[s].tstart;
#endif
     radius[n].flag = 1;
     radius[n].r += (t-radius[n].tstart)*vr; /*write initial radius of activation*/
     }
#ifndef PERFECT
     else if(distn>2*robs && t > (radius[s].tstart+2*robs/vr+(distn-2*robs)/vf) ) /*non overlap and activation needed*/
#else
     else if(distn>2*robs && t > (radius[s].tstart+(distn-2*robs)/vf) )
#endif
     {
     check = 1;   
     noff++;
#ifndef PERFECT
     radius[n].tstart = radius[s].tstart+2*robs/vr+(distn-2*robs)/vf; /*write time of activation what about initial radius?*/
#else  
     radius[n].tstart = radius[s].tstart+(distn-2*robs)/vf;
#endif
     radius[n].r += (t-radius[n].tstart)*vr; /*write initial radius of activation*/ 
     radius[n].flag = 1;
     }
 
     s++; 
     oncheck = radius[s].flag;
     
 
   }
  }
 } // loop over n
	/* sort on activation*/
 
      	qsort(radius,unum,sizeof(rad),compare_int);
    
   /* update position activated points*/
   for(n=0;n<noff;n++){
   if(radius[n].r<robs)
   radius[n].r += vr*dt; //in principle it is never the case if PERFECT is working
   else
   radius[n].r += vf*dt;
   }
  

  /* find front*/
  maxpos = xlin;
#ifndef PERFECT
  minpos = vr*t;
#else
  minpos = lenx;
#endif

  for(j=0;j<NY;j++){
  
  xmp = xlin;

  for(n=0;n<noff;n++){

  ypos = f[j].y;   
  shift = 0;
  /*periodic boundary conditions*/
  if((radius[n].y-radius[n].r)<0 && (ypos-radius[n].y)>0.5*leny)
  shift = -leny;
  if((radius[n].y+radius[n].r)>leny && (radius[n].y-ypos)>0.5*leny)
  shift = leny;  

  ypos = ypos+shift;
  xm = sqrt(radius[n].r*radius[n].r-(ypos-radius[n].y)*(ypos-radius[n].y)) + radius[n].x;       

  if(xm>xmp) 
  xmp = xm;

  }

  f[j].x = xmp;

  meanpos += xmp;

  if(xmp>maxpos) /* maximum front position*/
  maxpos = xmp;  
  if(xmp<minpos) /*minimum front positoin*/
  minpos = xmp;
  }
  meanpos = meanpos/((double) NY); 
  
  sum=0;
  for(j=0;j<NY;j++){  /*use mean to calculate std*/
  f[j].x=xmp;
  sum += (xmp-meanpos)*(xmp-meanpos);
  }
  stdpos=sqrt(sum/((double) NY - 1));
 
  effv[tn].max = (maxpos/t); 
  effv[tn].min = (minpos/t); 
  effv[tn].mean= (meanpos/t); 

  effpos[tn].max=maxpos;
  effpos[tn].min=minpos;
  effpos[tn].mean=meanpos;
  effpos[tn].std=stdpos;

  /*if((tn%tdump)==0)
  {
  dump_front(f,me,r,t,dirname);
  dump_radius(radius,me,r,t,dirname);
  }*/
 tn ++; 
 }

/* dump velocity file */
 dump_velocity(effv, me, r, tn, dirname); 
/* dump meanfront file */
 dump_meanfront(effpos, me, r, tn, dirname); 


#ifdef MPI
 MPI_Barrier(MPI_COMM_WORLD);
#endif

 }


#ifndef PERFECT
} //necessary if looping on eta
#endif

free(f);
free(radius);

#ifdef MPI 
 MPI_Finalize();
#endif

exit(0);

}


