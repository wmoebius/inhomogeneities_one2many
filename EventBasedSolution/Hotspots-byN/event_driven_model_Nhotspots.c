#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <sys/stat.h>
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

double edump,etam,etas,vmean,tend,dt,lenx,leny,robs;
int nsim,Nhs,tdump,tnum,nread,seed; 

typedef struct{
  double x,y,r,tstart;
  int flag;
}rad;

typedef struct{
  double x,y;
}front;

typedef struct{
  double max,min,mean;
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
  fscanf(fin,"%d %lf \n",&Nhs,&robs);

  fclose(fin);
  tnum = floor(tend/(dt));
}

void write_parameters(int me){
  FILE *fid;
  char fname[100];

 if(me==0){
  sprintf(fname,"param.out");
  fid=fopen(fname,"w");
  fprintf(fid,"nread = %d \n",nread);
  fprintf(fid,"estart = %lf edump = %lf eta_max = %lf \n",etas,edump,etam);
  fprintf(fid,"vmean = %lf \n",vmean);
  fprintf(fid,"nsim = %d \n",nsim);
  fprintf(fid,"tend = %lf dt = %lf \n",tend,dt);
  fprintf(fid,"tdump = %d \n",tdump);
  fprintf(fid,"lenx = %lf leny = %lf \n",lenx,leny);
  fprintf(fid,"Nhs = %d, robs = %lf \n",Nhs,robs);

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

void read_random_pos(rad *radius,int me,int r)
{
  int n;
  char fname[100];
  FILE *fin;
  double xr,yr;

  sprintf(fname,"position_highdif_areas.me.%d.isample.%d.dat",me,r);
  fin = fopen(fname,"r");
  
  n = 0; 
  while(fscanf(fin,"%lf %lf ",&xr,&yr) != EOF){
  radius[n].x = xr;
  radius[n].y = yr;
  n ++;
  }
  Nhs = n; 
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

 for(i=0; i<Nhs; i++)
  {
  rnx = ( rand() / ((double)RAND_MAX+1)) * (lenx-min) + min;
  rny = ( rand() / ((double)RAND_MAX+1)) * (leny-min) + min;
  radius[i].x=rnx;
  radius[i].y=rny;

  fprintf(fidu,"%g %g \n",rnx,rny);
  }

fclose(fidu);
}

#ifdef FISHER_ZERO
void first_on(rad *radius)
{
int i,imin;
double rnxmin;

 rnxmin = lenx;

 for(i=0; i<Nhs; i++)
 {
   if(radius[i].x<rnxmin)
   {
   rnxmin = radius[i].x;
   imin = i;
   }
 }

 radius[imin].flag=1;

} 

#endif


/* function to compute distances between circlesr*/
/*void dist(rad *radius, double distr[][Nhs]){ 

int k,l;
double distx,disty;

 for(k=0;k<Nhs;k++){
  for(l=0;l<Nhs;l++){
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

 for(i=0;i<Nhs;i++){ 
 radius[i].r = 0;
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
/*  fwrite(radius,sizeof(radius),Nhs,fid);
  fclose(fid); */

for(in=0;in<Nhs;in++){
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
 for(n = 0;n<tn;n++)
 {
 fprintf(fid,"%lf %lf %lf %lf \n",n*dt,effv[n].max,effv[n].min,effv[n].mean);
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
double maxpos,minpos,meanpos; 
double toverlap,tnooverlap;
char dirname[100];

#ifdef FISHER_ZERO
int check_prop;
#endif

vel *effv;
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

radius = (rad*) malloc(sizeof(rad)*Nhs);
f = (front*) malloc(sizeof(front)*NY);

/*allocate memory for effv and rad*/

 radius = (rad*) malloc(sizeof(rad)*Nhs);
 effv = (vel*) malloc(sizeof(vel)*tnum);

/* initialize random number */


#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

/* start loop for different perturbations eta*/
for(eta=etas;eta<etam;eta=eta+edump){

fprintf(stderr,"eta = %lf \n",eta);

/* create directory*/  
  sprintf(dirname,"data.N.%d",Nhs); 
  mkdir(dirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

/* loop over nsim realization per processor*/
 for(r=0;r<nsim;r++){

  seed = me*10+time(NULL)+r;
  srand48(seed);   
  srand(seed);

  fprintf(stderr,"r = %d \n",r);

/* set velocition insdie obstacles (vr) and outside obstacle(vf)*/
  vf  = vmean-eta;
  vr = vmean+eta;
  fprintf(stderr,"vr=%g vf=%g \n",vr,vf);

 /* compute (random) positions obstaclesr*/
  if(nread==1)
  {
  read_random_pos(radius,me,r);
  vf = 1.4405;
  vr = 3.5;
  }
  else
  gen_urandom(radius,me,r,dirname);
/* initialize y-front position*/
  init_front(f); 

  timelag = robs/vr; //compute timelag caused by acceleration inside obstacle
  distrad = timelag*vf;
  bandw = vr*2; //range of obstacles t check
 
  /* initialize circles and distance between circles*/ 
  init_circles(radius);
  //double distr[Nhs][Nhs];
  //dist(radius,distr);

  /*initialize linear front*/
  if(nread==1) xlin=10.0;
  else  xlin = 0;
  noff = 0; //start of unactivated points
  
  /*loop over time as long as maximum front position smaller then lenx*/
  maxpos = 0;
  minpos = 0;
  meanpos = 0;
  tn = 0;


#ifdef FISHER_ZERO
/*if vf = 0: manually acitvate first obstacle*/
  noff = 1;
  first_on(radius);
  qsort(radius,Nhs,sizeof(rad),compare_int);
#endif


  while(tn<tnum && maxpos<lenx){ 
  t = tn*dt; 
 
  /*update linear front, starting from x=0*/
  xlin = xlin+vf*dt;

  /* update position activated points*/

#ifdef FISHER_ZERO
 check_prop=0; //check if front is still propagating: otherwise stop simulation
#endif

  for(n=0;n<noff;n++){
  if(radius[n].r<robs)
  {
  radius[n].r += vr*dt;
  #ifdef FISHER_ZERO
  	check_prop=1;
  #endif
  }
  else
  radius[n].r += vf*dt;
  }
 
 /* check which points are going to be activated*/
  fmax = vr*t+bandw;      

  for(n=noff;n<Nhs;n++){    
   if(radius[n].x<fmax){
     
     check = 0;
     s = 0;
     oncheck = radius[s].flag;
  

     xradl = radius[n].x + distrad - robs;
     if(radius[n].x-robs<0)
     xradl = radius[n].x/vr*vf;

#ifndef FISHER_ZERO
     /* check activation by linear front*/ 
     if(xlin>xradl)
     {
     check = 1;
     noff++;
     radius[n].tstart = t-(xlin-xradl)/vr;
     radius[n].flag = 1;
     if(vf!=0)
     radius[n].r = (xlin-xradl);
     }
#endif
    /*check acitivation by other scattering points*/
     while(check<1 && oncheck==1 && s<Nhs){
     
     //distn=sqrt(distr[s][n]);
    /*compute distance between points*/
    distx = (radius[s].x-radius[n].x);
    disty = (radius[s].y-radius[n].y);
    if(disty < -0.5*leny)disty = disty+leny; 
    if(disty > 0.5*leny) disty = disty-leny; 
    distn = sqrt(distx*distx+disty*disty);

     /*two situations: circles are overlapping, or not */
     toverlap = (distn/vr+radius[s].tstart);
     tnooverlap = (radius[s].tstart+2*robs/vr+(distn-2*robs)/vf);
     if(distn<=2*robs && t > toverlap) 
     {
     check = 1;   
     noff++;
     radius[n].tstart = toverlap;
     radius[n].flag = 1;
     radius[n].r = (t-toverlap)*vr;
     }
#ifndef FISHER_ZERO
     else if(distn>2*robs && t > tnooverlap)
     {
     check = 1;   
     noff++;
     radius[n].tstart = tnooverlap;
     radius[n].flag = 1;
     radius[n].r = (t-tnooverlap)*vr;
     }
#endif 

     s++; 
     oncheck = radius[s].flag;
     
 
   }
  }
 } // loop over n
	/* sort on activation*/
 
      	qsort(radius,Nhs,sizeof(rad),compare_int);
    
 

  /* find front*/
  maxpos = xlin;
  minpos = vr*t;

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
  
  effv[tn].max = (maxpos/t); 
  effv[tn].min = (minpos/t); 
  effv[tn].mean= (meanpos/t); 

  if((tn%tdump)==0)
  {
  dump_front(f,me,r,t,dirname);
  dump_radius(radius,me,r,t,dirname);
  }
 tn ++; 

#ifdef FISHER_ZERO
if(check_prop==0) maxpos = lenx+1;
#endif

 }

/* dump velocity file */
 dump_velocity(effv, me, r, tn, dirname); 

#ifdef MPI
 MPI_Barrier(MPI_COMM_WORLD);
#endif

 }



}

free(f);
free(radius);

#ifdef MPI 
 MPI_Finalize();
#endif

exit(0);

}


