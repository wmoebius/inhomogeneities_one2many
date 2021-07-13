/* EVENT DRIVEN MODEL FOR DIAMOND SHAPED OBSTACLES WITH ZERO VELOCITY INSIDE */ 

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <sys/stat.h>
#ifdef MPI
#include<mpi.h> 
#endif

int check_overlap();

int NY;
#define NY (1600)

#ifdef MPI
int me,NP;
#endif

double edump,etam,etas,tend,dt,lenx,leny,meanl,robs,vf;
int nsim,tdump,tnum,nread,seed,N,Nobs,Nin; 

/* the struct rad consist all inforamtion about scattering points: 
   x, y and r coordinate
   are they on or off (flag)
   are they left, right or upper corner (leri)
   the struct center has all center positions of obstacles and is used when checking whether point or line is crossing obstacle
*/

typedef struct{
  double x,y,r;
  int id,flag,leri;
}rad;

typedef struct{
  double x,y;
}center;

typedef struct{
  double x,y;
  double xval, tval; /*to store front pos which is NOT nan and the corresponding time*/
  int id,leri;
}front;

typedef struct{
  double max,min,mean;
}vel;

void read_parameters(){

  FILE *fin;
  fin=fopen("param.in","r");
  fscanf(fin,"%d \n",&nread);
  fscanf(fin,"%lf \n",&vf);
  fscanf(fin,"%d \n",&nsim);
  fscanf(fin,"%lf %lf \n",&tend,&dt);
  fscanf(fin,"%d \n",&tdump);
  fscanf(fin,"%lf %lf \n",&lenx,&leny);
  fscanf(fin,"%d %lf \n",&Nin,&robs);

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
  fprintf(fid,"vf = %lf \n",vf);
  fprintf(fid,"nsim = %d \n",nsim);
  fprintf(fid,"tend = %lf dt = %lf \n",tend,dt);
  fprintf(fid,"tdump = %d \n",tdump);
  fprintf(fid,"lenx = %lf leny = %lf \n",lenx,leny);
  fprintf(fid,"N = %d, robs = %lf \n",Nin,robs);

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

void read_random_pos(rad *radius,center *obs, int me,int r, char dirname[100])
{
  int n,s;
  char fname[100];
  FILE *fin;
  double xr,yr;

  sprintf(fname,"%s/position_highdif_areas.me.%d.isample.%d.dat",dirname,me,r);
  fin = fopen(fname,"r");
  
  s = 0;
  n = 0;
  while(fscanf(fin,"%lf %lf ",&xr,&yr) != EOF){
  obs[n].x = xr;
  obs[n].y = yr;
  radius[s].id = s;
  radius[s+1].id = s+1;
  radius[s+2].id = s+2;
  radius[s+3].id = s+3;

  radius[s].x = xr;
  radius[s+1].x = xr+robs;
  radius[s+2].x = xr;
  radius[s+3].x = xr-robs;

  radius[s].y = yr+robs;
  radius[s+1].y = yr;
  radius[s+2].y = yr-robs;
  radius[s+3].y = yr;

  radius[s].leri = 0;
  radius[s+1].leri = 1;
  radius[s+2].leri = 2;
  radius[s+3].leri = 3;

  s+=4;
  n+=1;
  }
  Nobs = n-1; 
  N = s-4; 
  fclose(fin);
}

/* create random positions of robstacles*/
void gen_urandom(rad *radius, center *obs, int me, int r, char dirname[100])
{

int i,k,count;
char fnameu[100], fnameu2[100];
double rn,rnx,rny,ynew;
FILE *fidu, *fidu2;


sprintf(fnameu,"%s/position_obstacles.me.%d.r.%d.dat",dirname,me,r);
sprintf(fnameu2,"%s/corners_obstacles.me.%d.r.%d.dat",dirname,me,r);
fidu = fopen(fnameu,"w");
fidu2 = fopen(fnameu2,"w");

count = Nin;

 for(i=0; i<4*Nin; i+=4)
  {
  rnx = ( rand() / ((double)RAND_MAX+1)) * lenx;
  rny = ( rand() / ((double)RAND_MAX+1)) * leny;

  obs[i / 4].x = rnx;
  obs[i / 4].y = rny;

  for(k=0;k<4;k++)
  {
  ynew = rny+(1-k)*robs;
  radius[i+k].x  = (k % 2)*robs+rnx;
  if(k==3){
  ynew=rny;
  radius[i+k].x=rnx-robs;
  }
  /*correction for periodic bcs*/
  /* add second obstacle on other side of domain: taking into account when checking intersection*/
#ifdef PERIODIC
  if(ynew>leny)
  { 
   ynew = ynew - leny;
   count ++;
   obs[count-1].x = rnx; 
   obs[count-1].y = rny-leny; 
  }
  if(ynew<0)
  { 
   ynew = ynew + leny;
   count ++;
   obs[count-1].x = rnx; 
   obs[count-1].y = rny+leny; 
  }
#endif
  radius[i+k].y  = ynew;
  radius[i+k].leri = k;
  radius[i+k].id = i/4;
  }
  }

for(i=0;i<count;i++)
{
  fprintf(fidu,"%g %g \n",obs[i].x,obs[i].y);
}

N = 4*Nin;
for(i=0;i<N;i++)
{
  fprintf(fidu2,"%g %g\n",radius[i].x,radius[i].y);
}

Nobs = count;
fclose(fidu);
fclose(fidu2);

}

void find_activation_points(rad *radius,center *obs)
{
int n,i,count,dcount,in,inn,on;
double xpos,ypos,xtop,ytop,ynew,ytemp;

 count = N;
 dcount = 0;
 on = -2;

 /* set flag to 1, all are active*/
 for(i=0;i<count;i++)
 { 
 radius[i].flag = 1;
 }

 for(n=0;n<N;n++)
 {
 xpos = radius[n].x;
 ypos = radius[n].y;

#ifndef PERIODIC
 /*if point is outside border, deactivate it*/
 if(ypos>leny || ypos<0){
 dcount++;
 radius[n].flag = 0;
 }
#endif 
 
 /*check if activation point is inside other obstacle*/
 in = check_overlap(obs,xpos,ypos,radius[n].id);
 
 /* if point is inside obstacle, deactivate it*/ 
  if(in==1){
  dcount++;
  radius[n].flag = 0;

  /* if top activation point already created by other corner, don't create it again*/
    //if( on != n-1 || n % 2 == 0)
    //{
    //ynew = ypos + (1-2*radius[n].leri)*robs;
    /*correction for periodic bc's*/
    //if(ynew>leny) ynew = ynew-leny;
    //if(ynew<0)    ynew = ynew+leny;

    //xtop = xpos + robs; 
    //ytop = ynew; 
    /* check if top point is in obstacle yes or no*/
    //inn = check_overlap(obs,xtop,ytop);
    //if(inn==0)
    //{ 
    //count++;
    //radius[count-1].x = xtop;
    //radius[count-1].y = ytop;
    //radius[count-1].flag = 1;
    //on = n;
    //}
    //}
  } 
}

  /* sort deactivated points and remove them*/
  qsort(radius,count,sizeof(rad),compare_int);
  N = count-dcount;

}

int check_overlap(center *obs, double xpos, double ypos, int idnum)
{

int n,in;
double dx, dy;
  
    in = 0;
    n = 0;
    while(n<Nobs && in==0)
    {
    dx = fabs(xpos-obs[n].x);
    dy = fabs(ypos-obs[n].y);
    if(dy > 0.5*leny)  dy = fabs(dy-leny); 

    if(dx / robs + dy / robs < 1 && n!=idnum)
    in=1;

    n++;
    }

   return(in);
}


int check_intersection(double xd, double yd, double xu, double yu, center *obs)
{

double x1,x2,x3,x4,y1,y2,y3,y4;
double F1,F2,F3,F4;
double xp1,xp2,yp1,yp2,xpd,xpu,ypd,ypu;
double min;
int n,in;

min = 1e-06;
in = 0;
n = 0;
while(n<Nobs && in==0)
{
/*set corners of diamond*/
x1 = obs[n].x-robs;
x2 = obs[n].x;
x3 = obs[n].x+robs;
x4 = obs[n].x;
y1 = obs[n].y;
y2 = obs[n].y+robs;
y3 = obs[n].y;
y4 = obs[n].y-robs;

F1 = (yu-yd)*x1 + (xd-xu)*y1 + (xu*yd - xd*yu);
F2 = (yu-yd)*x2 + (xd-xu)*y2 + (xu*yd - xd*yu);
F3 = (yu-yd)*x3 + (xd-xu)*y3 + (xu*yd - xd*yu);
F4 = (yu-yd)*x4 + (xd-xu)*y4 + (xu*yd - xd*yu);

/*if activation point at corner: F# is zero: correct rounding errors*/
if( fabs(F1-0) < min ) F1 = 0;
if( fabs(F2-0) < min ) F2 = 0;
if( fabs(F3-0) < min ) F3 = 0;
if( fabs(F4-0) < min ) F4 = 0;

/*check if all corners are in or outside diamond */
if((F1>=0 && F2>=0 && F3>=0 && F4>=0) || (F1<=0 && F2<=0 && F3<=0 && F4<=0))
{
in = 0;
}
else
{
xp1 = 1/sqrt(2)*x1 + 1/sqrt(2)*y1;
xp2 = 1/sqrt(2)*x2 + 1/sqrt(2)*y2;
yp1 = -1/sqrt(2)*x3 + 1/sqrt(2)*y3; 
yp2 = -1/sqrt(2)*x1 + 1/sqrt(2)*y1;

xpd = 1/sqrt(2)*xd + 1/sqrt(2)*yd;
xpu = 1/sqrt(2)*xu + 1/sqrt(2)*yu;
ypd = -1/sqrt(2)*xd + 1/sqrt(2)*yd;
ypu = -1/sqrt(2)*xu + 1/sqrt(2)*yu;

/*if activation point at corner: correct rounding errors*/

if     ( (xpd-xp1)<min && (xpu-xp1)<min)     	in = 0;
else if( (xpd-xp2)>-1*min && (xpu-xp2)>-1*min ) in = 0;
else if( (ypd-yp1)<min && (ypu-yp1)<min)     	in = 0;
else if( (ypd-yp2)>-1*min && (ypu-yp2)>-1*min ) in = 0;
else 
in = 1;

}
n++;
}

return(in);

}

int check_self_intersection(double xd, double yd, double xu, double yu, center *obs,int id)
{

double x1,x2,x3,x4,y1,y2,y3,y4;
double F1,F2,F3,F4;
double xp1,xp2,yp1,yp2,xpd,xpu,ypd,ypu;
double min;
int n,in;

min = 1e-06;
in = 0;
n = 0;
while(n<Nobs && in==0)
{
/*set corners of diamond*/
x1 = obs[n].x-robs;
x2 = obs[n].x;
x3 = obs[n].x+robs;
x4 = obs[n].x;
y1 = obs[n].y;
y2 = obs[n].y+robs;
y3 = obs[n].y;
y4 = obs[n].y-robs;

F1 = (yu-yd)*x1 + (xd-xu)*y1 + (xu*yd - xd*yu);
F2 = (yu-yd)*x2 + (xd-xu)*y2 + (xu*yd - xd*yu);
F3 = (yu-yd)*x3 + (xd-xu)*y3 + (xu*yd - xd*yu);
F4 = (yu-yd)*x4 + (xd-xu)*y4 + (xu*yd - xd*yu);

/*if activation point at corner: F# is zero: correct rounding errors*/
if( fabs(F1-0) < min ) F1 = 0;
if( fabs(F2-0) < min ) F2 = 0;
if( fabs(F3-0) < min ) F3 = 0;
if( fabs(F4-0) < min ) F4 = 0;

/*check if all corners are in or outside diamond */
if((F1>=0 && F2>=0 && F3>=0 && F4>=0) || (F1<=0 && F2<=0 && F3<=0 && F4<=0))
{
in = 0;
}
else
{
xp1 = 1/sqrt(2)*x1 + 1/sqrt(2)*y1;
xp2 = 1/sqrt(2)*x2 + 1/sqrt(2)*y2;
yp1 = -1/sqrt(2)*x3 + 1/sqrt(2)*y3; 
yp2 = -1/sqrt(2)*x1 + 1/sqrt(2)*y1;

xpd = 1/sqrt(2)*xd + 1/sqrt(2)*yd;
xpu = 1/sqrt(2)*xu + 1/sqrt(2)*yu;
ypd = -1/sqrt(2)*xd + 1/sqrt(2)*yd;
ypu = -1/sqrt(2)*xu + 1/sqrt(2)*yu;

/*if activation point at corner: correct rounding errors*/

if     ( (xpd-xp1)<min && (xpu-xp1)<min)     	in = 0;
else if( (xpd-xp2)>-1*min && (xpu-xp2)>-1*min ) in = 0;
else if( (ypd-yp1)<min && (ypu-yp1)<min)     	in = 0;
else if( (ypd-yp2)>-1*min && (ypu-yp2)>-1*min ) in = 0;
else if(n==id) in = 0;
else 
in=1;

}
n++;
}

return(in);

}

//initialize y-front
void  init_front(front *f){ 
 
 int j;
 double dy;
  
 dy = leny/NY;

 for(j=0;j<NY;j++){
 f[j].y = dy*j;  
 f[j].x = 0;
 f[j].xval = 0;
 f[j].tval = 0;
 }
}

// initialize circular wave
void init_circles(rad *radius){
int i;

 for(i=0;i<N;i++){ 
 radius[i].r = 0;
 radius[i].flag = 0;
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
	fprintf(fid,"%lf %lf %d %d \n",f[iy].x,f[iy].y,f[iy].id,f[iy].leri);
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
/*  fwrite(radius,sizeof(radius),N,fid);
  fclose(fid); */

for(in=0;in<N;in++){
	fprintf(fid,"%d %d %lf %lf %lf %d \n",radius[in].id,radius[in].leri,radius[in].x,radius[in].y,radius[in].r,radius[in].flag);
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

int main(int argc, char **argv){

int noff,oncheck,check,in;
int fid, fin, in1, in2;
int tn,s,r,n,j,jp;
int NAV;
double vr,timelag,bandw,distrad;
double xlin,fmax,xradl,t,xcheck;
double distx,disty,distn,shift,ytemp;
double eta,xm, xmp, xmf, ypos;
int    idf,lerif;
double maxpos,minpos,meanpos; 
char dirname[100];


vel    *effv;
rad    *radius;
center *obs;
front  *f;

#ifdef MPI
  MPI_Init(&argc,&argv) ;
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
#endif

/* read and write paramters*/
read_parameters();
write_parameters(me);

/*allocate memory for effv and rad*/
 obs = (center*) malloc(sizeof(center)*4*Nin);
 effv =   (vel*) malloc(sizeof(vel)*tnum);
 f =      (front*) malloc(sizeof(front)*NY);
 radius = (rad*) malloc(sizeof(rad)*4*Nin);
/* initialize random number */


#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

/* start loop for different perturbations eta*/

fprintf(stderr,"vf = %lf \n",vf);

/* create directory*/  
  sprintf(dirname,"data.Nobs.%d",Nin); 
  mkdir(dirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

/* loop over nsim realization per processor*/
 for(r=0;r<nsim;r++){
 /* malloc readius inside loop, since N can vary each realization*/

  seed = me*10+time(NULL)+r;
  srand48(seed);   
  srand(seed);

  fprintf(stderr,"r = %d \n",r);

 /* compute (random) positions obstaclesr*/
  if(nread==1)
  {
   char Initdir[100];
   sprintf(Initdir,"Obs_pos"); 
   fprintf(stderr,"WARNING: you are going to use read_random_pos, this routine is not updated for periodic bcs yet! \n");

  read_random_pos(radius,obs,me,r,Initdir);
  find_activation_points(radius,obs);
  /* realloc arrays for obstacle centers and acitvation points*/
  obs    = (center*) realloc(obs,Nobs*sizeof(center));
  radius = (rad*) realloc(radius,N*sizeof(rad));
  }
  else
  {

/* generate geometry of obsatcles and check for overlap*/
  radius = (rad*) realloc(radius,4*Nin*sizeof(rad));
  obs    = (center*) realloc(obs,4*Nin*sizeof(center));

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  gen_urandom(radius,obs,me,r,dirname);
  find_activation_points(radius,obs);
  /* realloc arrays for obstacle centers and acitvation points*/
  obs    = (center*) realloc(obs,Nobs*sizeof(center));
  radius = (rad*) realloc(radius,N*sizeof(rad));

  }

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
/* initialize y-front position*/
  init_front(f); 

  /* set velocition insdie obstacles (vr) and outside obstacle(vf)*/
  vr = 0;
  if(nread==1)
  {
  vf = 1.44;
  vr = 0;
  }

  bandw = vf*2; //range of obstacles t check
 
  /* initialize circles and distance between circles*/ 
  init_circles(radius);

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
  
  /* update position activated points*/
  for(n=0;n<noff;n++){
  radius[n].r += vf*dt;
  }

 /* check which points are going to be activated*/
  fmax = vf*t+bandw;      

  for(n=noff;n<N;n++){    
   if(radius[n].x<fmax){
     
     check = 0;
     s = 0;
     oncheck = radius[s].flag;

     /* check activation by linear front*/ 
     if(radius[n].leri==3) 
     {
     in = check_self_intersection(0,radius[n].y,xlin,radius[n].y,obs,radius[n].id); /*in=0 also if itersection is only with itself*/
     if(in==0 && xlin>=radius[n].x)
     {
     check = 1;
     noff++;
     radius[n].flag = 1;
     radius[n].r += xlin-radius[n].x;
     }
     }
     else /*for leri=0,1 or 2 usual check*/
     {
     in = check_intersection(0,radius[n].y,xlin,radius[n].y,obs); 
     if(in==0 && xlin>=radius[n].x)
     {
     check = 1;
     noff++;
     radius[n].flag = 1;
     radius[n].r += xlin-radius[n].x;
     }
     }
     
     while(check<1 && oncheck==1 && s<N){

     /* compute distance between scattering points*/ 
     shift = 0;
     distx = (radius[s].x-radius[n].x);
     disty = (radius[s].y-radius[n].y);
#ifdef PERIODIC
     if(disty < -0.5*leny)shift = leny; 
     if(disty > 0.5*leny) shift = -leny; 
#endif
     disty = disty+shift;
     distn = sqrt(distx*distx+disty*disty);

    /* check if radial wave can reach acitvation point */

#ifdef PERIODIC    
    /*periodic boundary conditions*/ 
    if(shift!=0)
    {
    in1 = check_intersection(radius[s].x,radius[s].y,radius[n].x,radius[n].y+shift,obs); 
    in2 = check_intersection(radius[s].x,radius[s].y-shift,radius[n].x,radius[n].y,obs); 
    if(in1==0 && in2==0) in = 0;
    }
    else
#endif
    in = check_intersection(radius[s].x,radius[s].y,radius[n].x,radius[n].y,obs); 
 
     if(in==0 && radius[s].r >= distn) 
     {
     check = 1;   
     noff++;
     radius[n].flag = 1;
     radius[n].r =+ radius[s].r-distn;
     }
 
     s++; 
     oncheck = radius[s].flag;
   }
  }
 } // loop over n

	/* sort on activation*/
    	qsort(radius,N,sizeof(rad),compare_int);
    
  

  /* find front*/
  NAV = 0;       /*correct number for averaging (NY-front position wiht nan)*/
  meanpos = 0;
  maxpos = 0;    /*initialize to minimim front*/ 
  minpos = vf*t; /*initialize to maximum front pos*/

  for(j=0;j<NY;j++)
  {
  
  ypos = f[j].y;   
  xmp = 0;
 
  /* check if front can be reached by linear front */ 
  in = check_intersection(0,ypos,xlin,ypos,obs); 
  if(in==0)
  {
  xmp = xlin;
  idf = -1;
  lerif = -1;
  }
  else
  { 
  for(n=0;n<noff;n++)
  {

  shift = 0;
#ifdef PERIODIC
  /*periodic boundary conditions*/
  if((radius[n].y-radius[n].r)<0 && (ypos-radius[n].y)>0.5*leny)
  shift = -leny;
  if((radius[n].y+radius[n].r)>leny && (radius[n].y-ypos)>0.5*leny)
  shift = leny;  
#endif
  ytemp = ypos+shift;
  xm = sqrt(radius[n].r*radius[n].r-(ytemp-radius[n].y)*(ytemp-radius[n].y)) + radius[n].x;       

  /* check if wave can be reached by radial front*/ 
   if(xm>xmp) 
   {
#ifdef PERIODIC
    /* periodic boundary conditions*/
    if(shift!=0)
    {
    in1 = check_intersection(radius[n].x,radius[n].y,xm,ytemp,obs);
    in2 = check_intersection(radius[n].x,radius[n].y-shift,xm,ypos,obs);
    if(in1==0 && in2==0) in = 0;
    }
    else
#endif
    in = check_intersection(radius[n].x,radius[n].y,xm,ypos,obs);

    if(in==0) 
    {
    xmp = xm;
    idf = radius[n].id;
    lerif = radius[n].leri;
    }

   if(xmp>xlin) 
   fprintf(stderr,"maxpos = %lf t = %lf \n",maxpos,t); 
   } 
  }

  }

  /* check if front is inside obstacle */
  xcheck =f[j].xval+vf*(t-f[j].tval);
  if(xcheck>xlin) xcheck = xlin;
  in=check_intersection(xmp,ypos,xcheck,ypos,obs); 
  if(in==1) 
  {
   xmp = NAN; 
   idf = -2;
   lerif = -2; 
  }

  f[j].x = xmp;
  f[j].id = idf;
  f[j].leri = lerif;
  
   if(isnan(xmp)==0)
   {
   f[j].xval = xmp;
   f[j].tval = t;
   NAV++;
   if(xmp>xlin) 
   fprintf(stderr,"maxpos = %lf t = %lf \n",maxpos,t); 

   meanpos += xmp;
   if(xmp>maxpos) /* maximum front position*/
   maxpos = xmp;  
   if(xmp<minpos) /*minimum front positoin*/
   minpos = xmp;
   }

  } /*loop over j=0..NY*/

  meanpos = meanpos/((double) NAV); 
  if(NAV==0) meanpos=0;

  effv[tn].max = (maxpos/t); 
  effv[tn].min = (minpos/t); 
  effv[tn].mean= (meanpos/t); 

  if((tn%tdump)==0)
  {
  dump_front(f,me,r,t,dirname);
  dump_radius(radius,me,r,t,dirname);
  }

/* exit when front cannot propagate anymore*/
  if(NAV==0) tn = tnum;

 tn ++; 
 }

/* dump velocity file */
 dump_velocity(effv, me, r, tn, dirname); 


#ifdef MPI
 MPI_Barrier(MPI_COMM_WORLD);
#endif

 }

free(radius);
free(obs);
free(effv);
free(f);

#ifdef MPI 
 MPI_Finalize();
#endif

exit(0);

}


