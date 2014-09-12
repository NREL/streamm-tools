
/* this is cube class funtions definition */

#include "opencubman_typdef.h"

#include "opencubman.h"

Cube *cubepointer;


int comparator(const void * a, const void * b)
{
double test;
int back;
int e,f;

e=*(int*)a;
f=*(int*)b;

test= ( cubepointer->c[e].d - cubepointer->c[f].d );


if(test<0.0) back=1;
if(test>0.0) back=-1;
if(test==0.0) back=0;

return back;

}



/* read cube from filename, store orbital number orbn */
/* if orbital number = -1 then read density */
void Cube::read_cube(char *filename,int orbn)
{
fstream file1;
int i,x,y,z;
int norb;

double v;
int count,temp;

file1.open (filename);
if (file1.is_open())
  {

  file1.getline (ci.textline1,256);
  file1.getline (ci.textline2,256);

  file1 >> natoms;
  if(natoms<0) 
     {natoms=natoms*-1;
     ci.uc=0;
     } else {
     ci.uc=1;
     };

  file1 >> ci.o.x >> ci.o.y >> ci.o.z;

  file1 >> ci.nx;

  if(ci.nx<0)
     {
     ci.nx=ci.nx*-1;
     ci.um=0;
     } else {
     ci.um=1;
     };

  file1 >> ci.vx.x >> ci.vx.y >> ci.vx.z;

  file1 >> ci.ny;
  file1 >> ci.vy.x >> ci.vy.y >> ci.vy.z;
  file1 >> ci.nz;
  file1 >> ci.vz.x >> ci.vz.y >> ci.vz.z;

  ci.n=ci.nx*ci.ny*ci.nz;

  /* allocate memory for cube */

  allocate(ci.n,natoms);

  for(i=0;i<natoms;i++)
     {
     file1>> atoms[i].ian >> atoms[i].an >> atoms[i].x >> atoms[i].y >> atoms[i].z;
     }

 if(orbn>=0) {
     file1 >> norb;
     for(i=0;i<norb;i++)
        {
        file1 >> temp;
        };
     };

  count=0;     
  for(x=0;x<ci.nx;x++)
     {
     for(y=0;y<ci.ny;y++)
        {
        for(z=0;z<ci.nz;z++)
           {

           if(orbn>=0) {
             for(i=0;i<norb;i++)
	      {
	      file1 >> v;
	      
	      if(i==orbn) 
	         {
	         c[count].w=v;
		 c[count].d=v*v;
		 count++;
		 };
              };
           } else {
           file1 >> v;
	   c[count].d=v;
	   count++;
	   };


	   };
        };
     };


  file1.close();
  } else {  
  cout << "Error opening cube file " << filename << endl;
  };

/* ends read cube function */
}


/* this allocated memory for sorting */

void Cube::sort_allocate()
{
int x,y,z;
int i,j;
sort = new int [ci.n];

j=0;
for(x=0;x<ci.nx;x++)
   {
   for(y=0;y<ci.ny;y++)
      {
      for(z=0;z<ci.nz;z++)
         {
          sort[j]=j;
	  j++;
	 };
      };
   };

sort_flag=1;

/* end of sort_allocate */
}

/* this returns a cube point number for x,y,z coordinate of cube */

int Cube::get_n_forcubexyz(int ix,int iy,int iz)
{
int count;
int x,y,z;
int cuben;

  cuben=0;
  count=0;     
  for(x=0;x<ci.nx;x++)
     {
     for(y=0;y<ci.ny;y++)
        {
        for(z=0;z<ci.nz;z++)
           {
           if(ix==x&&iy==y&&iz==z) cuben=count;
	   count++;

           };
	};
     };
     
return cuben;
}


/* sorting funtion */

void Cube::sortf()
{
int i;
qsort (sort, ci.n, sizeof(int), comparator);


for(i=0;i<ci.n;i++)
   {
/* debug
printf("%d %lf\n",i,c[sort[i]].d);
*/
   if(i==0) 
      { c[sort[i]].i=c[sort[i]].d;
        ci.integ=c[sort[i]].d;
      } else {
        c[sort[i]].i=c[sort[i]].d+c[sort[i-1]].i;
	ci.integ=ci.integ+c[sort[i]].d;
      };
   };

/* end of sorting function */
}

/* get cube elements */
double Cube::get_cube_element_w(int i)
{
return c[i].w;
}

double Cube::get_cube_element_d(int i)
{
return c[i].d;
}



/* scaling density array by value val */
void Cube::scale_d(double val)
{
int i;

for(i=0;i<ci.n;i++)
   {
   c[i].d=c[i].d*val;
   };

}



/* zeroing by values */

void Cube::zero_above_w(double w)
{
int i;
for(i=0;i<ci.n;i++)
    {
    if(fabs(c[i].w)>fabs(w)) c[i].z_flag=1;
    };
}

void Cube::zero_below_w(double w)
{
int i;
for(i=0;i<ci.n;i++)
    {
    if(fabs(c[i].w)<fabs(w)) c[i].z_flag=1;
    };
}


void Cube::zero_above_d(double d)
{
int i;
for(i=0;i<ci.n;i++)
    {
    if(c[i].d>d) c[i].z_flag=1;
    };
}

void Cube::zero_below_d(double d)
{
int i;
for(i=0;i<ci.n;i++)
    {
    if(c[i].d<d) c[i].z_flag=1;
    };
}

/* get coordinates of point of cube for particular ints x,y,z */

XYZ Cube::get_cube_point(int x,int y,int z)
{
XYZ p;

p.x=ci.o.x+(double)x*ci.vx.x+(double)y*ci.vy.x+(double)z*ci.vz.x;
p.y=ci.o.y+(double)x*ci.vx.y+(double)y*ci.vy.y+(double)z*ci.vz.y;
p.z=ci.o.z+(double)x*ci.vx.z+(double)y*ci.vy.z+(double)z*ci.vz.z;

return p;
}

/* set values to zero above (d>0) or below (d<0) plane */
void Cube::zero_plane(PLANE p,int d)
{
int x,y,z;
double dist;
int j;
XYZ cpt;
double r;

r=sqrt(p.n.x*p.n.x+p.n.y*p.n.y+p.n.z*p.n.z);

  j=0;
  for(x=0;x<ci.nx;x++)
     {
     for(y=0;y<ci.ny;y++)
        {
        for(z=0;z<ci.nz;z++)
           {
           cpt=get_cube_point(x,y,z);
	   dist=p.n.x*(p.o.x-cpt.x)+p.n.y*(p.o.y-cpt.y)+p.n.z*(p.o.z-cpt.z);
	   dist=dist/r;
	   if(dist>0&&d>0) c[j].z_flag=1;
	   if(dist<0&&d<0) c[j].z_flag=1;
	   j++;
	   };
	};
     };


}





/* move cube by vector */
void Cube::movebyvector(double x,double y, double z)
{
int i;

ci.o.x=ci.o.x+x;
ci.o.y=ci.o.y+y;
ci.o.z=ci.o.z+z;

for(i=0;i<natoms;i++)
   {
   atoms[i].x=atoms[i].x+x;
   atoms[i].y=atoms[i].y+y;
   atoms[i].z=atoms[i].z+z;
   
   };

}



int Cube::n_forxelectron(double x)
{
int i;
int count;
int n;

count=0;
for(i=0;i<ci.n;i++)
   {
   if(c[sort[i]].i/ci.integ>x&&count==0)
      {
      n=i;
      if(printlevel>0) {
        printf("x= %lf n=%d percent of electron (indegrtd dens) =%lf orb_cont_value= %lf\n",x,n,c[sort[i]].i/ci.integ,(c[sort[i]].w));
        };
      count++;
      };
      
   };
return n;
}


int Cube::n_forw(double w)
{
int i;
int count;
int n;

w=sqrt(w*w);

count=0;
for(i=0;i<ci.n;i++)
   {
   if(sqrt(c[sort[i]].d)<w&&count==0)
      {
      n=i;
      if(printlevel>0) {
        printf("w= %lf n=%d percent of electron (indegrtd dens) =%lf orb_cont_value= %lf\n",w,n,c[sort[i]].i/ci.integ,(c[sort[i]].w));
        };
      count++;
      };
      
   };
return n;
}




//printf("\n Total density = %lf or %lf or %lf in Bohr\n",integ*v,integ*v/(0.5292)*(0.5292)*(0.5292),integ*v*(0.5292)*(0.5292)*(0.5292));


double Cube::d_forxelectron(double x)
{

return c[sort[n_forxelectron(x)]].d;
}

double Cube::w_forxelectron(double x)
{

return c[sort[n_forxelectron(x)]].w;
}


double Cube::x_forw(double w)
{

return c[sort[n_forw(w)]].i/ci.integ;
}


void Cube::reset_z_flag()
{
int i;

for(i=0;i<ci.n;i++)
   {
   c[i].z_flag=0;
   };

}

void Cube::write_cube(char *filename)
{
write_cube_new(filename,0,0.0);
}

void Cube::write_cube_density(char *filename)
{
write_cube_new(filename,-1,0.0);
}

/* Old routines that print cube file with orbital */
/* this is new saving option, int opt = 1 then function safes va instead of wavefunction value */
/* if opt=-1 then density is saved */
void Cube::write_cube_new(char *filename,int opt,double va)
{
int i,j,k,x,y,z;
FILE *fp1;
fp1=fopen(filename,"w");


fprintf(fp1,"\n");
if(opt==-1) {
  fprintf(fp1,"this is orbital density\n");
  } else {
  fprintf(fp1,"this is SOMO orbital\n");
  };

fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",natoms,ci.o.x,ci.o.y,ci.o.z);

fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",ci.nx,ci.vx.x,ci.vx.y,ci.vx.z);
fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",ci.ny,ci.vy.x,ci.vy.y,ci.vy.z);
fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",ci.nz,ci.vz.x,ci.vz.y,ci.vz.z);

for(i=0;i<natoms;i++)
   {
   fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf % 13.6lf\n",atoms[i].ian,atoms[i].an,atoms[i].x,atoms[i].y,atoms[i].z);
   };




/* Print SOMO density  */

if(opt!=-1) fprintf(fp1," 1    1\n");

j=0;
k=0;
for(x=0;x<ci.nx;x++)
   {
   for(y=0;y<ci.ny;y++)
      {
      for(z=0;z<ci.nz;z++)
         {
	 
	    {

	    if(opt!=-1) {
	    /* saving orbital */
	    if(c[j].z_flag==0) 
	      {
              if(opt==0) {
                  fprintf(fp1," % 13.6E ",(c[j].w));
                  } else {
	//	     printf("!");
                     if(c[j].w>0) {
                     fprintf(fp1," % 13.6E ",va);
                     } else {
		     fprintf(fp1," % 13.6E ",va*-1.0);};
                  };



              } else {
	      if(opt==0) {
	         fprintf(fp1," % 13.6E ",0.0);
                 };
              if(opt==1) {
	         fprintf(fp1," % 13.6E ",0.0);
                 };
		 /* opt==2 turns on "antialiasing" */
	      if(opt==2) {
	          if(z>1&&z<ci.nz-2)
		     {
		     if(c[j-2].z_flag==0||c[j-1].z_flag==0||c[j+1].z_flag==0||c[j+2].z_flag==0)
		        {
                        if(c[j].w>0) {
			   fprintf(fp1," % 13.6E ",0.5*va);
			   } else {
			   fprintf(fp1," % 13.6E ",va*-0.5);};
			} else {
			fprintf(fp1," % 13.6E ",0.0);
			};
		     } else {
		     fprintf(fp1," % 13.6E ",0.0);
		     };
	          };
	      };
             /* end of saving orbital */
            } else {
            /* saving density */
            if(c[j].z_flag==0) 
	       {
	       fprintf(fp1," % 13.6E ",(c[j].d));
               } else {
               fprintf(fp1," % 13.6E ",(0.0));
	       };

            /* end of saving density */
	    };

	    j++;
	    k++;
	    if(z==(ci.nz-1)) {fprintf(fp1,"\n");k=0;};
	    if(k==6) {
	             k=0;
		     fprintf(fp1,"\n");
		     };
	    };
	 };
      };
   };

/* end of cube wrtie */
}






