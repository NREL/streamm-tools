/* this is cube class definition */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAXATOM 200

#include <iostream>
#include <fstream>
using namespace std;

typedef struct
{
double x,y,z;
} XYZ;

typedef struct
{
XYZ o;
XYZ n;
} PLANE;

typedef struct
{
double w; //value of orbital
double d; //value of density
double i; //value of integrated quantity
bool z_flag; // this flag is set to 1 if element is set to zero
} cubept;

typedef struct
{
int order; // used to store order of sorted points
double i;//integrated density from highest point to this one
} sortedpt;

typedef struct
{
XYZ o;
XYZ vx,vy,vz;
int nx,ny,nz;
int n; //total points
int uc; // units for cube 0 for Bohr, 1 - for A
int um; // units for molecules
double integ; //indegrated density
char textline1[256];
char textline2[256];
} cubeinfo;

typedef struct
{
int ian;
double an;
double x,y,z;

}atm;



