/* Fe for CV */

#include "opencubman_typdef.h"
#include "opencubman.h"


int main(int argc,char *argv[])
{
Cube cube1;
double l;
int orbn;
double fe;

if(argc==3||argc==4)
{

fe=atof(argv[2]);

orbn=0;
if(argc==4) orbn=atoi(argv[3]);


cubepointer = &cube1;

cube1.read_cube(argv[1],orbn);

cube1.sort_allocate();
cube1.sortf();

l=cube1.w_forxelectron(fe);


printf("Fe= %lf CV= %lf\n",fe,l);

cube1.deinit();
}
else{

/* Print help */

cout << " O P E N C U B M A N   P A C K A G E   version 1.00 " << endl;
cout << endl;
cout << " Fraction of electron (Fe) to Contour Value (CV) Program " << endl << endl;

cout << " This program reads G98/G03 cube file with orbital(s) " << endl;
cout << " and calculates CV corresponding to selected Fe " << endl << endl;
cout << " Usage: fetocv input.cube Fe (orbital) " << endl << endl;
cout << " where Fe is a float number, and (orbital) is an optional " << endl;
cout << " parameter specifing, which orbital from cube file to consider " << endl;
cout << " (start counting from 0, defalut is 0) " << endl;

/* end if(argc==3... */
};

/* end program */
}




