/* This is a source code of an example program */
/* The example demonstrates how OpenCubMan functions
   can be used to manipulate electron density
   to obtain cross-sections etc. 
   
   All comments/help is provided below, 
   within the main() function.
   
 
   Once the program is modified according to your needs
   please compile it with:
   
   g++ opencubman.cpp example.cpp -lm -o example_binary
   
   Then, execute example_binary to perform modifications
   of your cude file(s).


   If you have any questions regarding this example,
   please contact me (maharan@chem.univ.gda.pl).

   */
   
   
   
#include "opencubman_typdef.h"
#include "opencubman.h"


int main(int argc,char *argv[])
{

/* The following program read orbital from cube file
 * creates cross sections using two planes
 * and saves the resulting orbital density
 * to be visualized in external program
*/

/* define cube in program memory */
Cube cube1;

/* define points in space 
 * They will be used to define crosssection plane */
XYZ p1,p2;

/* define vectors
 * They will be normal vectors defining planes */
XYZ v1,v2;

/* define plane */
PLANE pl1,pl2;



/* read-in first orbital (0) from the cube file */
cube1.read_cube("your_cube_file.cube",0);

/* define point one */
p1.x=0.0;
p1.y=0.0;
p1.z=0.0;
/* define vector one */
v1.x=0.0;
v1.y=1.0;
v1.z=0.0;
/* define point two */
p2.x=1.0;
p2.y=0.0;
p2.z=0.0;
/* define vector two */
v2.x=0.0;
v2.y=1.0;
v2.z=0.0;

/* build plane1 from point1 and normal vector1 */
pl1.o=p1;
pl1.n=v1;

/* build plane2 from point2 and normal vector2 */
pl2.o=p2;
pl2.n=v2;


/* make crosssection through orbital using plane1 (above plane) */
cube1.zero_plane(pl1,1);

/* make second corsssection through orbital using plane2(above plane) */
cube1.zero_plane(pl2,1);

/* write resulting cube file with density to file output.cube */
cube1.write_cube_density("output.cube");


/* in some cases the orbital density is very diffuse (values too small)
 * then you can encounter problems with visualization
 * using common programs
 *
 * you can use the trick to scale the density matrix
 * you can do it as follows */

/* scale by 1000.0 */
cube1.scale_d(1000.0);
/* save scaled density into the cube file */
cube1.write_cube_density("output_scaled.cube");

/* clear zeroing flags */
cube1.reset_z_flag();

/* clear memory */
cube1.deinit();

/* end program */
}




