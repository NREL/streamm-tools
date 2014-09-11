
/* this is cube class definition */


class Cube {

private:
int natoms;
atm *atoms;
bool sort_flag;
public:

cubept *c;
sortedpt *s;
cubeinfo ci;



/* flags */

int printlevel;

/* for sorting only */
int *sort;
int *inegrated;

/* conctructor */
Cube()
  {ci.n=0;
   ci.uc=0;
   ci.um=0;
   sort_flag=0;
   printlevel=0;
   };


/* allocate */

void allocate(int n,int na)
{
int i;
c = new cubept [n];
atoms = new atm [na];

for(i=0;i<n;i++)
   {
   c[i].w=0.0;
   c[i].d=0.0;
   c[i].z_flag=0;
   };

};


/* deinit */

void deinit()
{
delete[] c;
delete[] atoms;
//if(sort_flag==1) delete[] sort;
};

/* other functions */
void read_cube(char *,int);
void sort_allocate();
void sortf();
double get_cube_element_w(int);
double get_cube_element_d(int);
void zero_above_w(double);
void zero_below_w(double);
void zero_above_d(double);
void zero_below_d(double);
XYZ get_cube_point(int,int,int);
int get_n_forcubexyz(int,int,int);

void scale_d(double);
void zero_plane(PLANE,int);
void movebyvector(double ,double, double );
int n_forxelectron(double );
int n_forw(double );
double w_forxelectron(double);
double d_forxelectron(double);
double x_forw(double);
void reset_z_flag();
void write_cube(char *);
void write_cube_new(char *,int,double);
void write_cube_density(char *);
//int comparator(const void *a , const void *b );

/* end of cube class definition */
};


extern Cube *cubepointer;


