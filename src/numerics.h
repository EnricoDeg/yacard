// numerics parameters
#define alpha (4.0/9.0)
#define beta  (1.0/36.0)
#define aa    (20.0/27.0)
#define ab    (25.0/216.0)
#define alpha10 0.09166666665902128
#define alpha12 1.3500000000438646
#define beta13  0.2666666666794667
#define a10    -0.3506944444280639
#define a12     0.8249999999974599
#define a13     0.7444444444773477
#define a14     0.014583333334139971
#define alpha01 8.000000000449464
#define beta02  6.000000000542
#define a01    -6.666666667794732
#define a02     9.000000000052879
#define a03     1.333333333178555
#define a04    -0.0833333333455545
#define lmp    11

// public interface
void numerics_allocate(int lim, int nbsize[3], int lmx);

void numerics_free();

void numerics_read_input(char *filename);

void numerics_init(int lxi, int let, int lze, int nbc[2][3], int lim);

void numerics_deriv(double *rfield, int lmx, 
	                int lxi, int let, int lze, 
	                int ijk[3][3], int nn, int nz, int m,
	                int lim);

// private
static void init_coeff();

static void fcint(double fltk, double fltr, double *alphz,
	              double *betz, double *za, double *zb, double *zc);

static void fcbcm(double fltk, double fltrbc);

static void penta(double *xu, double *xl, 
	              int is, int ie, int ns, int ne, int nt, int lim);

static void sbcco();
