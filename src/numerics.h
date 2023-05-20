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

void numerics_init(int lxi, int let, int lze, int nbc[2][3], int lim, 
	               int lmx, int ijk[3][3], int mcd[2][3]);

/**
 * derivative with compact finite difference scheme
 *
 * @param[in,out] rfield     full data structure of conservative variables
 * @param[in]     nn         direction of derivative
 * @param[in]     nz         component of field
 * @param[in]     m          conservative variable to derive
 * 
 */
void numerics_deriv(double *rfield, int nn, int nz, int m);

/**
 * halo exchange before deriving or filtering
 *
 * @param[in,out] rfield     full data structure of conservative variables
 * @param[in]     nt         before derivative or filtering
 * @param[in]     nrt        one or all components
 * @param[in]     n45        extra for periodic BC
 * @param[in]     m          conservative variable for halo exchange
 * 
 */
void numerics_halo_exch(double *rfield, int nt, int nrt, int n45, int m);

// private
static void init_coeff();

static void fcint(double fltk, double fltr, double *alphz,
	              double *betz, double *za, double *zb, double *zc);

static void fcbcm(double fltk, double fltrbc);

static void penta(double *xu, double *xl, 
	              int is, int ie, int ns, int ne, int nt, int lim);

static void sbcco();
