// Halo exchange
#define NRALL 0
#define NRONE 1
#define N45NO 0
#define N45go 1
#define MAX_MPI_REQUESTS 24
#define BEFORE_DERIV 0
#define BEFORE_FILTE 1
#define EXCH_ALL 0
#define EXCH_ONE 1

// State variables
#define SV_DENSITY 0
#define SV_MOMENTUM_X 1
#define SV_MOMENTUM_Y 2
#define SV_MOMENTUM_Z 3
#define SV_ENERGY 4
#define N_SV 5

// Boundary contidions
#define BC_NON_REFLECTIVE   10
#define BC_WALL_INVISCID    20
#define BC_WALL_VISCOUS     25
#define BC_INTER_CURV       30
#define BC_INTER_STRAIGHT   35
#define BC_INTER_SUBDOMAINS 40
#define BC_PERIODIC         45