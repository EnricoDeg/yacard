#include "grid_generator_regular.h"
#include "config_parser.h"
#include "data_struct.h"
#include "check.h"

static char component[] = "grid_generator";
static int lxi0;
static int let0;
static int lze0;
static double span;
static double doml0;
static double doml1;
static double domh;

static struct t_grid_fields grid_points_patch;

void grid_generator_regular_read_input(char *filename) {
    config_open(filename, component);
    lxi0  = config_read_int("lxi0");
    let0  = config_read_int("let0");
    lze0  = config_read_int("lze0");
    span  = config_read_float("span");
    doml0 = config_read_float("doml0");
    doml1 = config_read_float("doml1");
    domh  = config_read_float("domh");
    config_close();
}

void grid_generator_regular_go(int mbk, int mb, int lmx, int *mo) {
    struct t_grid_fields grid_points_global;
    
    grid_points_patch.coordinate[0].ptr = safe_malloc((lmx+1)*sizeof(double));
    grid_points_patch.coordinate[1].ptr = safe_malloc((lmx+1)*sizeof(double));
    grid_points_patch.coordinate[2].ptr = safe_malloc((lmx+1)*sizeof(double));
    
    int np = (lxi0+1) * (let0+1) * (lze0+1) - 1;
    grid_points_global.coordinate[0].ptr = safe_malloc((np+1)*sizeof(double));
    grid_points_global.coordinate[1].ptr = safe_malloc((np+1)*sizeof(double));
    grid_points_global.coordinate[2].ptr = safe_malloc((np+1)*sizeof(double));
    
    int myid;
    safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );
    // master process in block generate full block grid and then send the partition
    // to the processes in the same block (and to itself)
    if (myid == mo[mb]) {
        // generate full block grid
        for (int k=0; k<=lze0; k++) {
            for (int j=0; j<=let0; j++) {
                for (int i=0; i<=lxi0; i++) {
                	double ra0 = (doml1+doml0) / lxi0;
                    double ra1 = 2.0 * domh / let0;
                    grid_points_global.coordinate[0].ptr[i+j*(lxi0+1)+k*(lxi0+1)*(let0+1)] = 
                                                                              -doml0 + i*ra0;
                    grid_points_global.coordinate[1].ptr[i+j*(lxi0+1)+k*(lxi0+1)*(let0+1)] = 
                                                                               -domh + j*ra1;
                    grid_points_global.coordinate[2].ptr[i+j*(lxi0+1)+k*(lxi0+1)*(let0+1)] = 
                                                              span*((double)(lze0-k)/lze0-0.5);
                }
            }
        }
        // send partition
        int id_end;
        if (mb == mbk) {
        	int mpro;
        	safe_mpi( MPI_Comm_size(MPI_COMM_WORLD, &mpro) );
            id_end = mpro-1;
        } else {
            id_end = mo[mb+1];
        }
        // 1. need lpos (to be computed in domdcomp init)
        // 2. compute lio here for each receiving process
        // 3. fill buffer to be sent
        // 4. MPI_Send each coordinate
        for (int id=myid; id<=id_end; id++) {
            
        }

    } else {
        // recv partition
        safe_mpi( MPI_Recv(grid_points_patch.coordinate[0].ptr, lmx, 
        	               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD,
        	               MPI_STATUS_IGNORE) );
    
        safe_mpi( MPI_Recv(grid_points_patch.coordinate[1].ptr, lmx, 
        	               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD,
        	               MPI_STATUS_IGNORE) );
    
        safe_mpi( MPI_Recv(grid_points_patch.coordinate[2].ptr, lmx, 
        	               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD,
        	               MPI_STATUS_IGNORE) );
    
    }
    
    
    
    
    
    safe_free(grid_points_global.coordinate[0].ptr);
    safe_free(grid_points_global.coordinate[1].ptr);
    safe_free(grid_points_global.coordinate[2].ptr);
}
