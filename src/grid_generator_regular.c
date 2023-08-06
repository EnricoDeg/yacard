#include "grid_generator_regular.h"
#include "config_parser.h"
#include "data_struct.h"
#include "check.h"
#include "utils.h"

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

double* grid_generator_regular_get_x() {
    return grid_points_patch.coordinate[0].ptr;
} 

double* grid_generator_regular_get_y() {
    return grid_points_patch.coordinate[1].ptr;
} 

double* grid_generator_regular_get_z() {
    return grid_points_patch.coordinate[2].ptr;
} 

void grid_generator_regular_go(int mbk, int mb, int lmx, int *mo, int *lpos, int *lxim,
                               int *letm, int *lzem, int lxio, int leto) {
    struct t_grid_fields grid_points_global;
    
    grid_points_patch.coordinate[0].ptr = safe_malloc((lmx+1)*sizeof(double));
    grid_points_patch.coordinate[1].ptr = safe_malloc((lmx+1)*sizeof(double));
    grid_points_patch.coordinate[2].ptr = safe_malloc((lmx+1)*sizeof(double));
    
    int myid;
    safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );
    // master process in block generate full block grid and then send the partition
    // to the processes in the same block (and to itself)
    if (myid == mo[mb]) {
    	int np = (lxi0+1) * (let0+1) * (lze0+1) - 1;
        grid_points_global.coordinate[0].ptr = safe_malloc((np+1)*sizeof(double));
        grid_points_global.coordinate[1].ptr = safe_malloc((np+1)*sizeof(double));
        grid_points_global.coordinate[2].ptr = safe_malloc((np+1)*sizeof(double));
        
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

        for (int id=myid; id<=id_end; id++) {

            int lmx_id = ( lxim[id] + 1 ) * ( letm[id] + 1 ) * ( lzem[id] + 1 ) - 1;
            int *lio = safe_malloc((letm[id]+1)*(lzem[id]+1)*sizeof(int));

            // compute lio for each receiving process
            for (int k=0; k<=lzem[id]; k++) {
                int kp = k * (leto+1) * (lxio+1);
                for (int j=0; j<= letm[id]; j++) {
                    int jp = j * (lxio+1);
                    lio[j+k*(letm[id]+1)] = jp+kp;
                }
            }

            // fill buffer
            double *buffer_x = safe_malloc((lmx_id+1)*sizeof(double));
            double *buffer_y = safe_malloc((lmx_id+1)*sizeof(double));
            double *buffer_z = safe_malloc((lmx_id+1)*sizeof(double));
            
            int lp = lpos[id];
            for (int k=0; k<=lzem[id]; k++) {
                for (int j=0; j<= letm[id]; j++) {
                    int lq = lp + lio[j+k*(letm[id]+1)];
                    for (int i=0; i<= lxim[id]; i++) {
                        int l = indx3(i, j, k, 0, lxim[id], letm[id]);
                        buffer_x[l] = grid_points_global.coordinate[0].ptr[lq+i];
                        buffer_y[l] = grid_points_global.coordinate[1].ptr[lq+i];
                        buffer_z[l] = grid_points_global.coordinate[2].ptr[lq+i];
                    }
                }
            }

            // send each coordinate
            if (id==myid) {
                for (int k=0; k<=lzem[id]; k++) {
                    for (int j=0; j<= letm[id]; j++) {
                        for (int i=0; i<= lxim[id]; i++) {
                            int l = indx3(i, j, k, 0, lxim[id], letm[id]);
                            grid_points_patch.coordinate[0].ptr[l] = buffer_x[l];
                            grid_points_patch.coordinate[1].ptr[l] = buffer_y[l];
                            grid_points_patch.coordinate[2].ptr[l] = buffer_z[l];
                        }
                    }
                }
            } else {
                safe_mpi(MPI_Send(buffer_x, lmx_id+1, MPI_DOUBLE, id, id, MPI_COMM_WORLD));
                safe_mpi(MPI_Send(buffer_y, lmx_id+1, MPI_DOUBLE, id, id, MPI_COMM_WORLD));
                safe_mpi(MPI_Send(buffer_z, lmx_id+1, MPI_DOUBLE, id, id, MPI_COMM_WORLD));
            }

            safe_free(buffer_x);
            safe_free(buffer_y);
            safe_free(buffer_z);

            safe_free(lio);

        }

        safe_free(grid_points_global.coordinate[0].ptr);
        safe_free(grid_points_global.coordinate[1].ptr);
        safe_free(grid_points_global.coordinate[2].ptr);

    } else {
        // recv partition
        safe_mpi( MPI_Recv(grid_points_patch.coordinate[0].ptr, lmx+1, 
        	               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD,
        	               MPI_STATUS_IGNORE) );
    
        safe_mpi( MPI_Recv(grid_points_patch.coordinate[1].ptr, lmx+1, 
        	               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD,
        	               MPI_STATUS_IGNORE) );
    
        safe_mpi( MPI_Recv(grid_points_patch.coordinate[2].ptr, lmx+1, 
        	               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD,
        	               MPI_STATUS_IGNORE) );
    
    }

}
