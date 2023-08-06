#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "IO.h"
#include "parameters.h"
#include "check.h"
#include "config_parser.h"
#include "IO_backend_tecplot.h"
#include "utils.h"

static int backend;
static char component[] = "IO";

void IO_read_input(char *filename) {
    static char function[] = "IO_read_input";

    config_open(filename, component);
    backend = config_read_int("backend");
    config_close();
}

void IO_write_output_file() {
	static char function[] = "IO_write_output_file";

    if (backend == IO_TECPLOT) {
//        IO_write_output_file_tecplot();
    } else if (backend == IO_VTK) {
        
    } else if (backend == IO_SILO) {
        
    } else {
        finish(component, function);
    }
}

void IO_write_restart_file(int mbk, int mb, int lmx,
	                       int *mo, int *lpos,
	                       int *lxim, int *letm, int *lzem,
	                       int *lximb, int *letmb, int *lzemb, 
	                       int lxio, int leto, 
	                       int *current_step,
	                       int *ndt, double *dt, 
	                       double *dts, double *dte, double *timo, 
	                       double *restart_local_buffer) {
    
    int myid;
    safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );
    
    if (myid == mo[mb]) {
        int np = (lximb[mb]+1) * (letmb[mb]+1) * (lzemb[mb]+1) - 1;
        double *global_buffer = safe_malloc(N_SV*(np+1)*sizeof(double));
        
        // recv local data
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
            double *buffer = safe_malloc((lmx_id+1)*sizeof(double));
            for (int v=0; v<N_SV; v++) {

            	if (id==myid) {

                    int lp = lpos[id];
                    for (int k=0; k<=lzem[id]; k++) {
                        for (int j=0; j<= letm[id]; j++) {
                            int lq = lp + lio[j+k*(letm[id]+1)];
                            for (int i=0; i<= lxim[id]; i++) {
                                int l = indx3(i, j, k, 0, lxim[id], letm[id]);
                                global_buffer[lq+i+v*(np+1)] = restart_local_buffer[l+v*(lmx+1)];
                            }
                        }
                    }

            	} else {
                    safe_mpi( MPI_Recv(buffer, lmx_id+1, 
                              MPI_DOUBLE, id, id, MPI_COMM_WORLD,
                              MPI_STATUS_IGNORE) );
                    int lp = lpos[id];
                    for (int k=0; k<=lzem[id]; k++) {
                        for (int j=0; j<= letm[id]; j++) {
                            int lq = lp + lio[j+k*(letm[id]+1)];
                            for (int i=0; i<= lxim[id]; i++) {
                                int l = indx3(i, j, k, 0, lxim[id], letm[id]);
                                global_buffer[lq+i+v*(np+1)] = buffer[l];
                            }
                        }
                    }
                
                }
                
            }

            safe_free(lio);
            safe_free(buffer);
            
        }

        // write restart file
        double current_step_double = (double)(*current_step);
        double ndt_double = (double)(*ndt);

        char restart_filename[500] = "restart";
        char mb_str[5];
        snprintf(mb_str, 5, "%d", mb);
        strcat(restart_filename, mb_str);
        strcat(restart_filename, ".dat");
        
        FILE *fptr = safe_open(restart_filename, "wb");
        // restart file header
        fwrite(&current_step_double, 1, sizeof(double), fptr);
        fwrite(&ndt_double, 1, sizeof(double), fptr);
        fwrite(dt, 1, sizeof(double), fptr);
        fwrite(dts, 1, sizeof(double), fptr);
        fwrite(dte, 1, sizeof(double), fptr);
        fwrite(timo, 1, sizeof(double), fptr);
        // restart file content
        fwrite(global_buffer, 1, N_SV*(np+1)*sizeof(double), fptr);
        safe_close(fptr);
        
        safe_free(global_buffer);
        
    } else {

    	// send local data
    	for (int v=0; v<N_SV; v++)
            safe_mpi( MPI_Send(&restart_local_buffer[v*(lmx+1)], lmx+1, 
                               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD) );

    }

}

void IO_read_restart_file(int mbk, int mb, int lmx,
	                       int *mo, int *lpos,
	                       int *lxim, int *letm, int *lzem,
	                       int *lximb, int *letmb, int *lzemb, 
	                       int lxio, int leto, 
	                       int *current_step,
	                       int *ndt, double *dt, 
	                       double *dts, double *dte, double *timo, 
	                       double *restart_local_buffer) {
    
    int myid;
    safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );
    
    if (myid == mo[mb]) {
        int np = (lximb[mb]+1) * (letmb[mb]+1) * (lzemb[mb]+1) - 1;
        double *global_buffer = safe_malloc(N_SV*(np+1)*sizeof(double));

        // read restart file
        double current_step_double;
        double ndt_double;

        char restart_filename[500] = "restart";
        char mb_str[5];
        snprintf(mb_str, 5, "%d", mb);
        strcat(restart_filename, mb_str);
        strcat(restart_filename, ".dat");
        
        FILE *fptr = safe_open(restart_filename, "rb");
        // restart file header
        fread(&current_step_double, sizeof(double), 1, fptr);
        fread(&ndt_double, sizeof(double), 1, fptr);
        fread(dt, sizeof(double), 1, fptr);
        fread(dts, sizeof(double), 1, fptr);
        fread(dte, sizeof(double), 1, fptr);
        fread(timo, sizeof(double), 1, fptr);
        // restart file content
        fread(global_buffer, sizeof(double), N_SV*(np+1), fptr);
        safe_close(fptr);
        *current_step = (int)current_step_double;
        *ndt = (int)ndt_double;
        printf("%d: current step = %d\n", myid, *current_step);
        printf("%d: ndt = %d\n", myid, *ndt);
        printf("%d: dt = %f\n", myid, *dt);
        printf("%d: dts = %f\n", myid, *dts);
        printf("%d: dte = %f\n", myid, *dte);
        printf("%d: timo = %f\n", myid, *timo);
        
        // send local data
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
            double *buffer = safe_malloc((lmx_id+1)*sizeof(double));
            for (int v=0; v<N_SV; v++) {

            	if (id==myid) {

                    int lp = lpos[id];
                    for (int k=0; k<=lzem[id]; k++) {
                        for (int j=0; j<= letm[id]; j++) {
                            int lq = lp + lio[j+k*(letm[id]+1)];
                            for (int i=0; i<= lxim[id]; i++) {
                                int l = indx3(i, j, k, 0, lxim[id], letm[id]);
                                restart_local_buffer[l+v*(lmx+1)] = global_buffer[lq+i+v*(np+1)];
                            }
                        }
                    }

            	} else {
                    int lp = lpos[id];
                    for (int k=0; k<=lzem[id]; k++) {
                        for (int j=0; j<= letm[id]; j++) {
                            int lq = lp + lio[j+k*(letm[id]+1)];
                            for (int i=0; i<= lxim[id]; i++) {
                                int l = indx3(i, j, k, 0, lxim[id], letm[id]);
                                buffer[l] = global_buffer[lq+i+v*(np+1)];
                            }
                        }
                    }
            		safe_mpi( MPI_Send(buffer, lmx_id+1, 
                              MPI_DOUBLE, id, id, MPI_COMM_WORLD) );
                
                }
                
            }

            safe_free(lio);
            safe_free(buffer);
            
        }
        
        
        safe_free(global_buffer);
    } else {
        // recv local data
    	for (int v=0; v<N_SV; v++)
            safe_mpi( MPI_Recv(&restart_local_buffer[v*(lmx+1)], lmx+1, 
                               MPI_DOUBLE, mo[mb], myid, MPI_COMM_WORLD, 
                               MPI_STATUS_IGNORE) );
    }
    
}
