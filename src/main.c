#include <stdio.h>
#include <stdlib.h>
#include "numerics.h"
#include "parameters.h"
#include "check.h"
#include "domdcomp.h"
#include "grid_generator.h"
#include "IO.h"

int main() {
    
    safe_mpi( MPI_Init(NULL,NULL) );
    
    char filename[] = "canard.cfg";
    
    int nblocks = 0;
    
    /* 
     * Domain decomposition module allocation
     */
    domdcomp_allocate(nblocks);
    
    /* 
     * Components read input
     */
    numerics_read_input(filename);
    domdcomp_read_input(filename);
    grid_generator_read_input(filename);
    
    /* 
     * Domain decomposition initialization
     */
    domdcomp_init();
    
    /* 
     * Numerics module allocation
     */
//    numerics_allocate(lim, nbsize, lmx);
    
    int lxi = domdcomp_get_subdomain_points_x();
    int let = domdcomp_get_subdomain_points_y();
    int lze = domdcomp_get_subdomain_points_z();
    int lmx = domdcomp_get_subdomain_points_full();
    struct t_subdomain_boundary nbc = domdcomp_get_my_block_bc();
    struct t_subdomain_boundary mcd = domdcomp_get_my_block_bp();
    
    /* 
     * Numerics initialization
     */
//    numerics_init(lxi, let, lze, nbc, lim, lmx, mcd);
    
    /* 
     * Generate mesh
     */
    int  mb   = domdcomp_get_my_block();
    int *mo   = domdcomp_get_blocks_master_procs();
    int *lpos = domdcomp_get_procs_ini_pos();
    int *lxim = domdcomp_get_all_subdomain_points_x();
    int *letm = domdcomp_get_all_subdomain_points_y();
    int *lzem = domdcomp_get_all_subdomain_points_z();
    int  lxio = domdcomp_get_my_block_points_x();
    int  leto = domdcomp_get_my_block_points_y();
    grid_generator_go(nblocks, mb, lmx, mo, lpos, lxim, letm, lzem, 
                      lxio, leto);

    /* 
     * Read restart files
     */
    int *lximb = domdcomp_get_blocks_points_x();
    int *letmb = domdcomp_get_blocks_points_y();
    int *lzemb = domdcomp_get_blocks_points_z();
    int current_step, ndt;
    double dt, dts, dte, timo;
    double *restart_local_buffer = safe_malloc(N_SV*(lmx+1)*sizeof(double));
    IO_read_restart_file(nblocks, mb, lmx, mo, lpos,
                         lxim, letm, lzem,
                         lximb, letmb, lzemb, 
                         lxio, leto, 
                         &current_step, &ndt, &dt, 
                         &dts, &dte, &timo, 
                         restart_local_buffer);
    
    /* 
     * Numerics calculate derivative
     */
//    numerics_deriv(rfield, nn, nz, m);
    
    /* 
     * Write restart files
     */
     current_step = 1.0;
     ndt  = 2.0;
     dt   = 3.0;
     dts  = 4.0;
     dte  = 5.0;
     timo = 6.0;
//    double *restart_local_buffer = safe_malloc(N_SV*(lmx+1)*sizeof(double));
    for (int v=0, l=0; v<N_SV; v++)
        for (int i=0; i<=lmx+1; i++, l++)
            restart_local_buffer[l] = i+v;
    IO_write_restart_file(nblocks, mb, lmx, mo, lpos,
                          lxim, letm, lzem,
                          lximb, letmb, lzemb, 
                          lxio, leto, 
                          &current_step, &ndt, &dt, 
                          &dts, &dte, &timo, 
                          restart_local_buffer);
    safe_free(restart_local_buffer);
    
    /* 
     * Components destroy
     */
//    numerics_free();
    domdcomp_free();
    
    safe_mpi( MPI_Finalize() );
    
    return 0;
    
}
