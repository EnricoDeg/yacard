#include <stdio.h>
#include <stdlib.h>
#include "numerics.h"
#include "parameters.h"
#include "check.h"
#include "domdcomp.h"
#include "grid_generator.h"

int main() {
    safe_mpi( MPI_Init(NULL,NULL) );

    char filename[] = "canard.cfg";

    int nblocks = 0;
    
    /* 
     * Components allocation
     */
//    numerics_allocate(lim, nbsize, lmx);
    domdcomp_allocate(nblocks);

    /* 
     * Components read input
     */
//    numerics_read_input(filename);
    domdcomp_read_input(filename);
    grid_generator_read_input(filename);

    /* 
     * Components initialization
     */
    domdcomp_init();
//    numerics_init(lxi, let, lze, nbc, lim, lmx, ijk, mcd);

    /* 
     * Generate mesh
     */
    int  mb   = domdcomp_get_my_block();
    int  lmx  = domdcomp_get_subdomain_points_full();
    int *mo   = domdcomp_get_blocks_master_procs();
    int *lpos = domdcomp_get_procs_ini_pos();
    int *lxim = domdcomp_get_all_subdomain_points_x();
    int *letm = domdcomp_get_all_subdomain_points_y();
    int *lzem = domdcomp_get_all_subdomain_points_z();
    int  lxio = domdcomp_get_my_block_points_x();
    int  leto = domdcomp_get_my_block_points_y();
    grid_generator_go(nblocks, mb, lmx, mo, lpos, lxim, letm, lzem, lxio, leto);

    /* 
     * Numerics calculate derivative
     */
//    numerics_deriv(rfield, nn, nz, m);
    
    /* 
     * Components destroy
     */
//    numerics_free();
    domdcomp_free();

    safe_mpi( MPI_Finalize() );
    
    return 0;
}