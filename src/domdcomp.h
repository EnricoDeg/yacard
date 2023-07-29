#include "data_struct.h"
/* 
 * Public Interface
 */

void domdcomp_allocate(int n_blocks);

void domdcomp_free();

void domdcomp_read_input(char *filename);

void domdcomp_init();

/* 
 * Functions to provide domain decomposition
 * info to other components
 */

int   domdcomp_get_my_block();

int   domdcomp_get_subdomain_points_full();

int * domdcomp_get_blocks_master_procs();

int * domdcomp_get_procs_ini_pos();

int   domdcomp_get_subdomain_points_x();

int   domdcomp_get_subdomain_points_y();

int   domdcomp_get_subdomain_points_z();

int * domdcomp_get_all_subdomain_points_x();

int * domdcomp_get_all_subdomain_points_y();

int * domdcomp_get_all_subdomain_points_z();

int   domdcomp_get_my_block_points_x();

int   domdcomp_get_my_block_points_y();

int   domdcomp_get_my_block_points_z();

struct t_subdomain_boundary domdcomp_get_my_block_bc();

struct t_subdomain_boundary domdcomp_get_my_block_bp();


/* 
 * Private functions
 */
static int idsd3(int i, int j, int k, int mm, int nn);