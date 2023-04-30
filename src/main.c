#include <stdio.h>
#include <stdlib.h>
#include "numerics.h"
#include "parameters.h"
#include "check.h"

int main() {
	// values to be initialize for the dwarf
    int lim;
    int nbsize[3];
    int lmx;
    int lxi, let, lze, nbc[2][3], ijk[3][3];
    char filename[] = "canard.cfg";
    double *rfield;
    int nn, nz, m;
    
    nn = 0;
    nz = 0;
    m  = 0;
    
    lxi = 12;
    let = 12;
    lze = 12;
    lim = ( lxi + 1 ) + ( let + 1 ) + ( lze + 1 ) - 1;
    lmx = ( lxi + 1 ) * ( let + 1 ) * ( lze + 1 ) - 1;
    
    ijk[0][0] = lxi;
    ijk[1][0] = let;
    ijk[2][0] = lze;
    ijk[0][1] = let;
    ijk[1][1] = lze;
    ijk[2][1] = lxi;
    ijk[0][2] = lze;
    ijk[1][2] = lxi;
    ijk[2][2] = let;
    
    nbsize[0] = 169;
    nbsize[1] = 169;
    nbsize[2] = 169;
    
    for (int i=0; i<2; i++)
        for (int j=0; j<3; j++)
            nbc[i][j] = BC_NON_REFLECTIVE;
    
    rfield = safe_malloc(3*(lmx+1)*sizeof(double));
    for (int j=0; j<3; j++)
        for (int i=0; i<lmx+1; i++)
            rfield[i+j*(lmx+1)] = 1.0;
    
    // numerics allocation
    numerics_allocate(lim, nbsize, lmx);
    // numerics input parameters
    numerics_read_input(filename);
    // numerics initialization
    numerics_init(lxi, let, lze, nbc, lim);
    // numerics calculate derivative
    numerics_deriv(rfield, lmx, lxi, let, lze, ijk, nn, nz, m, lim);
    
    printf("rfield\n");
    for (int i=0; i<lmx+1; i++)
        printf("%d, %f\n", i, rfield[i]);
    
    // numerics destroy
    numerics_free();
    
    safe_free(rfield);
    
    return 0;
}