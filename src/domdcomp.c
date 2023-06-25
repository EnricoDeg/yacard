#include "domdcomp.h"
#include "check.h"
#include "config_parser.h"
#include "parameters.h"

static int nblocks;
static int mpro;
static int lxio, leto, lzeo;
static int lxi, let, lze;
static int mb;
static int lmx;

static int ijkp[3];
static int nbsize[3];
static int mmcd[3][2];
static int nbc[3][2];
static int mcd[3][2];
static int ijk[3][3];

static int *nbpc;
static int *mo;
static int *lximb;
static int *letmb;
static int *lzemb;
static int *lxim;
static int *letm;
static int *lzem;
static int *nbbc;
static int *mbcd;
static int *imjp;
static int *jptag;
static int *imjl;
static int *jltag;
static int *jjp;
static int *jjl;
static char component[] = "domdcomp";

void domdcomp_allocate(int n_blocks) {

    safe_mpi( MPI_Comm_size(MPI_COMM_WORLD, &mpro) );
    mpro--;

	nblocks = n_blocks;
	int lpp = 8*nblocks+7;
	int lqq = 12*nblocks+11;
    nbpc    = safe_malloc(3*(nblocks+1)*sizeof(int));
    mo      = safe_malloc((nblocks+1)*sizeof(int));
    lximb   = safe_malloc((nblocks+1)*sizeof(int));
    letmb   = safe_malloc((nblocks+1)*sizeof(int));
    lzemb   = safe_malloc((nblocks+1)*sizeof(int));
    lxim    = safe_malloc((mpro+1)*sizeof(int));
    letm    = safe_malloc((mpro+1)*sizeof(int));
    lzem    = safe_malloc((mpro+1)*sizeof(int));
    nbbc    = safe_malloc(2*3*(nblocks+1)*sizeof(int));
    mbcd    = safe_malloc(2*3*(nblocks+1)*sizeof(int));
    imjp    = safe_malloc((lpp+1)*sizeof(int));
    jptag   = safe_malloc((lpp+1)*sizeof(int));
    imjl    = safe_malloc((lqq+1)*sizeof(int));
    jltag   = safe_malloc((lqq+1)*sizeof(int));
    jjp     = safe_malloc((lpp+1)*sizeof(int));
    jjl     = safe_malloc((lqq+1)*sizeof(int));
}

void domdcomp_free() {
    safe_free(nbpc);
    safe_free(mo);
    safe_free(lximb);
    safe_free(letmb);
    safe_free(lzemb);
    safe_free(lxim);
    safe_free(letm);
    safe_free(lzem);
    safe_free(nbbc);
    safe_free(mbcd);
    safe_free(imjp);
    safe_free(jptag);
    safe_free(imjl);
    safe_free(jltag);
    safe_free(jjp);
    safe_free(jjl);
}

void domdcomp_read_input(char *filename) {
    config_open(filename, component);
    config_read_array_int("domdcomp.nbpc", nbpc);
    config_read_array_int("domdcomp.lximb", lximb);
    config_read_array_int("domdcomp.letmb", letmb);
    config_read_array_int("domdcomp.lzemb", lzemb);
    config_read_array_int("domdcomp.nbbc", nbbc);
    config_read_array_int("domdcomp.mbcd", mbcd);
    config_close();
}

void domdcomp_init() {
    int myid;
    
    safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );

    mo[0] = 0;
    for (int i=1; i<=nblocks; i++)
        mo[i] = mo[i-1] + nbpc[i-1]
                        + nbpc[i-1+(nblocks+1)]
                        + nbpc[i-1+2*(nblocks+1)];

    for (int i=0; i<=nblocks; i++) {
        if (myid >= mo[i])
            mb = i;
    }

    lxio = lximb[mb];
    leto = letmb[mb];
    lzeo = lzemb[mb];

    ijkp[0] = (myid - mo[mb])%(nbpc[mb]);
    ijkp[1] = (myid - mo[mb] / 
    	              nbpc[mb])%(nbpc[mb+(nblocks+1)]);
    ijkp[2] = (myid - mo[mb] / 
    	             (nbpc[mb]*nbpc[mb+(nblocks+1)]))%
                     (nbpc[mb+2*(nblocks+1)]);

    for (int nn=0; nn<3; nn++) {
        int nstart = ((nn+1)%3);
        int nend   = ((nstart+1)%3);
        for (int ip=0; ip<2; ip++) {
            int mm = mbcd[mb+nn*(nblocks+1)+ip*3*(nblocks+1)];
            if (mm == -1) {
                mmcd[nn][ip] = -1;
            } else {
                mmcd[nn][ip] = idsd3((1 - ip) * (nbpc[mm+nn*(nblocks+1)] - 1),
                                     ijkp[nstart],
                                     ijkp[nend],
                                     mm,
                                     nn);
            }
        }
    }

    for (int nn=0; nn<3; nn++) {
    	int ll, mp;
        switch(nn) {
            case 0:
                ll = lxio;
                mp = 1;
                break;
            case 1:
                ll = leto;
                mp = nbpc[mb];
                break;
            case 2:
                ll = lzeo;
                mp = nbpc[mb] * nbpc[mb+(nblocks+1)];
                break;
        }
        int lp = ijkp[nn];
        int ma = nbpc[mb+nn*(nblocks+1)];
        int l;
        if (ma == 1) {
            l = ll;
            nbc[nn][0] = nbbc[mb+nn*(nblocks+1)];
            nbc[nn][1] = nbbc[mb+nn*(nblocks+1)+1*3*(nblocks+1)];
            mcd[nn][0] = mmcd[nn][0];
            mcd[nn][1] = mmcd[nn][1];
        }
        
        if (ma > 1) {
            if (lp == 0) {
                l = ll - ( ( ll + 1 ) / ma ) * ( ma - 1 );
                nbc[nn][0] = nbbc[mb+nn*(nblocks+1)];
                nbc[nn][1] = BC_INTER_SUBDOMAINS;
                mcd[nn][0] = mmcd[nn][0];
                mcd[nn][1] = myid + mp;
            }
            if (lp > 0 && lp < ma-1) {
                l = ( ll + 1 ) / ma - 1;
                nbc[nn][0] = BC_INTER_SUBDOMAINS;
                nbc[nn][1] = BC_INTER_SUBDOMAINS;
                mcd[nn][0] = myid - mp;
                mcd[nn][1] = myid + mp;
            }
            if (lp == ma-1) {
                l = ( ll + 1 ) / ma - 1;
                nbc[nn][0] = BC_INTER_SUBDOMAINS;
                nbc[nn][1] = nbbc[mb+nn*(nblocks+1)+3*(nblocks+1)];
                mcd[nn][0] = myid - mp;
                mcd[nn][1] = mmcd[nn][1];
            }
        }
        
        switch(nn) {
            case 0:
                lxi = l;
                break;
            case 1:
            	let = l;
                break;
            case 2:
                lze = l;
                break;
        }
    }

    if (myid == 0) {
        lxim[0] = lxi;
        letm[0] = let;
        lzem[0] = lze;
        for (int i=1; i<=mpro; i++) {
        	int itag = 1;
        	safe_mpi( MPI_Recv(&lxim[i], 1, MPI_INT, i, itag, 
        		               MPI_COMM_WORLD, MPI_STATUS_IGNORE) );
        	itag = 2;
        	safe_mpi( MPI_Recv(&letm[i], 1, MPI_INT, i, itag, 
        		               MPI_COMM_WORLD, MPI_STATUS_IGNORE) );
        	itag = 3;
        	safe_mpi( MPI_Recv(&lzem[i], 1, MPI_INT, i, itag, 
        		               MPI_COMM_WORLD, MPI_STATUS_IGNORE) );
        }
    } else {
        int itag = 1;
        safe_mpi( MPI_Send(&lxi, 1, MPI_INT, 0, itag, MPI_COMM_WORLD) );
        itag = 2;
        safe_mpi( MPI_Send(&let, 1, MPI_INT, 0, itag, MPI_COMM_WORLD) );
        itag = 3;
        safe_mpi( MPI_Send(&lze, 1, MPI_INT, 0, itag, MPI_COMM_WORLD) );
    }
    safe_mpi( MPI_Bcast(lxim, mpro+1, MPI_INT, 0, MPI_COMM_WORLD) );
    safe_mpi( MPI_Bcast(letm, mpro+1, MPI_INT, 0, MPI_COMM_WORLD) );
    safe_mpi( MPI_Bcast(lzem, mpro+1, MPI_INT, 0, MPI_COMM_WORLD) );
    
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
    
    nbsize[0] = (ijk[1][0] + 1) * (ijk[2][0] + 1);
    nbsize[1] = (ijk[1][1] + 1) * (ijk[2][1] + 1);
    nbsize[2] = (ijk[1][2] + 1) * (ijk[2][2] + 1);

    printf("lxi = %d\n", lxi);
    printf("let = %d\n", let);
    printf("lze = %d\n", lze);
    

}

static int idsd3(int i, int j, int k, int mm, int nn) {
    switch(nn) {
        case 0:
            return mo[mm] + ( k * nbpc[mm+1*(nblocks+1)] + j ) *
                            nbpc[mm] + i; 
            break;
        case 1:
            return mo[mm] + ( j * nbpc[mm+1*(nblocks+1)] + i ) *
                            nbpc[mm] + k;
            break;
        case 2:
            return mo[mm] + ( i * nbpc[mm+1*(nblocks+1)] + k ) *
                            nbpc[mm] + j;
            break;
    }
}