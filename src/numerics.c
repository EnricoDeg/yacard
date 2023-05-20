#include <math.h>
#include <stdio.h>
#include "numerics.h"
#include "utils.h"
#include "check.h"
#include "parameters.h"
#include "config_parser.h"
#include "check_mpi.h"

#define PI 3.1415926535897932384626433

struct t_blim {
    int x;
    int y;
    int z;
};

struct t_patch {
    int lmx;
    int lxi;
    int let;
    int lze;
    int lim;
    int ijk[3][3];
    int nbc[2][3];
    int mcd[2][3];
};

static struct t_blim blim;
static struct t_patch patch;
static double alphf, betf, fa, fb, fc, pfltk, pfltrbc;
static double pbci[2][2][lmp+1], pbco[2][2][lmp+1];
static double pbcot[2][2];
static double albef[2][3][5];
static double fbc[3][5];
static int ndf[2][2][3];
static double sap[lmp+1];

static double *xu, *yu, *xl, *yl;

static double *send01, *send02, *send03, *send11, *send12, *send13;
static double *recv01, *recv02, *recv03, *recv11, *recv12, *recv13;
static double *send, *recv;

static double *drva1, *drva2, *drva3;
static double *drva;

static double *sar, *sbr;

static char component[] = "numerics";

void numerics_allocate(int lim, int nbsize[3], int lmx) {
    int i = nbsize[0] - 1;
    int j = nbsize[1] - 1;
    int k = nbsize[2] - 1;

    blim.x = i;
    blim.y = j;
    blim.z = k;

    xu     = safe_malloc(3*(lim+1)*sizeof(double));
    yu     = safe_malloc(3*(lim+1)*sizeof(double));
    xl     = safe_malloc(2*(lim+1)*sizeof(double));
    yl     = safe_malloc(2*(lim+1)*sizeof(double));

    sar    = safe_malloc((lmx+1)*sizeof(double));
    sbr    = safe_malloc((lmx+1)*sizeof(double));

    send01 = safe_malloc(2*2*(i+1)*sizeof(double));
    send02 = safe_malloc(2*2*(j+1)*sizeof(double));
    send03 = safe_malloc(2*2*(k+1)*sizeof(double));

    recv01 = safe_malloc(5*2*2*(i+1)*sizeof(double));
    recv02 = safe_malloc(5*2*2*(j+1)*sizeof(double));
    recv03 = safe_malloc(5*2*2*(k+1)*sizeof(double));

    send11 = safe_malloc(2*3*(i+1)*sizeof(double));
    send12 = safe_malloc(2*3*(j+1)*sizeof(double));
    send13 = safe_malloc(2*3*(k+1)*sizeof(double));

    recv11 = safe_malloc(5*2*3*(i+1)*sizeof(double));
    recv12 = safe_malloc(5*2*3*(j+1)*sizeof(double));
    recv13 = safe_malloc(5*2*3*(k+1)*sizeof(double));

    drva1  = safe_malloc(2*5*(i+1)*sizeof(double));
    drva2  = safe_malloc(2*5*(j+1)*sizeof(double));
    drva3  = safe_malloc(2*5*(k+1)*sizeof(double));

}

void numerics_free() {

    safe_free(xu);
    safe_free(yu);
    safe_free(xl);
    safe_free(yl);

    safe_free(sar);
    safe_free(sbr);

    safe_free(send01);
    safe_free(send02);
    safe_free(send03);

    safe_free(recv01);
    safe_free(recv02);
    safe_free(recv03);

    safe_free(send11);
    safe_free(send12);
    safe_free(send13);

    safe_free(recv11);
    safe_free(recv12);
    safe_free(recv13);

    safe_free(drva1);
    safe_free(drva2);
    safe_free(drva3);

}

void numerics_read_input(char *filename) {
    config_open(filename, component);
    pfltk   = config_read_float("fltk");
    pfltrbc = config_read_float("fltrbc");
    config_close();
    pfltk *= PI;
}

void numerics_init(int lxi, int let, int lze, int nbc[2][3], int lim, 
	               int lmx, int ijk[3][3], int mcd[2][3]) {
	int nstart, nend, istart, iend;

	init_coeff();
    
    for (int nn=0; nn<3; nn++) {
        switch(nn) {
            case 0:
                istart = 0;
                iend   = istart + lxi;
                break;
            case 1:
            	istart = lxi + 1;
                iend   = istart + let;
                break;
            case 2:
                istart = lxi + let + 2;
                iend   = istart + lze;
                break;
        }

        for (int ip=0; ip<2; ip++) {
            int np = nbc[ip][nn];
            switch(np) {
                case BC_NON_REFLECTIVE:
                case BC_WALL_INVISCID:
                case BC_WALL_VISCOUS:
                case BC_INTER_CURV:
                    ndf[0][ip][nn] = 0;
                    ndf[1][ip][nn] = 0;
                    break;
                case BC_INTER_STRAIGHT:
                case BC_INTER_SUBDOMAINS:
                case BC_PERIODIC:
                    ndf[0][ip][nn] = 1;
                    ndf[1][ip][nn] = 1;
                	break;
            }
        }
        
        nstart = ndf[0][0][nn];
        nend   = ndf[0][1][nn];
        penta(xu, xl, istart, iend, nstart, nend, 0, lim);

        nstart = ndf[1][0][nn];
        nend   = ndf[1][1][nn];
        penta(yu, yl, istart, iend, nstart, nend, 1, lim);
    }

    // fill domain minimal decomposition information for internal use
    patch.lmx = lmx;
    patch.lxi = lxi;
    patch.let = let;
    patch.lze = lze;
    patch.lim = lim;
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            patch.ijk[i][j] = ijk[i][j];
    for (int nn=0; nn<3; nn++)
        for (int ip=0; ip<2; ip++)
            patch.nbc[ip][nn] = nbc[ip][nn];
    for (int nn=0; nn<3; nn++)
        for (int ip=0; ip<2; ip++)
            patch.mcd[ip][nn] = mcd[ip][nn];

}

void numerics_deriv(double *rfield, int nn, int nz, int m) {
    
    int ustart, uend;
    int ntk    = 0;
    int nstart = ndf[0][0][nn];
    int nend   = ndf[0][1][nn];
    int slim;
    int lmx = patch.lmx;
    int lxi = patch.lxi;
    int let = patch.let;
    int lze = patch.lze;
    int lim = patch.lim;
    int ijk[3][3];
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            ijk[i][j] = patch.ijk[i][j];
    int klim   = ijk[2][nn];
    int jlim   = ijk[1][nn];
    int ilim   = ijk[0][nn];
    
    switch(nn) {
        case 0:
            ustart = 0;
            uend   = ustart + lxi;
            recv   = recv01;
            drva   = drva1;
            slim   = blim.x;
            break;
        case 1:
            ustart = lxi + 1;
            uend   = ustart + let;
            recv   = recv02;
            drva   = drva2;
            slim   = blim.y;
            break;
        case 2:
            ustart = lxi + let + 2;
            uend   = ustart + lze;
            recv   = recv03;
            drva   = drva3;
            slim   = blim.z;
            break;
    }
    
    // fill sar with input field
    for (int k=0; k<=klim; k++) {
        for (int j=0; j<=jlim; j++) {
            for (int i=0; i<=ilim; i++) {
                int l = indx3(i, j, k, nn, lxi, let);
                int ij = i + 
                         j * (ilim+1) +
                         k * ((ilim+1)*(jlim+1));
                sar[ij] = rfield[l+nz*(lmx+1)];
            }
    	}
    }
    
    // fill sbr lower boundary
    if (nstart == 0) {
        for (int k=0; k<=klim; k++) {
            for (int j=0; j<=jlim; j++) {
                int istart = j * (ilim+1) +
                                k * ((ilim+1)*(jlim+1));
                int iend   = istart + ilim;
                int kp = k * (jlim + 1);
                int jk = kp + j;
                
                {
                    sbr[istart] = 0.0;
                    double tmp[4] = {a01, a02, a03, a04};
                    int idx[4] = {1, 2, 3, 4};
                    for (int t=0; t<4; t++)
                        sbr[istart] += ( tmp[t] * (sar[istart+idx[t]] - sar[istart]) );
                }

                {
                    sbr[istart+1] = 0.0;
                    double tmp[4] = {a10, a12, a13, a14};
                    int idx[4] = {0, 2, 3, 4};
                    for (int t=0; t<4; t++)
                        sbr[istart+1] += (tmp[t] * (sar[istart+idx[t]] - sar[istart+1]) );
                }
            }
        }
    } else if (nstart == 1) {
        for (int k=0; k<=klim; k++) {
            for (int j=0; j<=jlim; j++) {
                int istart = j * (ilim+1) +
                                k * ((ilim+1)*(jlim+1));
                int iend   = istart + ilim;
                int kp = k * (jlim + 1);
                int jk = kp + j;
                
                {
                    sbr[istart] = 0.0;
                    for (int t=0; t<lmp+1; t++)
                        sbr[istart] += (pbci[ntk][0][t] * sar[istart+t]);
                    int f4 = jk + 0*slim + 0*(slim*2) + m*(slim*2*2);
                    sbr[istart] += recv[f4];
                }
                
                {
                    sbr[istart+1] = 0.0;
                    for (int t=0; t<lmp+1; t++)
                        sbr[istart+1] += (pbci[ntk][1][t] * sar[istart+t]);
                    int f4 = jk + 1*slim + 0*(slim*2) + m*(slim*2*2);
                    sbr[istart+1] += recv[f4];
                }
            }
        }
    }
    
    // fill sbr inner
    for (int k=0; k<=klim; k++) {
        for (int j=0; j<=jlim; j++) {
            int istart = j * (ilim+1) +
                            k * ((ilim+1)*(jlim+1));
            int iend   = istart + ilim;
            for (int i=istart+2; i<=iend-2; i++) {
                sbr[i] = aa * ( sar[i+1] - sar[i-1] ) + 
                         ab * ( sar[i+2] - sar[i-2] );
            }
        }
    }
    
    // fill sbr upper boundary
    if (nend == 0) {
        for (int k=0; k<=klim; k++) {
            for (int j=0; j<=jlim; j++) {
                int istart = j * (ilim+1) +
                                k * ((ilim+1)*(jlim+1));
                int iend   = istart + ilim;
                int kp = k * (jlim + 1);
                int jk = kp + j;
                
                {
                    sbr[iend] = 0.0;
                    double tmp[4] = {a01, a02, a03, a04};
                    int idx[4] = {1, 2, 3, 4};
                    for (int t=0; t<4; t++)
                        sbr[iend] += ( tmp[t] * (sar[iend] - sar[iend-idx[t]]) );
                }
                
                {
                    sbr[iend-1] = 0.0;
                    double tmp[4] = {a10, a12, a13, a14};
                    int idx[4] = {0, 2, 3, 4};
                    for (int t=0; t<4; t++)
                        sbr[iend+1] += ( tmp[t] * (sar[iend-1] - sar[iend-idx[t]]) );
                }
            }
        }
    } else if (nend == 1) {
        for (int k=0; k<=klim; k++) {
            for (int j=0; j<=jlim; j++) {
                int istart = j * (ilim+1) +
                                k * ((ilim+1)*(jlim+1));
                int iend   = istart + ilim;
                int kp = k * (jlim + 1);
                int jk = kp + j;
                
                {
                    sbr[iend] = 0.0;
                    for (int t=0; t<lmp+1; t++)
                        sbr[iend] -= (pbci[ntk][0][t] * sar[iend-t]);
                    int f4 = jk + 0*slim + 1*(slim*2) + m*(slim*2*2);
                    sbr[iend] -= recv[f4];
                }
                
                {
                    sbr[iend-1] = 0.0;
                    for (int t=0; t<lmp+1; t++)                    	
                        sbr[iend-1] -= (pbci[ntk][1][t] * sar[iend-t]);
                    int f4 = jk + 1*slim + 1*(slim*2) + m*(slim*2*2);
                    sbr[iend-1] -= recv[f4];
                }
            }
        }
    }
    
    // fill sar
    for (int k=0; k<=klim; k++) {
        for (int j=0; j<=jlim; j++) {
            int istart = j * (ilim+1) +
                            k * ((ilim+1)*(jlim+1)); 
            int iend   = istart + ilim;
            sar[istart]   = sbr[istart];
            sar[istart+1] = sbr[istart+1] - xl[ustart+1+(lim+1)] * sar[istart];
            for (int i=istart+2; i<=iend; i++) {
                sar[i] = sbr[i] - xl[i-istart+ustart] * sar[i-2] 
                                - xl[i-istart+ustart+(lim+1)] * sar[i-1];
            }
        }
    }
    
    // fill sbr and drva
    for (int k=0; k<=klim; k++) {
        for (int j=0; j<=jlim; j++) {
            
            int istart = j * (ilim+1) + 
                            k * ((ilim+1)*(jlim+1)); 
            int iend   = istart + ilim;
            
            int kpp = k * (jlim+1);
            int jkk = kpp + j;

            sbr[iend]   = xu[uend] * sar[iend];
            sbr[iend-1] = xu[uend-1] * sar[iend-1] - 
                          xu[uend-1+(lim+1)] * sbr[iend];            

            for (int i=iend-2; i>istart-1; i--) {
                sbr[i] = xu[i-iend+uend]           * sar[i]   - 
                         xu[i-iend+uend+(lim+1)]   * sbr[i+1] - 
                         xu[i-iend+uend+2*(lim+1)] * sbr[i+2];
            }
            
            drva[jkk+m*(slim+1)]            = sbr[istart]; 
            drva[jkk+m*(slim+1)+5*(slim+1)] = sbr[iend];
            
        }
    }
    
    // fill rfield
    for (int k=0; k<=klim; k++) {
        for (int j=0; j<=jlim; j++) {
            for (int i=0; i<=ilim; i++) {
                int l = indx3(i, j, k, nn, lxi, let);
                int ij = i + 
                         j * (ilim+1) +
                         k * ((ilim+1)*(jlim+1));
                rfield[l+nn*(lmx+1)] = sbr[ij];
            }
    	}
    }
    
}

void numerics_halo_exch(double *rfield, int nt, int nrt, int n45, int m) {

    MPI_Request ireq[MAX_MPI_REQUESTS];
    MPI_Status  ista[MAX_MPI_REQUESTS];
    int ir = 0;
    int slim;
    int dim2;
    int lmx = patch.lmx;
    int lxi = patch.lxi;
    int let = patch.let;
    int ijk[3][3];
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            ijk[i][j] = patch.ijk[i][j];
    int nbc[2][3];
    for (int nn=0; nn<3; nn++)
        for (int ip=0; ip<2; ip++)
            nbc[ip][nn] = patch.nbc[ip][nn];
    int mcd[2][3];
    for (int nn=0; nn<3; nn++)
        for (int ip=0; ip<2; ip++)
            mcd[ip][nn] = patch.mcd[ip][nn];
    
    for (int nn=0; nn<3; nn++) {
    	int nz = ( 1 - nrt ) * nn;
        if (nt == 0) {
        	dim2 = 2;
            switch(nn) {
                case 0:
                    send = send01;
                    recv = recv01;
                    slim = blim.x;
                    break;
                case 1:
            	    send = send02;
                    recv = recv02;
                    slim = blim.y;
                    break;
                case 2:
                    send = send03;
                    recv = recv03;
                    slim = blim.z;
                    break;
            }
        } else {
        	dim2 = 3;
            switch(nn) {
                case 0:
                    send = send11;
                    recv = recv11;
                    slim = blim.x;
                    break;
                case 1:
            	    send = send12;
            	    recv = recv12;
            	    slim = blim.y;
            	    break;
                case 2:
                    send = send13;
                    recv = recv13;
                    slim = blim.z;
                    break;
            }
        }
        
        for (int ip=0; ip<2; ip++) {
            int iq     = 1 - ip;
            int istart = ip * ijk[0][nn];
            int iend   = 1 - 2 * ip;

            double ra0;
            int ii;
            switch(nbc[ip][nn]) {
                case(BC_INTER_STRAIGHT):
                    ra0 = 0.0;
                    ii = 1;
                    break;
                case(BC_INTER_SUBDOMAINS):
                    ra0 = 0.0;
                    ii = 0;
                	break;
                case(BC_PERIODIC):
                    ra0 = (double)n45;
                    ii = 1;
                	break;
            }
            
            if (ndf[nt][ip][nn]) {
            	// pack buffer
                for (int k=0; k<=ijk[2][nn]; k++) {
                    int kpp = k * (ijk[1][nn] + 1);
                    for (int j=0; j<=ijk[1][nn]; j++) {
                        int jkk = kpp + j;
                        int lll = indx3(istart, j, k, nn, lxi, let);
                        double res = ra0 * rfield[lll+nz*(lmx+1)];
                        for (int i=0; i<=lmp; i++) {
                            lll = indx3(istart+iend*(i+ii), j, k, nn, lxi, let);
                            sap[i] = rfield[lll+nz*(lmx+1)];
                        }

                        int f4 = jkk + 0*slim + ip*(slim*dim2);
                        send[f4] = 0.0;
                        for (int f=0; f<=lmp; f++)
                            send[f4] += (pbco[nt][0][f] * sap[f]);
                        send[f4] -= (res * pbcot[nt][0]);
                        
                        f4 = jkk + 1*slim + ip*(slim*dim2);
                        send[f4] = 0.0;
                        for (int f=0; f<=lmp; f++)
                            send[f4] += (pbco[nt][1][f] * sap[f]);
                        send[f4] -= (res * pbcot[nt][1]);
                        
                        f4 = jkk + (nt+1)*slim + ip*(slim*dim2);
                        send[f4] = send[f4] + nt * (sap[0] - res - send[f4]);

                    }
                }

                // halo exchange
                safe_mpi( MPI_Isend(send, slim*dim2, MPI_DOUBLE, mcd[ip][nn],
                                         iq, MPI_COMM_WORLD, &ireq[ir]) );
                ir++,
                safe_mpi( MPI_Irecv(recv, slim*dim2, MPI_DOUBLE, mcd[ip][nn],
                                         ip, MPI_COMM_WORLD, &ireq[ir]) );
                ir++;
                
            }

        }  
    }
    if (ir > 0)
        safe_mpi( MPI_Waitall(ir+1, ireq, ista) );

    // periodic boundary conditions
    if (n45 == N45go) {

        for (int nn=0; nn<3; nn++) {
            int nz = (1 - nrt) * nn;
            if (nt == 0) {
                dim2 = 2;
                switch(nn) {
                    case 0:
                        recv = recv01;
                        slim = blim.x;
                        break;
                    case 1:
                        recv = recv02;
                        slim = blim.y;
                        break;
                    case 2:
                        recv = recv03;
                        slim = blim.z;
                        break;
                }
            } else {
                dim2 = 3;
                switch(nn) {
                    case 0:
                        recv = recv11;
                        slim = blim.x;
                        break;
                    case 1:
            	        recv = recv12;
            	        slim = blim.y;
            	        break;
                    case 2:
                        recv = recv13;
                        slim = blim.z;
                        break;
                }
            }
            
            for (int ip=0; ip<2; ip++) {
                int istart = ip*ijk[0][nn];
                if (nbc[ip][nn] == BC_PERIODIC) {
                    for (int k=0; k<=ijk[2][nn]; k++) {
                        int kpp = k * (ijk[1][nn] + 1);
                        for (int j=0; j<=ijk[1][nn]; j++) {
                            int jkk = kpp + j;
                            int lll = indx3(istart, j, k, nn, lxi, let);

                            int f4 = jkk + 0*slim + ip*(slim*dim2) + m*(2*slim*dim2);
                            recv[f4] += (rfield[lll+nz*(lmx+1)] * pbcot[nt][0]);

                            f4 = jkk * 1*slim + ip*(slim*dim2) + m*(2*slim*dim2);
                            recv[f4] += (rfield[lll+nz*(lmx+1)] * pbcot[nt][1]);

                            f4 = jkk + (nt+1)*slim + ip*(slim*dim2) + m*(2*slim*dim2);
                            recv[f4] += (nt * rfield[lll+nz*(lmx+1)]);
                        }
                    }
                }
            }
            
        }
    }

}

static void init_coeff() {

    fcbcm(pfltk, pfltrbc);
    fcint(pfltk, 0.5, &alphf, &betf, &fa, &fb, &fc);

    albef[1][0][0] = 0.0;
    albef[1][0][1] = 0.0;
    albef[1][0][2] = 1.0;
    albef[1][0][3] = alphf;
    albef[1][0][4] = betf;

    albef[1][1][0] = 0.0;
    albef[1][1][1] = alphf;
    albef[1][1][2] = 1.0;
    albef[1][1][3] = alphf;
    albef[1][1][4] = betf;

    albef[1][2][0] = betf;
    albef[1][2][1] = alphf;
    albef[1][2][2] = 1.0;
    albef[1][2][3] = alphf;
    albef[1][2][4] = betf;

    for (int k=0; k<2; k++) {
        for (int j=0; j<2; j++) {
        	for (int i=0; i<lmp+1; i++) {
        		pbco[k][j][i] = 0.0;
        		pbci[k][j][i] = 0.0;
        	}
        }
    }

    sbcco();

    for (int k=0; k<2; k++) {
        for (int j=0; j<2; j++) {
            int dsum = 0.0;
            for (int i=0; i<lmp+1; i++) {
                dsum += pbco[k][j][i];
            }
            pbcot[k][j] = dsum;
        }
    }

}

static void fcint(double fltk, double fltr, double *alphz,
	              double *betz, double *za, double *zb, double *zc) {
    double cosf[3], fctrk;

    cosf[0] = cos(fltk  );
    cosf[1] = cos(2.0*fltk);
    cosf[2] = cos(3.0*fltk);

    fctrk = 1.0 / ( 30.0 + 5.0 * (7.0-16.0*fltr) * cosf[0] + 
            2.0 * (1.0+8.0*fltr) * cosf[1] - 3.0 * cosf[2]);
    *alphz = fctrk * (20.0 * (2.0*fltr-1) - 30.0 * cosf[0] + 
            12.0 * (2.0 * fltr - 1.0) * cosf[1] - 2.0 * cosf[2]);
    *betz  = 0.5 * fctrk * (2.0 * (13.0 - 8.0 * fltr) + 
           (33.0 - 48.0 * fltr) * cosf[0] + 6.0 * cosf[1] - cosf[2]);
    *za    = 60.0 * (1.0-fltr) * pow(cos(0.5*fltk),4.0) * fctrk;
    *zb    = -2.0 * (*za) / 5.0;
    *zc    = (*za) / 15.0;

}

static void fcbcm(double fltk, double fltrbc) {
    double alphz, betz, za, zb, zc, aok, fctrk, resk;

    aok = log(fltrbc);
    fcint(fltk, 0.5, &alphz, &betz, &za, &zb, &zc);
    fctrk = 1.0 / ( 1.0 + alphz * fltrbc + betz * pow(fltrbc,2.0));
    
    albef[0][0][0] = 0.0;
    albef[0][0][1] = 0.0;
    albef[0][0][2] = 1.0;
    albef[0][0][3] = alphz * fctrk;
    albef[0][0][4] = betz *  fctrk;

    resk = (fltrbc-1.0) * (za + zc + (fltrbc+1) * (zb+fltrbc*zc)) / aok;
    fbc[0][0] = (za - 5.0 * resk / 3.0) * fctrk;
    fbc[0][1] = (zb + 10.0 * resk / 21.0) * fctrk;
    fbc[0][2] = (zc - 5.0 * resk / 42.0) * fctrk;
    fbc[0][3] = (5.0 * resk / 252.0) * fctrk;
    fbc[0][4] = (-resk / 630.0) * fctrk;
    
    albef[0][1][0] = 0.0;
    albef[0][1][1] = alphz + betz * fltrbc;
    albef[0][1][2] = 1.0;
    albef[0][1][3] = alphz;
    albef[0][1][4] = betz;
    
    resk = (fltrbc-1.0) * (zb + zc * (fltrbc+1.0)) / aok;
    fbc[1][0] = za + zb + zc + 1627.0 * resk / 1260.0;
    fbc[1][1] = za + 10.0 * resk / 21.0;
    fbc[1][2] = zb - 5.0 * resk / 42.0;
    fbc[1][3] = zc + 5.0 * resk / 252.0;
    fbc[1][4] = -resk / 630.0;

    albef[0][2][0] = betz;
    albef[0][2][1] = alphz;
    albef[0][2][2] = 1.0;
    albef[0][2][3] = alphz;
    albef[0][2][4] = betz;
    
    resk = zc * (fltrbc-1.0) / aok;
    fbc[2][0] = zb + zc + 1627.0 * resk / 1260.0;
    fbc[2][1] = za - 5.0 * resk / 3.0;
    fbc[2][2] = za - 5.0 * resk / 42.0;
    fbc[2][3] = zb + 5.0 * resk / 252.0;
    fbc[2][4] = zc - resk / 630.0;

}

static void sbcco() {
    double *ax, *bx, *rx, *sx, *tmp;
    int l, istart, iend;
    
    for(int nt=0; nt<2; nt++) {
        l = lmp;
        iend = 2 * (l + 1);
        ax  = safe_malloc((iend)*(iend)*sizeof(double));
        bx  = safe_malloc((iend)*(iend)*sizeof(double));
        sx  = safe_malloc((iend)*(iend)*sizeof(double));
        rx  = safe_malloc((iend)*(iend)*sizeof(double));
        tmp = safe_malloc((iend)*(iend)*sizeof(double));

        for (int j=0; j<iend; j++) {
            for (int i=0; i<iend; i++) {
                ax[i+j*iend] = 0.0;
                bx[i+j*iend] = 0.0;
            }
        }

        int istart = 0;

        if (nt == 0) {
            
            ax[istart+istart*iend]     = 1.0;
            ax[istart+(istart+1)*iend] = alpha01;
            ax[istart+(istart+2)*iend] = beta02;

            bx[istart+istart*iend]     = - (a01 + a02 + a03 + a04);
            bx[istart+(istart+1)*iend] = a01;
            bx[istart+(istart+2)*iend] = a02;
            bx[istart+(istart+3)*iend] = a03;
            bx[istart+(istart+4)*iend] = a04;

            ax[(istart+1)+istart*iend]     = alpha10;
            ax[(istart+1)+(istart+1)*iend] = 1.0;
            ax[(istart+1)+(istart+2)*iend] = alpha12;
            ax[(istart+1)+(istart+3)*iend] = beta13;

            bx[(istart+1)+istart*iend]     = a10;
            bx[(istart+1)+(istart+1)*iend] = -(a10 + a12 + a13 + a14);
            bx[(istart+1)+(istart+2)*iend] = a12;
            bx[(istart+1)+(istart+3)*iend] = a13;
            bx[(istart+1)+(istart+4)*iend] = a14;

            for (int i=istart+2; i<iend-2; i++) {
                ax[i+(i-2)*iend] = beta;
                ax[i+(i-1)*iend] = alpha;
                ax[i+i*iend]     = 1.0;
                ax[i+(i+1)*iend] = alpha;
                ax[i+(i+2)*iend] = beta;

                bx[i+(i-2)*iend] = -ab;
                bx[i+(i-1)*iend] = -aa;
                bx[i+i*iend]     = 0.0;
                bx[i+(i+1)*iend] = aa;
                bx[i+(i+2)*iend] = ab;
            }
            
            ax[(iend-2)+(iend-1)*iend] = ax[(istart+1)+istart*iend];
            ax[(iend-2)+(iend-2)*iend] = ax[(istart+1)+(istart+1)*iend];
            ax[(iend-2)+(iend-3)*iend] = ax[(istart+1)+(istart+2)*iend];
            ax[(iend-2)+(iend-4)*iend] = ax[(istart+1)+(istart+3)*iend];
            
            bx[(iend-2)+(iend-1)*iend] = -bx[(istart+1)+istart*iend];
            bx[(iend-2)+(iend-2)*iend] = -bx[(istart+1)+(istart+1)*iend];
            bx[(iend-2)+(iend-3)*iend] = -bx[(istart+1)+(istart+2)*iend];
            bx[(iend-2)+(iend-4)*iend] = -bx[(istart+1)+(istart+3)*iend];
            bx[(iend-2)+(iend-5)*iend] = -bx[(istart+1)+(istart+4)*iend];
            
            ax[(iend-1)+(iend-1)*iend] = ax[istart+istart*iend];
            ax[(iend-1)+(iend-2)*iend] = ax[istart+(istart+1)*iend];
            ax[(iend-1)+(iend-3)*iend] = ax[istart+(istart+2)*iend];
            
            bx[(iend-1)+(iend-1)*iend] = -bx[istart+istart*iend];
            bx[(iend-1)+(iend-2)*iend] = -bx[istart+(istart+1)*iend];
            bx[(iend-1)+(iend-3)*iend] = -bx[istart+(istart+2)*iend];
            bx[(iend-1)+(iend-4)*iend] = -bx[istart+(istart+3)*iend];
            bx[(iend-1)+(iend-5)*iend] = -bx[istart+(istart+4)*iend];

        }
        
        if (nt == 1) {

            ax[istart+istart*iend]     = albef[0][0][2];
            ax[istart+(istart+1)*iend] = albef[0][0][3];
            ax[istart+(istart+2)*iend] = albef[0][0][4];

            for (int f=0; f<5; f++)
                bx[istart+(istart+1+f)*iend] = fbc[0][f];

            bx[istart+istart*iend] = -fbc[0][0];
            for (int f=1; f<5; f++)
                bx[istart+istart*iend] -= fbc[0][f];

            ax[(istart+1)+istart*iend]     = albef[0][1][1];
            ax[(istart+1)+(istart+1)*iend] = albef[0][1][2];
            ax[(istart+1)+(istart+2)*iend] = albef[0][1][3];
            ax[(istart+1)+(istart+3)*iend] = albef[0][1][4];

            {
                int idx[5] = {0,2,3,4,5};
                for (int f=0; f<5; f++)
                    bx[(istart+1)+(istart+idx[f])*iend] = fbc[1][f];
            }

            bx[(istart+1)+(istart+1)*iend] = -fbc[1][0];
            for (int f=1; f<5; f++)
                bx[(istart+1)+(istart+1)*iend] -= fbc[1][f];

            ax[(istart+2)+istart*iend]     = albef[0][2][0];
            ax[(istart+2)+(istart+1)*iend] = albef[0][2][1];
            ax[(istart+2)+(istart+2)*iend] = albef[0][2][2];
            ax[(istart+2)+(istart+3)*iend] = albef[0][2][3];
            ax[(istart+2)+(istart+4)*iend] = albef[0][2][4];

            {
                int idx[5] = {0,1,3,4,5};
                for (int f=0; f<5; f++)
                    bx[(istart+2)+(istart+idx[f])*iend] = fbc[2][f];
            }

            bx[(istart+2)+(istart+2)*iend] = -fbc[2][0];
            for (int f=1; f<5; f++)
                bx[(istart+2)+(istart+2)*iend] -= fbc[2][f];

            for (int i=istart+3; i<iend-3; i++) {
                ax[i+(i-2)*iend] = betf;
                ax[i+(i-1)*iend] = alphf;
                ax[i+i*iend]     = 1.0;
                ax[i+(i+1)*iend] = alphf;
                ax[i+(i+2)*iend] = betf;

                bx[i+(i-3)*iend] = fc;
                bx[i+(i-2)*iend] = fb;
                bx[i+(i-1)*iend] = fa;
                bx[i+i*iend]     = -2 * (fa + fb + fc);
                bx[i+(i+1)*iend] = fa;
                bx[i+(i+2)*iend] = fb;
                bx[i+(i+3)*iend] = fc;
            }

            ax[(iend-3)+(iend-1)*iend] = ax[(istart+2)+istart*iend];
            ax[(iend-3)+(iend-2)*iend] = ax[(istart+2)+(istart+1)*iend];
            ax[(iend-3)+(iend-3)*iend] = ax[(istart+2)+(istart+2)*iend];
            ax[(iend-3)+(iend-4)*iend] = ax[(istart+2)+(istart+3)*iend];
            ax[(iend-3)+(iend-5)*iend] = ax[(istart+2)+(istart+4)*iend];

            bx[(iend-3)+(iend-1)*iend] = bx[(istart+2)+istart*iend];
            bx[(iend-3)+(iend-2)*iend] = bx[(istart+2)+(istart+1)*iend];
            bx[(iend-3)+(iend-3)*iend] = bx[(istart+2)+(istart+2)*iend];
            bx[(iend-3)+(iend-4)*iend] = bx[(istart+2)+(istart+3)*iend];
            bx[(iend-3)+(iend-5)*iend] = bx[(istart+2)+(istart+4)*iend];
            bx[(iend-3)+(iend-6)*iend] = bx[(istart+2)+(istart+5)*iend];

            ax[(iend-2)+(iend-1)*iend] = ax[(istart+1)+istart*iend];
            ax[(iend-2)+(iend-2)*iend] = ax[(istart+1)+(istart+1)*iend];
            ax[(iend-2)+(iend-3)*iend] = ax[(istart+1)+(istart+2)*iend];
            ax[(iend-2)+(iend-4)*iend] = ax[(istart+1)+(istart+3)*iend];

            bx[(iend-2)+(iend-1)*iend] = bx[(istart+1)+istart*iend];
            bx[(iend-2)+(iend-2)*iend] = bx[(istart+1)+(istart+1)*iend];
            bx[(iend-2)+(iend-3)*iend] = bx[(istart+1)+(istart+2)*iend];
            bx[(iend-2)+(iend-4)*iend] = bx[(istart+1)+(istart+3)*iend];
            bx[(iend-2)+(iend-5)*iend] = bx[(istart+1)+(istart+4)*iend];
            bx[(iend-2)+(iend-6)*iend] = bx[(istart+1)+(istart+5)*iend];

            ax[(iend-1)+(iend-1)*iend] = ax[istart+istart*iend];
            ax[(iend-1)+(iend-2)*iend] = ax[istart+(istart+1)*iend];
            ax[(iend-1)+(iend-3)*iend] = ax[istart+(istart+2)*iend];

            bx[(iend-1)+(iend-1)*iend] = bx[istart+istart*iend];
            bx[(iend-1)+(iend-2)*iend] = bx[istart+(istart+1)*iend];
            bx[(iend-1)+(iend-3)*iend] = bx[istart+(istart+2)*iend];
            bx[(iend-1)+(iend-4)*iend] = bx[istart+(istart+3)*iend];
            bx[(iend-1)+(iend-5)*iend] = bx[istart+(istart+4)*iend];
            bx[(iend-1)+(iend-6)*iend] = bx[istart+(istart+5)*iend];

        }

        mtrxi(ax, sx, iend*iend);

        for (int j=0; j<iend; j++) {
            for (int i=0; i<iend; i++) {
                rx[i+j*iend] = ax[i+j*iend];
            }
        }

        {
            int i = iend / 2 - 2;
            rx[i+(i+2)*iend] = 0.0;
            rx[(i+1)+(i+2)*iend] = 0.0;
            rx[(i+1)+(i+3)*iend] = 0.0;
        }

        {
            int i = iend / 2 + 1;
            rx[i+(i-2)*iend] = 0.0;
            rx[(i-1)+(i-2)*iend] = 0.0;
            rx[(i-1)+(i-3)*iend] = 0.0;
        }

        matmul_square(sx, bx, tmp, iend*iend);
        matmul_square(rx, tmp, ax, iend*iend);

        {
            int i = iend / 2;
            for (int lp=0; lp<lmp+1; lp++)
            	pbco[nt][0][lmp-lp] = ax[i+(istart+lp)*iend];
            for (int lp=0; lp<lmp+1; lp++)
                pbci[nt][0][lp]       = ax[i+(istart+lmp+1+lp)*iend];
        }
        
        {
            int i = iend / 2 + 1;
            for (int lp=0; lp<lmp+1; lp++)
                pbco[nt][1][lmp-lp] = ax[i+(istart+lp)*iend];
            for (int lp=0; lp<lmp+1; lp++)
                pbci[nt][1][lp] = ax[i+(istart+lmp+1+lp)*iend];
        }
        
        safe_free(ax);
        safe_free(bx);
        safe_free(rx);
        safe_free(sx);
        safe_free(tmp);
    }

}

// CHOLESKY DECOMPOSITION OF PENTADIAGONAL MATRICES
static void penta(double *xu, double *xl, 
	              int is, int ie, int ns, int ne, int nt, int lim) {
    double albe[2][3][5];
    double alpho, beto;
    int iik;

	if (nt==0) {

        double tmp10[5] = {0.0 , 0.0    , 1.0, alpha01, beta02};
        double tmp20[5] = {0.0 , alpha10, 1.0, alpha12, beta13};
        double tmp30[5] = {beta, alpha  , 1.0, alpha  ,   beta};
        double tmp11[5] = {0.0 , 0.0    , 1.0, alpha  ,   beta};
        double tmp21[5] = {0.0 , alpha  , 1.0, alpha  ,   beta};
        double tmp31[5] = {beta, alpha  , 1.0, alpha  ,   beta};
        alpho = alpha;
        beto   = beta;
        for (int i=0; i<5; i++) {
            albe[0][0][i] = tmp10[i];
            albe[0][1][i] = tmp20[i];
            albe[0][2][i] = tmp30[i];
            albe[1][0][i] = tmp11[i];
            albe[1][1][i] = tmp21[i];
            albe[1][2][i] = tmp31[i];
        }

	} else {

        for (int i=0; i<2; i++) {
            for (int j=0; j<3; j++) {
                for (int k=0; k<5; k++) {
                    albe[i][j][k] = albef[i][j][k];
                }
            }
        }
        alpho = alphf;
        beto  = betf;

	}
    
    for (int j=0; j<2; j++) {
        for (int i=is; i<=ie; i++) {
            xl[i+j*(lim+1)] = 1.0;
        }
    }

    for (int j=0; j<3; j++) {
        for (int i=is; i<=ie; i++) {
            xu[i+j*(lim+1)] = 1.0;
        }
    }

    iik = is;
    xu[iik]           = 1.0;
    xu[iik+(lim+1)]   = albe[ns][0][3];
    xu[iik+2*(lim+1)] = albe[ns][0][4];

    iik = is + 1;
    xl[iik+(lim+1)]   = albe[ns][1][1] * xu[(iik-1)];
    xu[iik]           = 1.0 / (1.0 - xu[(iik-1)+(lim+1)] * xl[iik+(lim+1)]);
    xu[iik+(lim+1)]   = albe[ns][1][3] - xu[(iik-1)+2*(lim+1)] * xl[iik+(lim+1)];
    xu[iik+2*(lim+1)] = albe[ns][1][4];

    iik = is + 2;
    xl[iik]           = albe[ns][2][0] * xu[iik-2];
    xl[iik+(lim+1)]   = (albe[ns][2][1] - xu[(iik-2)+(lim+1)] * xl[iik]) * xu[iik-1];
    xu[iik]           = 1.0 / (1.0 - xu[(iik-2)+2*(lim+1)] * xl[iik] - xu[(iik-1)+(lim+1)] * xl[iik+(lim+1)]);
    xu[iik+(lim+1)]   = albe[ns][2][3] - xu[(iik-1)+2*(lim+1)] * xl[iik+(lim+1)];
    xu[iik+2*(lim+1)] = albe[ns][2][4];

    for (int i=is+3; i<=ie-3; i++) {
        xl[i]           = beto * xu[i-2];
        xl[i+(lim+1)]   = ( alpho - xu[(i-2)+(lim+1)] * xl[i] ) * xu[i-1];
        xu[i]           = 1.0 / (1.0 - xu[(i-2)+2*(lim+1)] * xl[i] - xu[(i-1)+(lim+1)] * xl[i+(lim+1)]);
        xu[i+(lim+1)]   = alpho - xu[(i-1)+2*(lim+1)] * xl[i+(lim+1)];
        xu[i+2*(lim+1)] = beto;
    }

    iik = ie - 2;
    xl[iik] = albe[ne][2][4] * xu[iik-2];
    xl[iik+(lim+1)] = (albe[ne][2][3] - xu[(iik-2)+(lim+1)] * xl[iik]) * xu[iik-1];
    xu[iik] = 1.0 / (1.0 - xu[(iik-2)+2*(lim+1)] * xl[iik] - xu[(iik-1)+(lim+1)] * xl[iik+(lim+1)]);
    xu[iik+(lim+1)] = albe[ne][2][1] - xu[(iik-1)+2*(lim+1)] * xl[iik+(lim+1)];
    xu[iik+2*(lim+1)] = albe[ne][2][0];

    iik = ie - 1;
    xl[iik]         = albe[ne][1][4] * xu[iik-2];
    xl[iik+(lim+1)] = (albe[ne][1][3] - xu[(iik-2)+(lim+1)] * xl[iik]) * xu[iik-1];
    xu[iik]         = 1.0 / (1.0 - xu[(iik-2)+2*(lim+1)] * xl[iik] - xu[(iik-1)+(lim+1)] * xl[iik+(lim+1)] );
    xu[iik+(lim+1)] = albe[ne][1][1] - xu[(iik-1)+2*(lim+1)] * xl[iik+(lim+1)];

    iik = ie;
    xl[iik]         = albe[ne][0][4] * xu[iik-2];
    xl[iik+(lim+1)] = ( albe[ne][0][3] - xu[(iik-2)+(lim+1)] * xl[iik] ) * xu[iik-1];
    xu[iik]         = 1.0 / (1.0 - xu[(iik-2)+2*(lim+1)] * xl[iik] - xu[(iik-1)+(lim+1)] * xl[iik+(lim+1)]);

    for (int i=is; i<=ie; i++) {
        xu[i+(lim+1)]   *= xu[i];
        xu[i+2*(lim+1)] *= xu[i];
    }

}