#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

int indx3(int i, int j, int k, int nn, int lxi, int let) {
    switch(nn)
    {
        case 0:
            return (k*(let+1)+j)*(lxi+1)+i;
            break;
        case 1:
            return (j*(let+1)+i)*(lxi+1)+k;
            break;
        case 2:
            return (i*(let+1)+k)*(lxi+1)+j;
            break;
    }
}

int maxloc(double *p, int size) {
    int imax;
    double dmax;

    imax = 0;
    dmax = p[0];
    for(int i=1; i<size; i++) {
        if (p[i] > dmax) {
            dmax = p[i];
            imax = i;
        }
    }
    return imax;
}

void mtrxi(double *mtrx, double *imtrx, int size) {
	int size1d = sqrt(size);
    int ipvt[size1d];
    double rx[size];
    double arx[size];
    double sx[size];
    double temp[size1d];
    double dum;
    
    // initialize work array
    for (int i=0; i<size; i++) {
        rx[i] = mtrx[i];
    }
    // initialize index work array
    for (int i=0; i<size1d; i++)
        ipvt[i] = i;
    
    for (int i=0; i<size1d; i++) {
    	// absolute value of work array
        for (int k=0; k<size; k++) {
            arx[k] = abs(rx[k]);
        }
        int loc = i*size1d+i;
        int imax = maxloc(arx+loc, size1d-i-1);
        int mk = i + imax;
        
        // swap elements of ipvt and rx
        if (mk!=i) {
            dum = ipvt[mk];
            ipvt[mk] = ipvt[i];
            ipvt[i] = dum;
            for (int l=0; l<size1d; l++) {
            	int mmk = mk + l * size1d;
            	int ii  = i  + l * size1d;
                dum = rx[mmk];
                rx[mmk] = rx[ii];
                rx[ii] = dum;
            }
        }

        double ra0 = 1.0 / rx[i + i*size1d];
        // fill temporary array
        for (int k=0; k<size1d; k++)
            temp[k] = rx[k+i*size1d];
        
        for (int k=0; k<size1d; k++) {
            double ra1 = ra0 * rx[i+k*size1d];
            for (int l=0; l<size1d; l++)
                rx[l+k*size1d] = rx[l+k*size1d] - ra1 * temp[l];
            rx[i+k*size1d] = ra1;
        }

        for (int k=0; k<size1d; k++)
            rx[k+i*size1d] = -ra0 * temp[k];
        rx[i+i*size1d] = ra0;
    }
    
    // inverse matrix
    for (int j=0; j<size1d; j++)
        for (int i=0; i<size1d; i++)
            imtrx[i+ipvt[j]*size1d] = rx[i+j*size1d];
    
}

void matmul_square(double *mat1, double *mat2, double *rslt, int size)
{
    int size1d = sqrt(size);
    
    for (int j = 0; j < size1d; j++) {
        for (int i = 0; i < size1d; i++) {
            rslt[i+j*size1d] = 0.0;
            for (int k = 0; k < size1d; k++)
                rslt[i+j*size1d] += mat1[i+k*size1d] * mat2[k+j*size1d];
        }
    }
    
}