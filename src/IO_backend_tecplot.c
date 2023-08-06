#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "IO_backend_tecplot.h"
#include "check.h"

static int *lhmb;

void IO_write_output_file_tecplot(int nblocks,
	                              int lxio, int leto, int lzeo,
                                  int lmx, int mq, int mb, int *mo,
                                  struct t_blocks nbpc,
                                  int *lxim, int *letm, int *lzem,
                                  int *lpos,
                                  float *vart) {
    
    lhmb = safe_malloc(nblocks*sizeof(int));

    int myid;
    safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );

    int ltomb = (lxio + 1) + (leto + 1) + (lzeo + 1);

    int lje = -1;
    int ljs = lje + 1;
    lje = ljs + mq * (lmx + 1) - 1;
    int llmb = mq * ltomb - 1;

    float *vara = safe_malloc((llmb+1)*sizeof(float));
    float *varb = safe_malloc((llmb+1)*sizeof(float));
    
    if (myid == mo[mb]) {
        int mps = mo[mb];
        int mpe = mps + nbpc.coordinate[0].ptr[mb] * 
                        nbpc.coordinate[1].ptr[mb] * 
                        nbpc.coordinate[2].ptr[mb] - 1;
        int lis = 0;
        int lie = mq * (lmx + 1) + 1;
        for (int i=0; i<lie-lis; i++)
            vara[lis+i] = vart[ljs+i];
        for (int mp=mps+1; mp<=mpe; mp++) {
            lis = lie + 1;
            lie = lis + mq * (lxim[mp] + 1) * 
                             (letm[mp] + 1) *
                             (lzem[mp] + 1) - 1;
            int lmpi = lie - lis + 1;
            int itag = 1;
//            safe_mpi( MPI_Recv() );
        }
        lis = 0;
        for (int mp=mps; mp<=mpe; mp++) {
            for (int m=0; m<mq; m++) {
                for (int k=0; k<=lzem[mp]; k++) {
                    for (int j=0; j<=letm[mp]; j++) {
                        ljs = lpos[mp] + m * ltomb + k * (leto + 1) *
                                                         (lxio + 1) +
                                                     j * (lxio + 1);
                        for (int i=0; i<=lxim[mp]; i++)
                            varb[ljs+i] = vara[lis+i];
                        lis = lis + lxim[mp] + 1;
                    }
                }
            }
        }
        safe_free(vara);
//        IO_write_output_file_tecplot_block();
    } else {
        int lmpi = lje-ljs + 1;
        int itag = 1;
//        safe_mpi( MPI_Send() );
        safe_free(vara);
    }
    safe_free(varb);
    
    safe_free(lhmb);
}

void IO_write_output_file_tecplot_block(int nblocks,
                                        int lxio, int leto, int lzeo,
                                        int lmx, int mq, int mb) {
    char fhead_cname[] = "thead_block";
    char mb_str[5];
    int myid;
    safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );
    
    int ltomb = (lxio + 1) + (leto + 1) + (lzeo + 1);
    
    int lje = -1;
    int ljs = lje + 1;
    lje = ljs + mq * (lmx + 1) - 1;

    snprintf(mb_str, 5, "%d", mb);
    FILE *fhead_ptr;
    fhead_ptr = fopen(strcat(fhead_cname, mb_str), "wb");
    // call techead
    int headsize = IO_write_techead(fhead_ptr, mb, mq);
    fclose(fhead_ptr);
    unsigned char head_buffer[headsize];
    fhead_ptr = fopen(strcat(fhead_cname, mb_str), "wb");
    fread(head_buffer,sizeof(head_buffer),1,fhead_ptr);
    fclose(fhead_ptr);
    
    
    
}

int IO_write_techead(FILE * fptr, int mb, int mq) {
    int nbites = 0;
    if (mb == 0) {
        
        fprintf(fptr, "#!TDV112\n");
        nbites += 9;
        
        int head = 1;
        fwrite(&head, sizeof(int), 1 , fptr);
        nbites += sizeof(int);
        
        int ftype = 0;
        fwrite(&ftype, sizeof(int), 1, fptr);
        nbites += sizeof(int);
        
        char * filename = "n0001.plt\n";
        fwrite(filename, sizeof(char), strlen(filename), fptr);
        nbites += sizeof(char)*strlen(filename);
        
        fwrite(&mq, sizeof(int), 1, fptr);
        nbites += sizeof(int);
        
        char * xname = "x\n";
        fwrite(xname, sizeof(char), strlen(xname), fptr);
        nbites += sizeof(char)*strlen(xname);
        
        char * yname = "y\n";
        fwrite(yname, sizeof(char), strlen(yname), fptr);
        nbites += sizeof(char)*strlen(yname);
        
        char * zname = "z\n";
        fwrite(zname, sizeof(char), strlen(zname), fptr);
        nbites += sizeof(char)*strlen(zname);
        
    }
}