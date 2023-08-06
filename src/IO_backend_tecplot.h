#include "data_struct.h"

void IO_allocate(int nblocks);

void IO_free();

void IO_write_output_file_tecplot(int nblocks,
	                              int lxio, int leto, int lzeo,
                                  int lmx, int mq, int mb, int *mo,
                                  struct t_blocks nbpc,
                                  int *lxim, int *letm, int *lzem,
                                  int *lpos,
                                  float *vart);

void IO_write_output_file_tecplot_block(int nblocks,
                                        int lxio, int leto, int lzeo,
                                        int lmx, int mq, int mb);

int IO_write_techead(FILE * fptr, int mb, int mq);