#include <stdio.h>
#include <mpi.h>

void safe_mpi(int error_code);
void *safe_malloc(int size);
void safe_free (void * ptr);
void finish(char * component, char * function);
FILE *safe_open(const char *filename, const char *mode);
void safe_close(FILE *fp);