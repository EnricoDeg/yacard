#include <mpi.h>

void safe_mpi(int error_code);
void *safe_malloc(int size);
void safe_free (void * ptr);