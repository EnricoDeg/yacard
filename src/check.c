#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "check.h"

void *safe_malloc(int size) {
    void *ptr = malloc(size);

    if (!ptr && (size > 0)) {
      fprintf(stderr, "malloc failed!\n");
      exit(EXIT_FAILURE);
    }

    return ptr;
}

void safe_free (void * ptr)
{
    if (ptr != NULL) {
        free (ptr);
    } else {
    	int myid;
        safe_mpi( MPI_Comm_rank(MPI_COMM_WORLD, &myid) );
        fprintf(stderr, "%d: null pointer to free!\n", myid);
        exit(EXIT_FAILURE);
    }
}

FILE *safe_open(const char *filename, const char *mode) {
    FILE *fp = fopen(filename, mode);
    if (!fp) {
        printf("Cannot open file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    return fp;
}

void safe_close(FILE *fp) {
    if (!fp) {
        fprintf(stderr, "Error closing file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    // /fclose(fp);
    if(fclose(fp) != 0) {
        fprintf(stderr, "Error closing file: %s", strerror(errno));
    }
}

void safe_mpi(int error_code) {
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    char error_string[MPI_MAX_ERROR_STRING];
    int length_of_error_string, error_class;

    if (error_code != MPI_SUCCESS) {
        MPI_Comm_rank(comm, &rank);
        MPI_Error_class(error_code, &error_class);
        MPI_Error_string(error_class, error_string, &length_of_error_string);
        fprintf(stderr, "%3d: %s\n", rank, error_string);
        MPI_Error_string(error_code, error_string, &length_of_error_string);
        fprintf(stderr, "%3d: %s\n", rank, error_string);
        MPI_Abort(comm, error_code);
    }

}

void finish(char * component, char * function) {
    fprintf(stderr,"Error in %s component: %s\n", component, function);
    exit(EXIT_FAILURE);

}