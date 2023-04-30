#include <mpi.h>

#define check_mpi_call(call, comm)              \
  do {                                          \
    int error_code = (call);                    \
    if (error_code != MPI_SUCCESS)              \
      check_mpi_error(error_code, comm);        \
  } while(0)

void check_mpi_error(int error_code, MPI_Comm comm);
