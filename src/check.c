#include <stdio.h>
#include <stdlib.h>
#include "check.h"

void *safe_malloc(int size) {
    void *ptr = malloc(size);

    if (!ptr && (size > 0)) {
      perror("malloc failed!");
      exit(EXIT_FAILURE);
    }

    return ptr;
}

void safe_free (void * ptr)
{
    if (ptr != NULL)
        free (ptr);
}