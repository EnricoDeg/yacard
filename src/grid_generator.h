// Supported grid backends
#define GRID_REGULAR 0

// public interface
void grid_generator_read_input(char *filename);

void grid_generator_go(int mbk, int mb, int lmx, int *mo, int *lpos, int *lxim,
                       int *letm, int *lzem, int lxio, int leto);
