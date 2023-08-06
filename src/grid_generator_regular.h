// public interface
void grid_generator_regular_read_input(char *filename);

void grid_generator_regular_go(int mbk, int mb, int lmx, int *mo, int *lpos, int *lxim,
                               int *letm, int *lzem, int lxio, int leto);

double* grid_generator_regular_get_x();

double* grid_generator_regular_get_y();

double* grid_generator_regular_get_z();