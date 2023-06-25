#include <libconfig.h>

void config_open(char *filename, char *component);
double config_read_float(char *variable);
void config_read_array_int(char *variable, int *array);
void config_close();