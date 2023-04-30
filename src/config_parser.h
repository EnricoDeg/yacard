#include <libconfig.h>

void config_open(char *filename, char *component);
double config_read_float(char *variable);
void config_close();