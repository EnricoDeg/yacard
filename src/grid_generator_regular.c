#include "grid_generator_regular.h"
#include "config_parser.h"

static char component[] = "grid_generator";
static int lxi0;

void grid_generator_regular_read_input(char *filename) {
    config_open(filename, component);
    lxi0 = config_read_int("lxi0");
    config_close();
}