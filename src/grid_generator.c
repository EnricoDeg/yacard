#include "grid_generator.h"
#include "check.h"
#include "config_parser.h"
#include "grid_generator_regular.h"

static int grid_type;
static char component[] = "grid_generator";

void grid_generator_read_input(char *filename) {
    static char function[] = "grid_generator_read_input";

    config_open(filename, component);
    grid_type = config_read_int("grid_type");
    config_close();
    if (grid_type == GRID_REGULAR) {
        grid_generator_regular_read_input(filename);
    } else {
        finish(component, function);
    }
}