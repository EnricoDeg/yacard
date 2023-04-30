#include <stdlib.h>
#include "config_parser.h"

static config_t cfg;
static config_setting_t *setting;

void config_open(char *filename, char *component) {
    
    config_init(&cfg);
    
    if(! config_read_file(&cfg, filename))
    {
        perror("Error reading config file\n");
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    setting = config_lookup(&cfg, component);
    if(setting == NULL) {
        perror("Error reading component in config file\n");
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }
}

double config_read_float(char *variable) {
    double value;
    if(!config_setting_lookup_float(setting, variable, &value)) {
        perror("Error reading variable component in config file\n");
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }
    return value;
}

void config_close() {
    config_destroy(&cfg);
}