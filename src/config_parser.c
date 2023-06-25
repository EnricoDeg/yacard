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

void config_read_array_int(char *variable, int *array) {
    
    setting = config_lookup(&cfg, variable);
    if(setting == NULL) {
        perror("Error reading variable in config file\n");
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    if(config_setting_is_array(setting) != CONFIG_TRUE) {
        perror("Variable is NOT array\n");
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    int array_length = config_setting_length(setting);
    for (int i=0; i<array_length; i++)
        array[i] = config_setting_get_int_elem(setting, i);
}

void config_close() {
    config_destroy(&cfg);
}