struct t_field {
    double *ptr;
};

struct t_physics_fields {
    struct t_field variable[5];
};

struct t_grid_fields {
    struct t_field coordinate[3];
};

struct t_coordinates {
    int *ptr;
};

struct t_blocks {
    struct t_coordinates coordinate[3]; 
};

struct t_directions {
    struct t_coordinates coordinate[3];
};

struct t_blocks_boundary {
    struct t_directions direction[2];
};

struct t_coordinate {
    int value;
};

struct t_direction {
    struct t_coordinate coordinate[3];
};

struct t_subdomain_boundary {
    struct t_direction direction[2];
};
