struct t_field {
    double *ptr;
};

struct t_physics_fields {
    struct t_field variable[5];
};

struct t_grid_fields {
    struct t_field coordinate[3];
};