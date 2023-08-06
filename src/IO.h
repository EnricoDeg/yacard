// Supported IO backends
#define IO_TECPLOT 0
#define IO_VTK 1
#define IO_SILO 2

void IO_read_input(char *filename);

void IO_write_output_file();

void IO_write_restart_file(int mbk, int mb, int lmx,
	                       int *mo, int *lpos,
	                       int *lxim, int *letm, int *lzem,
	                       int *lximb, int *letmb, int *lzemb, 
	                       int lxio, int leto, 
	                       int *current_step,
	                       int *ndt, double *dt, 
	                       double *dts, double *dte, double *timo, 
	                       double *restart_local_buffer);

void IO_read_restart_file(int mbk, int mb, int lmx,
	                       int *mo, int *lpos,
	                       int *lxim, int *letm, int *lzem,
	                       int *lximb, int *letmb, int *lzemb, 
	                       int lxio, int leto, 
	                       int *current_step,
	                       int *ndt, double *dt, 
	                       double *dts, double *dte, double *timo, 
	                       double *restart_local_buffer);