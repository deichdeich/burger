/* Initialization routine for homework #2 (C version)
 *
 * Arguments:
 *
 *    rho(:,:)     (real) array to receive initial density field
 *    P(:,:)       (real) array to receive initial pressure field
 *    Etot(:,:)    (real) array to receive initial total energy density field
 *    ux(:,:)      (real) array to receive initial x-velocity field
 *    uy(:,:)      (real) array to receive initial y-velocity field
 *    Nx, Ny       (integer) number of interior zones in x- and y-directions
 *    Nbz          (integer) number of boundary zones
 *    dx, dy       (real) cell sizes in x- and y-directions
 *    gamma        (real) gas ratio of specific heats
 *    E            (real) total explosion energy
 *    rho_ambient  (real) initial ambient gas density
 *    p_ambient    (real) initial ambient gas pressure
 *
 * Arrays are stored in Fortran order, ie. first index cycles fastest.
 * In C they look like buffers; the buffer index must be computed thus:
 *
 *   rho(i,j) <--> rho[(j+Nbz-1)*(Nx+2*Nbz) + i+Nbz-1]
 *
 *   i = 1-Nbz ... Nx+Nbz
 *   j = 1-Nbz ... Ny+Nbz
 *
 * Only the interior zones need to be initialized.
 *
 */

#ifdef IBM
#define FTOC(x) x
#else
#define FTOC(x) x##_
#endif

void FTOC(init)(double *rho, double *P, double *Etot, double *ux, double *uy,
            int *Nx, int *Ny, int *Nbz, double *dx, double *dy, double *gamma,
            double *E, double *rho_ambient, double *p_ambient) {

	for(int i =1 - *Nbz; i < (*Nx + *Nbz); i++){
		for(int j = 1 - *Nbz; j < (*Ny + *Nbz); j++){
			int index = (j + (*Nbz) - 1) * (*Nx + 2 * (*Nbz)) + i + (*Nbz) - 1;
			rho[index] = *rho_ambient;
			P[index] = *p_ambient;
				
			ux[index] = 0.;
			uy[index] = 0.;
		}
	}

	int xcenter = *Nx/2;
	int ycenter = *Ny/2;
	int pert_xdim = 2;
	int pert_ydim = 2;

	for(int i = xcenter - (pert_xdim / 2); i < (xcenter + (pert_xdim / 2)); i++){
		for(int j = (ycenter - (pert_ydim / 2)); j < (ycenter + (pert_ydim / 2)); j++){
			int index = (j + (*Nbz) - 1) * (*Nx + 2 * (*Nbz)) + i + (*Nbz)-1;
			P[index] = (*gamma - 1.) * (*E) /(pert_xdim * pert_ydim * (*dx) * (*dy));
		}
	} 
	for(int i =1 - *Nbz; i < (*Nx + *Nbz) ; i++){
		for(int j = 1 - *Nbz; j < (*Ny + *Nbz); j++){
			int index = (j + (*Nbz) - 1) * (*Nx + 2 * (*Nbz)) + i + (*Nbz) - 1;
			Etot[index] = P[index] / rho[index] * 1. / (*gamma - 1.);
		}
	}

}
