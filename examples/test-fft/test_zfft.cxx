/*
 * Test ZFFT 
 *
 */

#include <bout.hxx>
#include <boutmain.hxx>

#include <bout/constants.hxx>
#include <cyclic_reduction.hxx>
#include <dcomplex.hxx>
#include <fft.hxx>

#include <msg_stack.hxx>
#include <bout/assert.hxx>
#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

// Change this to dcomplex to test complex matrix inversion
typedef BoutReal T;

int physics_init(bool restarting) {

	Options *options = Options::getRoot();
	ASSERT2(mesh->xstart > 0); // Need at least one guard cell

	Field3D result;
	Field3D f = 1.0;
	BoutReal zsmooth = -1;		// < 0 for no smoothing
	result.allocate();

	BoutReal ***fd, ***rd;
	fd = f.getData();
	rd = result.getData();

	int ncz = mesh->ngz-1;

	for(int jx=0;jx<mesh->ngx;jx++){
		for(int jz=0;jz<=ncz;jz++) {
			BoutReal kwave=jx*jz*2.0*PI/ncz;
			fd[jx][0][jz] = cos(kwave);
		}
	}

	//.allocate memory
	static dcomplex **ft = (dcomplex**) NULL;
	ft = cmatrix(mesh->ngx, ncz/2 + 1);

	// Loop over all y indices
	for(int jx=0;jx<mesh->ngx;jx++){
		ZFFT(fd[jx][0], mesh->zShift(jx, 0), ft[jx]);
	}

	BoutReal tol;
	OPTION(options, tol, 1e-5);


	// Check result

	int passed = 1;
	// Only check up to ncz/2... after that it gets non-trivial to sort out box order
	for(int jx=0;jx<min(mesh->ngx,ncz/2);jx++){
		for(int jz=0;jz<=ncz/2;jz++) {
			dcomplex val = 0.0;
			if( (jx%ncz) == 0 && (jx%ncz) == jz ) val = 1.0;
			else if( (jx%ncz) == jz ) val = 0.5;
			if(abs(val - ft[jx][jz]) > tol){
				output << "\t(" << jx << "," << jz << ") : " << ft[jx][jz] << " ?= " << val << endl;
				passed = 0;
			}
		}
	}

	int allpassed;
	MPI_Allreduce(&passed, &allpassed, 1, MPI_INT, MPI_MIN, BoutComm::get());

	// Saving state to file
	SAVE_ONCE(allpassed);

	output << "******* ZFFT test case: ";
	if(allpassed) {
		output << "PASSED" << endl;
	}else 
		output << "FAILED" << endl;

	// Write data to file
	dump.write();
	dump.close();

	MPI_Barrier(BoutComm::get());

	// Send an error code so quits
	return 1;
}

int physics_run(BoutReal t) {
	// Doesn't do anything
	return 1;
}
