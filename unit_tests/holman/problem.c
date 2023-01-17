/**
 * This example demonstrates how to use ASSIST to integrate a test particle in 
 * the field of the Sun, planets, moon, and a set of massive asteroids, whose
 * positions come from JPL's DE440/441 ephemeris.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    // Initial conditions of asteroid
    int n_particles = 1;    
    double* state = malloc(n_particles*6*sizeof(double));

    state[0] = -2.724183384883979E+00 ;  // x in AU
    state[1] = 3.613497537656375E-03 ;  // y
    state[2] = 9.692679446383898E-02 ;  // z
    state[3] = -1.374545432301129E-04 ;  // vx in AU/day
    state[4] = -1.109218154000629E-02 ;  // vy
    state[5] = 2.360000345736158E-04 ;  // vz
   
    double jd_ref = 2451545.0;  // All times are emasured relative to this Julian date 
    double tstart = 8416.5;     // Initial time 
    double tend = 8446.5 ;       // Final time
    double tstep = 20;          // Integration step size. 
    
    
    // Interpolate between integration steps on these substeps 
    int n_substeps = 4;
    double hg[n_substeps+1];
    for(int i=0; i<=n_substeps; i++){
        hg[i]=(1.0/n_substeps)*i;
    }

    // Allocate memory for output.
    int n_alloc = ((tend-tstart)/tstep + 1)*n_substeps;
    double* outstate = (double *) malloc((n_alloc)*6*sizeof(double));
    double* outtime  = (double *) malloc((n_alloc)*sizeof(double));

    int n_steps_done;
    int status = integration_function(
            jd_ref,                 // I do not understand what this variable does
            tstart, tend, tstep,    // Time range of integration 
            0,                      // 1=geocentric, 0=barycentric
            1e-9,                   // epsilon
            n_particles,            // number of particles
            state,                  // initial conditions
            NULL,                   // additional particle parameters for non-gravitational effects (not used in this example)
            0, NULL, NULL, NULL,    // No variational particles in this example
            n_alloc,                // Allocated space in output array
            &n_steps_done,          // Number of steps actually done
            n_substeps,             // Number of substeps 
            hg,                     // Location of substeps (in fractions of an integration time step)
            outtime,                // Output times
            outstate,               // Output states
            0.0                     // Minimum dt in integration
            ); 
    if (status != REB_EXIT_SUCCESS) {
       printf("Warning! Simulation did *not* run successfully.\n");
    } 

    // Number of outputs generated
    int n_outs = n_steps_done*n_substeps + 1;
    
    // Check final time
    assert(outtime[n_outs-1] == tend);

    // Batch file to generate the JPL data (last queried on Jan 16th 2023)
    // !$$SOF
    // MAKE_EPHEM=YES
    // COMMAND='DES=2003666;'
    // EPHEM_TYPE=VECTORS
    // CENTER='500@0'
    // START_TIME='2023-01-17'
    // STOP_TIME='2023-02-16'
    // STEP_SIZE='1 DAYS'
    // VEC_TABLE='3'
    // REF_SYSTEM='ICRF'
    // REF_PLANE='ECLIPTIC'
    // VEC_CORR='NONE'
    // OUT_UNITS='AU-D'
    // VEC_LABELS='YES'
    // VEC_DELTA_T='NO'
    // CSV_FORMAT='NO'
    // OBJ_DATA='YES'


    double* jpl_horizons_output = malloc(n_particles*6*sizeof(double));
    jpl_horizons_output[0] = -2.710320457933958E+00;
    jpl_horizons_output[1] = -3.284425995373661E-01;
    jpl_horizons_output[2] = 1.033508308499156E-01;
    jpl_horizons_output[3] = 1.059255302926290E-03;
    jpl_horizons_output[4] = -1.102056611134586E-02;
    jpl_horizons_output[5] = 1.918473889773432E-04;
    
    int offset = (n_outs-1)*6;

    for (int i=0;i<6;i++){
        printf("diff: %.20f  %.20f   %.20f\n", state[i], outstate[i+offset],  jpl_horizons_output[i] - outstate[i+offset]);
    }    
        
}
