/**
 * WHFast
 *
 * An ASSIST integration using WHFast
 */
#include <stdio.h>
#include <stdlib.h>
#include "rebound.h"
#include "assist.h"
#include "const.h"

int main(int argc, char* argv[]){
    double t_start = 7304.5;
    double dt = 1;
    if (argc>1){
        dt = atof(argv[1]);
    }
    double t_end = t_start + 100.;  
    struct reb_particle p_start = {.x=3.3388753502614090e+00, .y=-9.1765182678903168e-01, .z=-5.0385906775843303e-01,
                                   .vx=2.8056633153049852e-03, .vy=7.5504086883996860e-03, .vz=2.9800282074358684e-03}; 

    // Reference solution using IAS15
    struct reb_particle p_ias15; 
    {
        struct reb_simulation* r = reb_create_simulation();
        struct assist_ephem* ephem = assist_ephem_create( "../../data/linux_p1550p2650.440", "../../data/sb441-n16.bsp");
        struct assist_extras* ax = assist_attach(r, ephem);
        
        // Normal integration 
        ax->forces ^= ASSIST_FORCE_GR_EIH;
        r->t = t_start;
        reb_add(r, p_start);
        reb_integrate(r, t_end);
        p_ias15 = r->particles[0];
        assist_free(ax);
        assist_ephem_free(ephem);
        reb_free_simulation(r);
    }
    
    // Solution using WHFast
    struct reb_particle p_whfast; 
    {
        struct reb_simulation* r = reb_create_simulation();
        struct assist_ephem* ephem = assist_ephem_create( "../../data/linux_p1550p2650.440", "../../data/sb441-n16.bsp");
        struct assist_extras* ax = assist_attach(r, ephem);

        // Setup WHFast
        r->dt = dt;
        r->integrator = REB_INTEGRATOR_WHFAST;
        r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_ASSIST;
        r->ri_whfast.m0assist = JPL_EPHEM_GMS;
        r->ri_whfast.safe_mode = 0;
        if (argc>2){
            r->ri_whfast.corrector =  atoi(argv[2]);
        }

        // Normal integration 
        ax->forces ^= ASSIST_FORCE_GR_EIH;
        r->t = t_start;
        reb_add(r, p_start);
        reb_integrate(r, t_end);
        p_whfast = r->particles[0];
        assist_free(ax);
        assist_ephem_free(ephem);
        reb_free_simulation(r);
    }
    printf("p_start:   \tx = %.12f \ty = %.12f \tz = %.12f\n", p_start.x, p_start.y, p_start.z);
    printf("p_ias15:   \tx = %.12f \ty = %.12f \tz = %.12f\n", p_ias15.x, p_ias15.y, p_ias15.z);
    printf("p_whfast:  \tx = %.12f \ty = %.12f \tz = %.12f\n", p_whfast.x, p_whfast.y, p_whfast.z);
    printf("%.16e\t%.16e\t%d\t%d\n",dt, reb_particle_distance(&p_ias15, &p_whfast), 0, 0);

    // Clean up memory
}

