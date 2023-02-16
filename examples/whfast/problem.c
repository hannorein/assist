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
void reb_whfast_kepler_solver(const struct reb_simulation* const r, struct reb_particle* const restrict p_j, const double M, unsigned int i, double _dt); 
void assist_additional_forces(struct reb_simulation* r);

int main(int argc, char* argv[]){
    double t_start = 7304.5;
    double dt = 1;
    int N_steps = 100;
    double t_end = t_start + dt*N_steps;  
    struct reb_particle p_start = {.x=3.3388753502614090e+00, .y=-9.1765182678903168e-01, .z=-5.0385906775843303e-01,
                                   .vx=2.8056633153049852e-03, .vy=7.5504086883996860e-03, .vz=2.9800282074358684e-03}; 

    // Reference solution using IAS15
    struct reb_particle p_ias15; 
    {
        struct reb_simulation* r = reb_create_simulation();
        struct assist_ephem* ephem = assist_ephem_create( "../../data/linux_p1550p2650.440", "../../data/sb441-n16.bsp");
        struct assist_extras* ax = assist_attach(r, ephem);
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
        ax->forces ^= ASSIST_FORCE_GR_EIH;
        // ax->forces ^= ASSIST_FORCE_SUN; // We should excluding Sun's gravity in assist_additional_forces to calculate the acceleration term (1/r - 1/r0) more accurately.
        r->t = t_start;
        reb_add(r, p_start);

        for (int i=0; i<N_steps; i++){
            // Drift
            struct reb_particle sun = assist_get_particle(ephem, ASSIST_BODY_SUN, r->t);
            reb_whfast_kepler_solver(r, r->particles, sun.m, 0, dt/2.);
            r->t += dt/2.;

            // Kick
            sun = assist_get_particle(ephem, ASSIST_BODY_SUN, r->t);
            r->particles[0].ax = 0;
            r->particles[0].ay = 0;
            r->particles[0].az = 0;
            assist_additional_forces(r);
            double rb = sqrt(r->particles[0].x*r->particles[0].x + r->particles[0].y*r->particles[0].y + r->particles[0].z*r->particles[0].z);
            r->particles[0].ax += JPL_EPHEM_GMS/(rb*rb*rb) *r->particles[0].x;
            r->particles[0].ay += JPL_EPHEM_GMS/(rb*rb*rb) *r->particles[0].y;
            r->particles[0].az += JPL_EPHEM_GMS/(rb*rb*rb) *r->particles[0].z;

            r->particles[0].vx += r->particles[0].ax * dt;
            r->particles[0].vy += r->particles[0].ay * dt;
            r->particles[0].vz += r->particles[0].az * dt;

            // Drift
            reb_whfast_kepler_solver(r, r->particles, sun.m, 0, dt/2.);
            r->t += dt/2.;
        }


        p_whfast = r->particles[0];
        assist_free(ax);
        assist_ephem_free(ephem);
        reb_free_simulation(r);
    }
    printf("p_start:   \tx = %.12f \ty = %.12f \tz = %.12f\n", p_start.x, p_start.y, p_start.z);
    printf("p_ias15:   \tx = %.12f \ty = %.12f \tz = %.12f\n", p_ias15.x, p_ias15.y, p_ias15.z);
    printf("p_whfast:  \tx = %.12f \ty = %.12f \tz = %.12f\n", p_whfast.x, p_whfast.y, p_whfast.z);

    // Clean up memory
}

