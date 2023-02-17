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

extern int ec;

int main(int argc, char* argv[]){
    double t_start = 7304.5;
    double dt = 1;
    if (argc>1){
        dt = atof(argv[1]);
    }
    int N_steps = floor(100.0/dt);
    double t_end = t_start + dt*N_steps;  
    struct reb_particle p_start = {.x=3.3388753502614090e+00, .y=-9.1765182678903168e-01, .z=-5.0385906775843303e-01,
                                   .vx=2.8056633153049852e-03, .vy=7.5504086883996860e-03, .vz=2.9800282074358684e-03}; 

    int ec_ias15 = 0;
    int ec_whfast = 0;
    // Reference solution using IAS15
    struct reb_particle p_ias15; 
    {
        ec = 0;
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
        ec_ias15 = ec;
        reb_free_simulation(r);
    }
    
    // Solution using WHFast
    struct reb_particle p_whfast; 
    {
        ec = 0;
        struct reb_simulation* r = reb_create_simulation();
        struct assist_ephem* ephem = assist_ephem_create( "../../data/linux_p1550p2650.440", "../../data/sb441-n16.bsp");
        struct assist_extras* ax = assist_attach(r, ephem);
        ax->forces ^= ASSIST_FORCE_GR_EIH;
        ax->forces ^= ASSIST_FORCE_SUN; // We should excluding Sun's gravity in assist_additional_forces to calculate the acceleration term (1/r - 1/r0) more accurately.
        r->t = t_start;
        reb_add(r, p_start);

        for (int i=0; i<N_steps; i++){
            // Drift
            struct reb_particle sun = assist_get_particle(ephem, ASSIST_BODY_SUN, r->t);
            reb_whfast_kepler_solver(r, r->particles, sun.m, 0, dt/2.);
            r->t += dt/2.;

            // Kick
            sun = assist_get_particle(ephem, ASSIST_BODY_SUN, r->t);
            struct reb_particle r0;
            r0.x = r->particles[0].x-sun.x;
            r0.y = r->particles[0].y-sun.y;
            r0.z = r->particles[0].z-sun.z;	    
            double r0_2 = r0.x*r0.x + r0.y*r0.y + r0.z*r0.z;
            double q = (sun.x*sun.x + sun.y*sun.y + sun.z*sun.z 
                    + 2*(sun.x*r0.x + sun.y*r0.y + sun.z*r0.z))/r0_2;
            double q1 = 1.0 + q;
            double fq = q*(3.0 + q*(3.0 + q))/(1.0 + q1*sqrt(q1));
            r->particles[0].ax = 0;
            r->particles[0].ay = 0;
            r->particles[0].az = 0;

            assist_additional_forces(r);
            double rb = sqrt(r->particles[0].x*r->particles[0].x + r->particles[0].y*r->particles[0].y + r->particles[0].z*r->particles[0].z);
            double GMsun_over_r3 = JPL_EPHEM_GMS/(rb*rb*rb);

            r->particles[0].ax -= (fq*r0.x - sun.x)*GMsun_over_r3;
            r->particles[0].ay -= (fq*r0.y - sun.y)*GMsun_over_r3;
            r->particles[0].az -= (fq*r0.z - sun.z)*GMsun_over_r3;

            r->particles[0].vx += r->particles[0].ax * dt;
            r->particles[0].vy += r->particles[0].ay * dt;
            r->particles[0].vz += r->particles[0].az * dt;

            // Drift
            reb_whfast_kepler_solver(r, r->particles, sun.m, 0, dt/2.);
            r->t += dt/2.;
        }


        ec_whfast = ec;
        p_whfast = r->particles[0];
        assist_free(ax);
        assist_ephem_free(ephem);
        reb_free_simulation(r);
    }
    //printf("p_start:   \tx = %.12f \ty = %.12f \tz = %.12f\n", p_start.x, p_start.y, p_start.z);
    //printf("p_ias15:   \tx = %.12f \ty = %.12f \tz = %.12f\n", p_ias15.x, p_ias15.y, p_ias15.z);
    //printf("p_whfast:  \tx = %.12f \ty = %.12f \tz = %.12f\n", p_whfast.x, p_whfast.y, p_whfast.z);
    printf("%.16e\t%.16e\t%d\t%d\n",dt, reb_particle_distance(&p_ias15, &p_whfast), ec_ias15, ec_whfast);

    // Clean up memory
}

