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
const static double reb_whfast_corrector_a_1 = 0.41833001326703777398908601289259374469640768464934;
const static double reb_whfast_corrector_a_2 = 0.83666002653407554797817202578518748939281536929867;
const static double reb_whfast_corrector_a_3 = 1.2549900398011133219672580386777812340892230539480;
const static double reb_whfast_corrector_a_4 = 1.6733200530681510959563440515703749787856307385973;
const static double reb_whfast_corrector_a_5 = 2.0916500663351888699454300644629687234820384232467;
const static double reb_whfast_corrector_a_6 = 2.5099800796022266439345160773555624681784461078960; 
const static double reb_whfast_corrector_a_7 = 2.9283100928692644179236020902481562128748537925454;
const static double reb_whfast_corrector_a_8 = 3.3466401061363021919126881031407499575712614771947;
const static double reb_whfast_corrector_b_178 = 0.093056103771425958591541059067553547100903397724386; 
const static double reb_whfast_corrector_b_177 = -0.065192863576377893658290760803725762027864651086787; 
const static double reb_whfast_corrector_b_176 = 0.032422198864713580293681523029577130832258806467604; 
const static double reb_whfast_corrector_b_175 = -0.012071760822342291062449751726959664253913904872527; 
const static double reb_whfast_corrector_b_174 = 0.0033132577069380655655490196833451994080066801611459; 
const static double reb_whfast_corrector_b_173 = -0.00063599983075817658983166881625078545864140848560259; 
const static double reb_whfast_corrector_b_172 = 0.000076436355227935738363241846979413475106795392377415; 
const static double reb_whfast_corrector_b_171 = -0.0000043347415473373580190650223498124944896789841432241; 


extern int ec;

void kepler_step(struct reb_simulation* r, double dt){
    struct assist_ephem* ephem = ((struct assist_extras*)r->extras)->ephem;
    struct reb_particle sun = assist_get_particle(ephem, ASSIST_BODY_SUN, r->t);
    reb_whfast_kepler_solver(r, r->particles, sun.m, 0, dt);
    r->t += dt;
}

void interaction_step(struct reb_simulation* r, double dt){
    struct assist_ephem* ephem = ((struct assist_extras*)r->extras)->ephem;
    struct reb_particle sun = assist_get_particle(ephem, ASSIST_BODY_SUN, r->t);
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
}


static void reb_whfast_corrector_Z(struct reb_simulation* r, const double a, const double b){
    kepler_step(r, a);
    interaction_step(r, -b);
    kepler_step(r, -2.*a);
    interaction_step(r, b);
    kepler_step(r, a);
}

void reb_whfast_apply_corrector(struct reb_simulation* r, double inv){
    const double dt = r->dt;
    // Seventeenth order corrector
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_8*dt,-inv*reb_whfast_corrector_b_171*dt);
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_7*dt,-inv*reb_whfast_corrector_b_172*dt);
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_6*dt,-inv*reb_whfast_corrector_b_173*dt);
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_5*dt,-inv*reb_whfast_corrector_b_174*dt);
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_4*dt,-inv*reb_whfast_corrector_b_175*dt);
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_3*dt,-inv*reb_whfast_corrector_b_176*dt);
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_177*dt);
    reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_178*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_178*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_177*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_3*dt,inv*reb_whfast_corrector_b_176*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_4*dt,inv*reb_whfast_corrector_b_175*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_5*dt,inv*reb_whfast_corrector_b_174*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_6*dt,inv*reb_whfast_corrector_b_173*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_7*dt,inv*reb_whfast_corrector_b_172*dt);
    reb_whfast_corrector_Z(r, reb_whfast_corrector_a_8*dt,inv*reb_whfast_corrector_b_171*dt);
}

int main(int argc, char* argv[]){
    double t_start = 7304.5;
    double dt = 1;
    int correctors = 0;
    if (argc>1){
        dt = atof(argv[1]);
    }
    if (argc>2){
        correctors = atoi(argv[2]);
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
        r->dt = dt;
        reb_add(r, p_start);

        if (correctors) reb_whfast_apply_corrector(r, 1.);
        
        for (int i=0; i<N_steps; i++){
            // Drift
            kepler_step(r, r->dt/2.);

            // Kick
            interaction_step(r, r->dt);

            // Drift
            kepler_step(r, r->dt/2.);
        }

        if (correctors) reb_whfast_apply_corrector(r, -1.);

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

