/**
 * A unit test integrating an asteroid forward, then backwards in time.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

const double au2meter = 149597870700;
#include "planets.h"
double maxdt = 0;
void hb(struct reb_simulation* r){
    //printf("%.20e %.20e\n",r->particles[0].x,r->particles[0].y);
    if (fabs(r->dt_last_done)>maxdt){
        maxdt = fabs(r->dt_last_done);
    }
}



double roundtrip(struct assist_ephem* ephem, double trange){
    struct reb_simulation* r = reb_create_simulation();
    r->heartbeat = hb;
    struct assist_extras* ax = assist_attach(r, ephem);
    //ax->ephem_cache = NULL;
    double t0 = 2458849.5-ephem->jd_ref+30000;// - 100000;
    double x0 = 3.3388753502614090e+00;
    double y0 = -9.1765182678903168e-01;
    double z0 = -5.0385906775843303e-01;
//    x0 = 3.50251071336706409909e+00;
     x0 = 2.88016863601199135658e+00;
    y0 = -1.73235639607547065033e+00;
z0 = -8.19370724376431436298e-01;
//    y0 = -5.43581485991787571876e-02;
//    z0 = -1.57076606690183528947e-01;
    //x0 = -6.63181069897987640616e-01;
    //y0 = 2.57322362393846670869e+00;
    //z0 = 1.07098547141593525289e+00;
//    x0 += 1e-6;
    ax->forces ^= ASSIST_FORCE_GR_EIH;
    ax->forces ^= ASSIST_FORCE_EARTH_HARMONICS;
    ax->forces ^= ASSIST_FORCE_SUN_HARMONICS;
    ax->forces ^= ASSIST_FORCE_ASTEROIDS;
    ax->forces ^= ASSIST_FORCE_PLANETS;
    r->t = t0; 
    // Initial conditions of asteroid Holman
    reb_add_fmt(r, "x y z vx vy vz",
        x0, y0, z0,
    5.32000300435103972568e-03, 6.46760310074497024591e-03, 2.44038701006633581073e-03
//      2.8056633153049852e-03,  7.5504086883996860e-03,  2.9800282074358684e-03
  //  1.84388098399610356366e-04, 7.92982328408079091553e-03, 3.23608315592052850004e-03
//-1.00294168284254924667e-02, -3.26658210229528386814e-03, -9.66328929492350297351e-04
    );
   // reb_integrate(r, t0+30000);
   // printf("%.20e; \n%.20e;\n%.20e; \n%.20e;\n%.20e; \n%.20e;\n",
   //         r->particles[0].x,
   //         r->particles[0].y,
   //         r->particles[0].z,
   //         r->particles[0].vx,
   //         r->particles[0].vy,
   //         r->particles[0].vz);
   // exit(0);
   
    if (1){
        // Out..
        //r->ri_ias15.epsilon = 1e-10;
        r->exact_finish_time = 0;
        reb_integrate(r,  t0 + trange);
        //reb_integrator_reset(r);
        //reb_step(r); 
        r->exact_finish_time = 1;
        reb_integrate(r,  t0);
        //printf("%d\n",r->ri_ias15.iterations_max_exceeded);
    }else{
        r->ri_ias15.epsilon=0;
        r->dt = 20;
        int steps = 0;
        while (r->t<t0+trange){
            reb_step(r);
            steps++;
        }
        int newsteps = (double)steps *1.1;
        r->dt *=-(double)steps/(double)(newsteps);
        while (newsteps>0){
            newsteps--;
            reb_step(r);
        }
    }
    printf("%.20e\t",maxdt);


//    assert(r->t == t0);

    double dx = r->particles[0].x - x0;
    double dy = r->particles[0].y - y0;
    double dz = r->particles[0].z - z0;

    double d = sqrt(dx*dx + dy*dy + dz*dz)*au2meter;
   
    assist_free(ax);
    reb_free_simulation(r);
    return d;
}

int main(int argc, char* argv[]){

    struct assist_ephem* ephem = assist_ephem_init(
            "../../data/linux_p1550p2650.440",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
    if (0){
        for (int i=0;i<1000;i++){
            double t = 4.85994979608233625186e+04+(double)i/200000;
            t = 3.35274985812470258679e+04+(double)i/200000;
          //  t = 1.60874979191289039591e+04+(double)i/200000;
         //   t = 2.29354985515976004535e+04+(double)i/200000;
           // t = 2.64554977358239084424e+04+(double)i/200000;
            int blk = (t+ephem->jd_ref - ephem->pl->beg)/ephem->pl->inc;
            struct reb_particle p = assist_get_particle(ephem, 0, t);
            printf("%.20f %.20f %d \n", t, p.x, blk);
        }
        exit(0);
    }
    //{
    //    double r0 = roundtrip(ephem, 100);
    //    printf("distance: %em\n",r0);
    //    assert(r0 < 1e-4); // required accuracy in m
    //}
    //{
    //    double r0 = roundtrip(ephem, 1000);
    //    printf("distance: %em\n",r0);
    //    assert(r0 < 1e-3); // required accuracy in m
    //}
    //{
   //     double r0 = roundtrip(ephem, 43219.423752);
   //     exit(1);
   //     printf("distance: %em\n",r0);
   //     assert(r0 < 5e-3); // required accuracy in m
   // //}
    {
        double rt = 1000;
        for (int i=0;i<1000;i++){
            rt *= 1.01;
            double r0 = roundtrip(ephem, rt);
            printf("%.20e %.20e\n",rt, r0);
        }
            //printf("distance: %em\n",r0);
//        assert(r0 < 5e-2); // required accuracy in m
    }
    assist_ephem_free(ephem);
}

