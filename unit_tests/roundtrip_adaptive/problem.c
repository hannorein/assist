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

double roundtrip(struct assist_ephem* ephem, double trange){
    struct reb_simulation* r = reb_create_simulation();
    struct assist_extras* ax = assist_attach(r, ephem);
    double t0 = 2458849.5-ephem->jd_ref;
    double x0 = 3.3388753502614090e+00;
    double y0 = -9.1765182678903168e-01;
    double z0 = -5.0385906775843303e-01;

    r->t = t0; 

    // Initial conditions of asteroid Holman
    reb_add_fmt(r, "x y z vx vy vz",
        x0, y0, z0,
        2.8056633153049852e-03,  7.5504086883996860e-03,  2.9800282074358684e-03);
   
    // Out..
    reb_integrate(r,  t0 + trange);
    // ..and back
    reb_integrate(r,  t0);

    assert(r->t == t0);

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
    {
        double r0 = roundtrip(ephem, 43219.423752);
    //    printf("distance: %em\n",r0);
        assert(r0 < 5e-3); // required accuracy in m
    }
    //{
    //    double r0 = roundtrip(ephem, 100000);
    //    printf("distance: %em\n",r0);
    //    assert(r0 < 5e-2); // required accuracy in m
    //}
    assist_ephem_free(ephem);
}

