#ifndef EULER_LIMITER_H
#define EULER_LIMITER_H

#include <euler/classes.h>



var r_limiter(
    varVec2 g0,
    var q_min,
    var q_max,
    var q,
    vec2 de,
    real area
) {
    // http://tetra.mech.ubc.ca/ANSLab/publications/michalak2008.pdf
    const var qij = dot(g0, de) + q;
    var r;
    for (int i=0; i<q.size(); ++i) {
        if ((qij[i] - q[i]) > 1.e-14) {
            r[i] = (q_max[i] - q[i])/(qij[i] - q[i]);
        } else if ((qij[i] - q[i]) < -1.e-14) {
            r[i] = (q_min[i] - q[i])/(qij[i] - q[i]);
        } else {r[i] = 1.;}
    }
    return r;
}

/*
var limiterVenkatakrishnan(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    // Venkatakrishnan’s limiter
    var r = r_limiter(g0, q_min, q_max, q, de, area);

    const var dqij = dot(g0, de);
    const real scale = std::pow(5.*sqrt(area), 1.5);

    var lim = (r*r + r*2.)/(r*r + r + 2.);

    for (int i=0; i<q.size(); ++i) {
        if (std::abs(dqij[i]) <= scale) {
            lim[i] = 1.;
        }
    }
    return lim;
}
*/

var limiterVenkatakrishnan(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    // Venkatakrishnan’s limiter
    const real K = 5.;

    const real e2 = std::pow(K*sqrt(area), 3);
    const var d2 = dot(g0, de);
    const var dmax = q_max - q;
    const var dmin = q_min - q;

    var lim;
    const real tol = 1e-14;
    for (int i=0; i<q.size(); ++i) {
        if (d2[i] > tol) {
            lim[i] = 1./d2[i]*((dmax[i]*dmax[i] + e2)*d2[i] + 2.*d2[i]*d2[i]*dmax[i])/(dmax[i]*dmax[i] + 2*d2[i]*d2[i] + dmax[i]*d2[i] + e2);
        } else if (d2[i] < -tol) {
            lim[i] = 1./d2[i]*((dmin[i]*dmin[i] + e2)*d2[i] + 2.*d2[i]*d2[i]*dmin[i])/(dmin[i]*dmin[i] + 2*d2[i]*d2[i] + dmin[i]*d2[i] + e2);
        } else {
            lim[i] = 1.;
        }
    }
    return lim;
}

var limiterVanLeer(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    var r = r_limiter(g0, q_min, q_max, q, de, area);
    return (std::abs(r) + r)/(std::abs(r) + 1.);
}

var limiterSuperBee(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    var r = r_limiter(g0, q_min, q_max, q, de, area);
    return std::max(std::max(var(0.), std::min(r*2., var(1.))), std::min(r, var(2.)));
}

var limiterOverBee(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    var r = r_limiter(g0, q_min, q_max, q, de, area);
    return std::max(var(0.), std::min(r*2., var(2.)));
}

var limiterVanAlbada(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    var r = r_limiter(g0, q_min, q_max, q, de, area);
    return (r*r + r)/(r*r + 1.);
}

var limiterMinMod(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    var r = r_limiter(g0, q_min, q_max, q, de, area);
    return std::max(var(0.), std::min(var(1.), r));
}


var limiterOne(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    return var(1.);
}


var limiterZero(
    const varVec2& g0,
    const var& q_min,
    const var& q_max,
    const var& q,
    const vec2& de,
    const real& area,
    const vec2& n
) {
    return var(0.);
}




#endif

