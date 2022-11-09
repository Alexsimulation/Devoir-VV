#ifndef EULER_FLUX_H
#define EULER_FLUX_H

#include <euler/classes.h>

// Function for pressure calc
inline real calc_p(var v, real gamma)  {
    return (gamma-1.)*(v.rhoe - 0.5/v.rho*(v.rhou*v.rhou + v.rhov*v.rhov));
}

// Function for speed of sound calc
inline real calc_c(var v, real gamma)  {
    return sqrt(gamma * calc_p(v, gamma)/v.rho);
}

inline real entropy_correction(real l, real d) {
    if (l > d) { return l; }
    else { return (l*l + d*d)/(2.*d); }
}


var convectiveFlux(
    const var& q,
    const vec2& n,
    const constants& consts
) {
    // Compute convective flux without artifical diffusion
    // Blazek page 98, eqn 4.69

    const real V = q.rhou/q.rho*n.x + q.rhov/q.rho*n.y;
    const real p = calc_p(q, consts.gamma);

    return var(
        V * q.rho,
        V * q.rhou + n.x*p,
        V * q.rhov + n.y*p,
        V * (q.rhoe + p)
    );
}

var artificialDiffusion(
    const var& qR,
    const var& qL,
    const vec2& n,
    const constants& consts
) {
    // Compute artifical diffusion
    // Blazed pages 104-105, eqns 4.93, 4.94, 4.95

    const real pR = calc_p(qR, consts.gamma);
    const real pL = calc_p(qL, consts.gamma);

    // Roe variables
    const real uR = qR.rhou/qR.rho;
    const real uL = qL.rhou/qL.rho;
    const real vR = qR.rhov/qR.rho;
    const real vL = qL.rhov/qL.rho;

    const real srhoR = sqrt(qR.rho);
    const real srhoL = sqrt(qL.rho);
    const real rho = srhoR*srhoL;
    const real u = (uL*srhoL + uR*srhoR)/(srhoL + srhoR);
    const real v = (vL*srhoL + vR*srhoR)/(srhoL + srhoR);
    const real h = ((qL.rhoe + pL)/qL.rho*srhoL + (qR.rhoe + pR)/qR.rho*srhoR)/(srhoL + srhoR);
    const real q2 = u*u + v*v;
    const real c = sqrt( (consts.gamma - 1.) * (h - 0.5*q2) );
    const real V = u*n.x + v*n.y;
    const real VR = uR*n.x + vR*n.y;
    const real VL = uL*n.x + vL*n.y;

    const real delta = 0.05*c;
    real lambda_cm = entropy_correction(std::abs(V-c), delta);
    real lambda_c  = entropy_correction(std::abs(V), delta);
    real lambda_cp = entropy_correction(std::abs(V+c), delta);

    const real kF1 = lambda_cm*((pR-pL) - rho*c*(VR-VL))/(2.*c*c);
    const real kF234_0 = lambda_c*((qR.rho - qL.rho) - (pR-pL)/(c*c));
    const real kF234_1 = lambda_c*rho;
    const real kF5 = lambda_cp*((pR-pL) + rho*c*(VR-VL))/(2*c*c);

    // Roe flux
    return var(
        kF1           + kF234_0                                                      + kF5,
        kF1*(u-c*n.x) + kF234_0*u      + kF234_1*(uR - uL - (VR-VL)*n.x)             + kF5*(u+c*n.x),
        kF1*(v-c*n.y) + kF234_0*v      + kF234_1*(vR - vL - (VR-VL)*n.y)             + kF5*(v+c*n.y),
        kF1*(h-c*V)   + kF234_0*q2*0.5 + kF234_1*(u*(uR-uL) + v*(vR-vL) - V*(VR-VL)) + kF5*(h+c*V)
    );
}


var emptyBound(
    const var& q,
    const var& b,
    const vec2& n,
    const constants& consts
) {
    return q;
}


var varWall(
    const var& q,
    const var& b,
    const vec2& n,
    const constants& consts
) {
    // Variable on a wall boundary
    vec2 U = vec2(q.rhou/q.rho, q.rhov/q.rho);
    // Flip the velocity around the edge
    U -= n*2.*dot(U, n);
    // Pressure stays the same
    const auto p = calc_p(q, consts.gamma);
    return var(
        q.rho,
        U.x*q.rho,
        U.y*q.rho,
        p/(consts.gamma - 1.) + 0.5*q.rho*(U.x*U.x + U.y*U.y)
    );
}


var varFarfield(
    const var& q,
    const var& b,
    const vec2& n,
    const constants& consts
) {
    const auto U = vec2(q.rhou/q.rho, q.rhov/q.rho);
    const auto V = dot(U, n);
    const auto c = calc_c(q, consts.gamma);
    const auto mach = norm(U)/c;

    var qb;

    if (V < 0) {
        // Flux enters, inlflow
        if (mach >= 1) {
            // Supersonic
            // Blazek p. 263 eqn 8.20
            qb = b;
        } else {
            // Subsonic
            // Blazek p. 264 eqn 8.22

            // Outside
            const real pa = calc_p(b, consts.gamma);
            const real ua = b.rhou/b.rho;
            const real va = b.rhov/b.rho;
            const real rhoa = b.rho;

            // Inside
            const real ud = q.rhou/q.rho;
            const real vd = q.rhov/q.rho;
            const real rhod = q.rho;
            const real cd = calc_c(q, consts.gamma);
            const real pd = calc_p(q, consts.gamma);

            const real pb = 0.5*(pa + pd - rhod*cd*(n.x*(ua - ud) + n.y*(va - vd)));
            const real rhob = rhoa + (pb - pa)/(cd*cd);
            const real ub = ua - n.x*(pa - pb)/(rhod*cd);
            const real vb = va - n.y*(pa - pb)/(rhod*cd);

            qb = var(
                    rhob,
                    ub*rhob,
                    vb*rhob,
                    pb/(consts.gamma-1.) + 0.5*rhob*(ub*ub + vb*vb)
                );
        }
    } else {
        // Flux exists, outflow
        if (mach >= 1) {
            // Supersonic
            // Blazek p. 264 eqn 8.21
            qb =  q;
        } else {
            // Subsonic
            // Blazek p. 264 eqn 8.23

            // Outside
            const real pa = calc_p(b, consts.gamma);
            const real ua = b.rhou/b.rho;
            const real va = b.rhov/b.rho;
            const real rhoa = b.rho;

            // Inside
            const real ud = q.rhou/q.rho;
            const real vd = q.rhov/q.rho;
            const real rhod = q.rho;
            const real cd = calc_c(q, consts.gamma);
            const real pd = calc_p(q, consts.gamma);

            //const real pb = pa;
            const real pb = pa;
            const real rhob = rhod + (pb - pd)/(cd*cd);
            const real ub = ud + n.x*(pd - pb)/(rhod*cd);
            const real vb = vd + n.y*(pd - pb)/(rhod*cd);
            qb = var(
                    rhob,
                    ub*rhob,
                    vb*rhob,
                    pb/(consts.gamma-1.) + 0.5*rhob*(ub*ub + vb*vb)
                );
        }
    }

    return qb;
}


std::map<std::string, var (*)(const var&, const var&, const vec2&, const constants&)> boundaryFuncMap() {
    return {
        {"null", &emptyBound},
        {"wall", &varWall},
        {"farfield", &varFarfield}
    };
}






#endif
