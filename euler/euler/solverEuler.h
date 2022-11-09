#ifndef EULER_SOLVER_EULER_H
#define EULER_SOLVER_EULER_H

#include <euler/solver.h>



class eulerSolver : public solver {
public:

    void calc_qt() {
        internal_calc_qt(qt, q);
    }

    virtual bool isTransient() {}
    virtual void dtFromCfl() {}
    virtual void update_cells() {}

};

class transientEulerSolver : public eulerSolver {
public:

    transientEulerSolver(
        const mesh& m_in,
        const std::vector<inputVar>& q0, 
        const constants cin
    ) {
        tolerance = 1.e-14;
        c = cin;
        q.resize(q0.size());
        for (int i=0; i<q0.size(); ++i) {
            q[i] = convertVariable(q0[i], c);
        }
        qt.resize(q.size());
        dtScale.resize(q.size());

        set_mesh(m_in);

        solver_type = "Euler";
    }

    bool isTransient() {return true;}

    void dtFromCfl() {
        // Reset dt scale
        #pragma omp simd
        for (int i=0; i<dtScale.size(); ++i) {
            dtScale[i] = 0.;
        }
        // Loop over edges
        for (int i=0; i<m.edgesLengths.size(); ++i) {
            // Blazek p. 175 eqn. 6.21
            const auto fL = m.edgesFaces0[i];
            const auto fR = m.edgesFaces1[i];
            const auto n = m.edgesNormals[i];

            const var qe = (q[fL] + q[fR])*0.5;
            const real dtScale_i = m.edgesLengths[i]*(std::abs(qe.rhou/qe.rho*n.x + qe.rhov/qe.rho*n.y) + calc_c(qe, c.gamma));
            dtScale[fL] += dtScale_i;
            dtScale[fR] += dtScale_i;
        }
        for (int i=0; i<m.facesAreas.size(); ++i) {
            // Blazek p. 175 eqn. 6.20
            const real dti = cfl * m.facesAreas[i] / dtScale[i];
            if ((dti < dt)||(i == 0)) {
                dt = dti;
            }
        }
    }

    void update_cells() {
        for (int i=0; i<m.facesAreas.size(); ++i) {
            q[i] += qt[i] * dt;
        }
    }

};

class steadyEulerSolver : public eulerSolver {
public:

    steadyEulerSolver(
        const mesh& m_in,
        const std::vector<inputVar>& q0, 
        const constants cin
    ) {
        tolerance = 1.0e-4;
        c = cin;
        q.resize(q0.size());
        for (int i=0; i<q0.size(); ++i) {
            q[i] = convertVariable(q0[i], c);
        }
        qt.resize(q.size());
        dtScale.resize(q.size());
        dtv.resize(q.size());

        set_mesh(m_in);

        solver_type = "Euler";
    }

    bool isTransient() {return false;}

    void dtFromCfl() {
        // Local time stepping
        // Reset dt scale
        #pragma omp simd
        for (int i=0; i<dtScale.size(); ++i) {
            dtScale[i] = 0.;
        }
        // Loop over edges
        for (int i=0; i<m.edgesLengths.size(); ++i) {
            // Blazek p. 175 eqn. 6.21
            const auto fL = m.edgesFaces0[i];
            const auto fR = m.edgesFaces1[i];
            const auto n = m.edgesNormals[i];

            const var qe = (q[fL] + q[fR])*0.5;
            const real dtScale_i = m.edgesLengths[i]*(std::abs(qe.rhou/qe.rho*n.x + qe.rhov/qe.rho*n.y) + calc_c(qe, c.gamma));
            dtScale[fL] += dtScale_i;
            dtScale[fR] += dtScale_i;
        }

        #pragma omp simd
        for (int i=0; i<m.facesAreas.size(); ++i) {
            // Blazek p. 175 eqn. 6.20
            dtv[i] = cfl * m.facesAreas[i] / dtScale[i];
        }
    }

    void update_cells() {
        #pragma omp simd
        for (int i=0; i<m.facesAreas.size(); ++i) {
            q[i] += qt[i] * dtv[i];
        }
    }

};



#endif


