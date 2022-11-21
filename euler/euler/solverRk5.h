#ifndef EULER_SOLVER_RK5_H
#define EULER_SOLVER_RK5_H

#include <euler/solver.h>




class rk5Solver : public solver {
public:
    std::vector<var> k;
    std::vector<var> qk;
    std::array<real, 5> alphas = {0.0533, 0.1263, 0.2375, 0.4414, 1.0};

    virtual void calc_qt() {}
    virtual bool isTransient() {}
    virtual void dtFromCfl() {}
    virtual void update_cells() {}

};


class transientRk5Solver : public rk5Solver {
public:

    transientRk5Solver(
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

        k.resize(q.size());

        qk.resize(q.size());

        set_mesh(m_in);

        solver_type = "Rk5";
    }

    void calc_qt() {

        for (int i=0; i<qk.size(); ++i) {
            qk[i] = q[i];
        }

        for (int i=0; i<(alphas.size() - 1); ++i) {
            internal_calc_qt(k, qk);

            for (int j=0; j<qk.size(); ++j) {
                qk[j] = q[j] + k[j] * dt * alphas[i];
            }
        }
        internal_calc_qt(k, qk);

        for (int i=0; i<qt.size(); ++i) {
            qt[i] = k[i] * alphas[alphas.size()-1] ;
        }
    }

    bool isTransient() {return true;}

    void dtFromCfl() {
        // Reset dt scale
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

class steadyRk5Solver : public rk5Solver {
public:

    steadyRk5Solver() {}

    steadyRk5Solver(
        const mesh& m_in,
        const std::vector<inputVar>& q0, 
        const constants cin
    ) {
        tolerance = 1.e-4;
        c = cin;
        q.resize(q0.size());
        for (int i=0; i<q0.size(); ++i) {
            q[i] = convertVariable(q0[i], c);
        }
        qt.resize(q.size());
        dtScale.resize(q.size());
        dtv.resize(q.size());

        k.resize(q.size());
        
        qk.resize(q.size());

        set_mesh(m_in);

        solver_type = "Rk5";
    }

    void calc_qt() {

        for (int i=0; i<qk.size(); ++i) {
            qk[i] = q[i];
        }

        for (int i=0; i<(alphas.size() - 1); ++i) {
            internal_calc_qt(k, qk);

            for (int j=0; j<qk.size(); ++j) {
                qk[j] = q[j] + k[j] * dtv[j] * alphas[i];
            }
        }
        internal_calc_qt(k, qk);

        for (int i=0; i<qt.size(); ++i) {
            qt[i] = k[i] * alphas[alphas.size()-1] ;
        }
    }

    bool isTransient() {return false;}

    void dtFromCfl() {
        // Reset dt scale
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
            dtv[i] = cfl * m.facesAreas[i] / dtScale[i];
        }
    }

    void update_cells() {
        for (int i=0; i<m.facesAreas.size(); ++i) {
            q[i] += qt[i] * dtv[i];
        }
    }

};






#endif

