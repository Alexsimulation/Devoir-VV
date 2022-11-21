#ifndef EULER_SOLVER_H
#define EULER_SOLVER_H

#include <euler/classes.h>
#include <euler/flux.h>
#include <euler/limiter.h>
#include <map>
#include <chrono>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <csignal>


#ifdef EULER_PYTHON_MODULE
    #include <pybind11/embed.h>

    void catch_signals() {
        auto handler = [](int code) { throw std::runtime_error("SIGNAL " + std::to_string(code)); };
        signal(SIGINT, handler);
        signal(SIGTERM, handler);
        signal(SIGKILL, handler);
    }
#else
    void catch_signals() {}
#endif


// Format string to scientific notation
std::string formatScientific(real v) {
    std::ostringstream streamObj;
    streamObj << std::scientific;
    streamObj << v;
    std::string s = streamObj.str();
    return s;
}

// Format string to scientific notation
std::string formatDefault(real v) {
    std::ostringstream streamObj;
    streamObj << v;
    std::string s = streamObj.str();
    return s;
}

void emptySmoother(
    std::vector<var>& qt_, 
    std::vector<var>& smoother, 
    std::vector<var>& smoother_qt, 
    const mesh& m,
    const real& epsilon,
    const uint& iters
) {}

void explicitSmoother(
    std::vector<var>& qt_, 
    std::vector<var>& smoother, 
    std::vector<var>& smoother_qt, 
    const mesh& m,
    const real& epsilon,
    const uint& iters
) {
    // Explicit residual smoothing
    for (int iter=0; iter<iters; ++iter) {
        for (int i=0; i<m.edgesFaces0.size(); ++i) {
            auto f0 = m.edgesFaces0[i];
            auto f1 = m.edgesFaces1[i];
            const auto v = (qt_[f1] - qt_[f0]) * epsilon;
            qt_[f0] += v;
            qt_[f1] -= v;
        }
    }
}

void implicitSmoother(
    std::vector<var>& qt_, 
    std::vector<var>& smoother, 
    std::vector<var>& smoother_qt, 
    const mesh& m,
    const real& epsilon,
    const uint& iters
) {
    for (int i=0; i<smoother_qt.size(); ++i) {
        smoother_qt[i] = qt_[i];
    }
    // Implicit residual smoothing
    // Jacobi diagonal dominant solver
    for (int jacobi=0; jacobi<iters; ++jacobi) {
        for (int i=0; i<smoother.size(); ++i) {
            smoother[i] = 0.;
        }
        for (int i=0; i<m.edgesFaces0.size(); ++i) {
            auto f0 = m.edgesFaces0[i];
            auto f1 = m.edgesFaces1[i];
            smoother[f0] += qt_[f1] * -1. * epsilon;
            smoother[f1] += qt_[f0] * -1. * epsilon;
        }
        for (int i=0; i<qt_.size(); ++i) {
            real ne;
            if (m.faceIsTriangle[i]) {
                ne = 3.;
            } else {
                ne = 4.;
            }
            qt_[i] = (smoother_qt[i] - smoother[i])/(1. + ne * epsilon);
        }
    }
}


struct solutionInfo {
    real final_time;
    int solver_steps;
    real solver_time;
    var final_residuals;
};



// General solver definition
class solver {
public:

    std::string solver_type = "root";

    mesh m;
    std::vector<var> q;
    std::vector<var> qt;
    std::vector<varVec2> g;
    std::vector<var> q_min;
    std::vector<var> q_max;
    std::vector<var> limiters;
    std::vector<var> q_nodes;
    std::vector<real> dtScale;
    std::vector<var> edgesFlux;
    std::vector<var> boundsFlux;
    std::vector<var> residuals;

    var (*limiterFunc)(const varVec2&, const var&, const var&, const var&, const vec2&, const real&, const vec2&) = &limiterVenkatakrishnan;
    void (*smootherFunc)(std::vector<var>&, std::vector<var>&, std::vector<var>&, const mesh&, const real&, const uint&) = &emptySmoother;
    std::string limiter_type = "venkatakrishnan";
    std::string smoother_type = "none";

    constants c;
    std::map<std::string, std::string> bcTypes;
    std::map<std::string, var> bcValues;

    uint printInterval = 1;
    uint maxSteps = 4294967295;
    real cfl = 0.8;
    real maxTime = 1e10;
    uint rescaleNSteps = 1;

    real dt;
    std::vector<real> dtv;

    bool isSecondOrder = false;
    real smoother_coefficient = 0.3;
    uint smoother_iters = 1;
    std::vector<var> smoother;
    std::vector<var> smoother_qt;

    solutionInfo info;

    bool saveCoefficients = false;
    std::string coeffField;
    inputVar coeffVar;
    std::vector<real> CD, CL, CM;

    real tolerance = 1.e-14;

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> stop;
    int step = 0;
    real time = 0.;
    var Runscaled = var(1., 1., 1., 1.);
    var R = var(1., 1., 1., 1.);
    var Rn = var(0., 0., 0., 0.);

    virtual bool isTransient() {}

    void set_tolerance(real in) { tolerance = in; }

    void set_smoother(std::string s_type, real coeff = -1., uint iters = 2) {
        if (s_type == "explicit") {
            smootherFunc = &explicitSmoother;
        } else if (s_type == "implicit") {
            smootherFunc = &implicitSmoother;
            smoother.resize(q.size());
            smoother_qt.resize(q.size());
        } else {
            throw std::invalid_argument( "unknown smoother type " + s_type );
        }
        if (coeff < 0) {
            if (s_type == "explicit") {coeff = 0.1;}
            else if (s_type == "implicit") {coeff = 0.6;}
        }
        smoother_coefficient = coeff;
        smoother_iters = iters;
        smoother_type = s_type;
    }
    void set_limiter(std::string s_type) {
        if (s_type == "van-leer") {
            limiterFunc = &limiterVanLeer;
        } else if (s_type == "venkatakrishnan") {
            limiterFunc = &limiterVenkatakrishnan;
        } else if (s_type == "van-albada") {
            limiterFunc = &limiterVanAlbada;
        } else if (s_type == "superbee") {
            limiterFunc = &limiterSuperBee;
        } else if (s_type == "overbee") {
            limiterFunc = &limiterOverBee;
        } else if (s_type == "min-mod") {
            limiterFunc = &limiterMinMod;
        } else if (s_type == "zero") {
            limiterFunc = &limiterZero;
        } else {
            throw std::invalid_argument( "unknown limiter type " + s_type );
        }
        limiter_type = s_type;
    }

    std::vector<var> get_limiters() {
        return limiters;
    }

    void set_bc(std::string bound, std::string b_type, inputVar b_value = inputVar()) {
        bcTypes[bound] = b_type;
        bcValues[bound] = convertVariable(b_value, c);
        for (int i=0; i<m.boundsLengths.size(); ++i) {
            if (m.boundsTags[i] == bound) {
                m.boundsFuncs[i] = boundaryFuncMap().at(b_type);
            }
        }
    }

    void writeVtk(const std::string filename);

    void set_mesh(const mesh& m_in) {
        m = m_in;
        edgesFlux.resize(m.edgesLengths.size());
        boundsFlux.resize(m.boundsLengths.size());
    }
    void set_cfl(real cfl_in) { cfl = cfl_in; }
    void set_print_interval(uint in) { printInterval = in; }
    void set_max_steps(uint in) { maxSteps = in; }
    uint get_max_steps() {return maxSteps;}
    void set_max_time(real in) { maxTime = in; }
    void set_order(int order) {
        if (order == 1) {
            isSecondOrder = false;
        } else if (order == 2) {
            isSecondOrder = true;
            if (g.size() != q.size()) {
                g.resize(q.size());
                q_nodes.resize(m.nodes.size());
                q_min.resize(q.size());
                q_max.resize(q.size());
                limiters.resize(q.size());
            }
        } else {
            throw std::invalid_argument( "order " + std::to_string(order) + " not supported");
        }
    }
    void set_coefficient_save_field(std::string s, inputVar v) {
        saveCoefficients = true;
        coeffField = s;
        coeffVar = v;
    }

    var computeResidual() {
        var R; R = 0.;
        for (int i=0; i<qt.size(); ++i) {
            R += qt[i] * qt[i];
        }
        return sqrt(R);
    }

    vec2 forceOnField(std::string field);

    std::tuple<real, real, real> computeCoefficients(std::string field, inputVar qinf_);

    void writeResiduals(std::string filename);

    std::vector<var> get_residuals() {
        return residuals;
    }

    virtual void dtFromCfl() {}

    virtual void calc_qt() {}

    virtual void update_cells() {}


    void internal_calc_qt(std::vector<var>& qt_, std::vector<var>& q_) {
        // Reset to 0
        #pragma omp simd
        for (int i=0; i<qt_.size(); ++i) {
            qt_[i] = 0.;
        }

        // If second order, compute gradients
        if (isSecondOrder) {
            computeNodalValues(q_);
            computeGradients(q_);
            compute_q_min_max(q_);
            compute_limiters(q_);
        }

        #pragma omp parallel default(none) shared(m, edgesFlux, boundsFlux, q_, g, limiters)
        {
            
            // Compute edges flux
            #pragma omp for
            for (int i=0; i<m.edgesLengths.size(); ++i) {
                
                // Add edge contribution to qt
                const auto fL = m.edgesFaces0[i];
                const auto fR = m.edgesFaces1[i];
                const auto n = m.edgesNormals[i];
                auto qfL = q_[fL];
                auto qfR = q_[fR];
                if (isSecondOrder) {
                    qfL += dot(
                            g[fL] * limiters[fL],
                            m.edgesCenters[i] - m.facesCenters[fL]
                        );

                    qfR += dot(
                            g[fR] * limiters[fR],
                            m.edgesCenters[i] - m.facesCenters[fR]
                        );
                }

                // Compute central flux
                auto fluxR = convectiveFlux(qfR, n, c);
                auto fluxL = convectiveFlux(qfL, n, c);
                auto flux = (fluxR + fluxL)*0.5;
                
                //auto flux = convectiveFlux((qfL + qfR)*0.5, n, c);

                // Compute diffusion
                auto diffusion = artificialDiffusion(qfR, qfL, n, c);
                // Add artificial diffusion
                flux -= diffusion * 0.5;
                // Substract flux to fR, add to fL

                edgesFlux[i] = flux * m.edgesLengths[i];
            }

            // Compute boundary fluxes
            // Boundaries are currently first order
            #pragma omp for
            for (int i=0; i<m.boundsLengths.size(); ++i) {

                // Add bound contribution to qt
                const auto tag = m.boundsTags[i];
                const auto n = m.boundsNormals[i];

                auto qf = q_[m.boundsFaces[i]];

                auto df = m.boundsCenters[i] - m.facesCenters[m.boundsFaces[i]];

                auto gf = varVec2(0.);
                if (isSecondOrder) {
                    gf = g[m.boundsFaces[i]] * limiters[m.boundsFaces[i]];
                }

                auto bvar = m.boundsFuncs[i](
                    qf, bcValues.at(tag), 
                    gf, 
                    df, n, c
                );
                
                if (isSecondOrder) {
                    qf += dot(
                            g[m.boundsFaces[i]] * limiters[m.boundsFaces[i]],
                            m.boundsCenters[i] - m.facesCenters[m.boundsFaces[i]]
                        );
                }

                auto flux = convectiveFlux((qf + bvar)*0.5, n, c);

                // Compute diffusion
                auto diffusion = artificialDiffusion(bvar, qf, n, c);
                // Add artificial diffusion
                flux -= diffusion * 0.5;

                // Substract flux to f
                boundsFlux[i] = flux * m.boundsLengths[i];
            }
        }   // Parallel code

        // Add fluxes to face variations
        for (int i=0; i<m.edgesLengths.size(); ++i) {
            const auto fL = m.edgesFaces0[i];
            const auto fR = m.edgesFaces1[i];
            qt_[fL] -= edgesFlux[i];
            qt_[fR] += edgesFlux[i];
        }

        // Add boundary fluxes to face variations
        for (int i=0; i<m.boundsLengths.size(); ++i) {
            const auto f = m.boundsFaces[i];
            qt_[f] -= boundsFlux[i];
        }

        for (int i=0; i<m.facesAreas.size(); ++i) {
            // Normalize by area
            qt_[i] /= m.facesAreas[i];
        }

        // At the end of all of this, smooth the residuals
        smootherFunc(qt_, smoother, smoother_qt, m, smoother_coefficient, smoother_iters);
    }

    bool convergence_check() {
        return (step < maxSteps)&(time < maxTime)&(std::max(R) > tolerance);
    }

    void init_solver(bool verb) {
        start = std::chrono::high_resolution_clock::now();
        step = 0;
        time = 0.;
        Runscaled = var(1., 1., 1., 1.);
        R = var(1., 1., 1., 1.);
        Rn = var(0., 0., 0., 0.);
        residuals.clear();
        CD.clear();
        CL.clear();
        CM.clear();

        if (verb) {
            std::cout << "Solver Settings " << "\n";
            std::cout << " - CFL       = " << cfl << "\n";
            std::cout << " - Max steps = " << maxSteps << "\n";
            std::cout << " - Tolerance = " << tolerance << "\n";
            std::cout << std::endl;
        }
    }

    void internal_solver_step(
        bool verb
    ) {
        // Compute dt from cfl
        dtFromCfl();
        if (isTransient()) {
            if ((time + dt) > maxTime) {
                dt = maxTime - time;
            }
        }

        // Compute qt
        calc_qt();

        // Update all cells
        update_cells();

        // Save coefficients
        if (saveCoefficients) {
            real cd, cl, cm;
            std::tie(cd, cl, cm) = computeCoefficients(
                coeffField, 
                coeffVar
            );
            CD.push_back(cd);
            CL.push_back(cl);
            CM.push_back(cm);
        }

        // Update time, step
        if (isTransient()) {
            time += dt;
        }
        if (step < rescaleNSteps) {
            Runscaled = computeResidual();
            Rn += Runscaled / ((real) rescaleNSteps);
            R = Runscaled / Rn;
            residuals.push_back(R);
        } else if (!isTransient()) {
            R = computeResidual() / Rn;
            residuals.push_back(R);
        }
        step += 1;
        // Print results
        if (verb) {
            if ((step % printInterval) == 0) {
                if (isTransient()) {
                    R = computeResidual() / Rn;
                }
                std::cout << "Step = " << step;
                if (isTransient()) {
                    std::cout << ", Time = " << time;
                }
                std::cout << "\n";
                std::cout << " - R(rho)   = " << R.rho  << "\n";
                std::cout << " - R(rho*u) = " << R.rhou << "\n";
                std::cout << " - R(rho*v) = " << R.rhov << "\n";
                std::cout << " - R(rho*e) = " << R.rhoe << "\n";
                std::cout << std::endl;
            }
        }
        
        catch_signals();
    }

    void end_solver(bool verb) {
        stop = std::chrono::high_resolution_clock::now();
        if (verb) {
            if (isTransient()) {
                R = computeResidual() / Rn;
            }
            std::cout << "Step = " << step;
            if (isTransient()) {
                std::cout << ", Time = " << time;
            }
            std::cout << "\n";
            std::cout << " - R(rho)   = " << R.rho  << "\n";
            std::cout << " - R(rho*u) = " << R.rhou << "\n";
            std::cout << " - R(rho*v) = " << R.rhov << "\n";
            std::cout << " - R(rho*e) = " << R.rhoe << "\n";
            std::cout << std::endl;
        }
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        if (verb) {
            if (isTransient()) {
                std::cout << "Final time    = " << time << std::endl;
            }
            std::cout <<     "Solution time = " << ((double) duration.count())/1000. << " s" << std::endl;
        }
        info.solver_steps = step;
        info.solver_time = ((double) duration.count())/1000.;
        info.final_residuals = R;
        info.final_time = time;
    }

    void simulate(bool verb=true) {
        // Init solver
        init_solver(verb);

        // Main loop
        while (convergence_check()) {
            internal_solver_step(verb);
        }
        
        end_solver(verb);
    }

    void computeNodalValues(std::vector<var>& q_) {
        // Reset q_nodes
        for (int i=0; i<q_nodes.size(); ++i) {
            q_nodes[i] = 0.;
        }
        // Loop over all faces and add nodal values
        for (int i=0; i<m.facesAreas.size(); ++i) {
            std::array<uint, 3> nodes = {
                m.facesNodes0[i],
                m.facesNodes1[i],
                m.facesNodes2[i],
            };
            for (auto node : nodes) {
                q_nodes[node] += q_[i] * 1./norm(m.facesCenters[i] - m.nodes[node]);
            }
            if (!m.faceIsTriangle[i]) {
                uint node = m.facesNodes3[i];
                q_nodes[node] += q_[i] * 1./norm(m.facesCenters[i] - m.nodes[node]);
            }
        }
        for (int i=0; i<q_nodes.size(); ++i) {
            q_nodes[i] /= m.nodesGradWeights[i];
        }
    }

    void computeGradients(std::vector<var>& q_) {
        // Reset gradients
        for (int i=0; i<g.size(); ++i) {
            g[i].x = 0.;
            g[i].y = 0.;
        }
        // Compute gradient using divergence theorem
        for (int i=0; i<m.edgesLengths.size(); ++i) {
            // Add edge contribution to faces gradients
            const auto n0 = m.edgesNodes0[i];
            const auto n1 = m.edgesNodes1[i];
            const auto edgeV = varVec2( (q_nodes[n0] + q_nodes[n1])*0.5 * m.edgesLengths[i] );
            //const auto f0 = m.edgesFaces0[i];
            //const auto f1 = m.edgesFaces1[i];
            //const auto edgeV = varVec2( (q_[f0] + q_[f1])*0.5 * m.edgesLengths[i] );
            g[m.edgesFaces0[i]] += edgeV * m.edgesNormals[i];
            g[m.edgesFaces1[i]] -= edgeV * m.edgesNormals[i];
        }
        for (int i=0; i<m.boundsLengths.size(); ++i) {
            // Add boundary contribution to faces gradients
            const auto n0 = m.boundsNodes0[i];
            const auto n1 = m.boundsNodes1[i];
            const auto edgeV = varVec2( (q_nodes[n0] + q_nodes[n1])*0.5 * m.boundsLengths[i] );
            //const auto f0 = m.boundsFaces[i];
            //const auto f1 = m.boundsFaces[i];
            //const auto edgeV = varVec2( (q_[f0] + q_[f1])*0.5 * m.edgesLengths[i] );
            g[m.boundsFaces[i]] += edgeV * m.boundsNormals[i];
        }
        // Divide by cell areas
        for (int i=0; i<m.facesAreas.size(); ++i) {
            g[i] /= m.facesAreas[i];
        }
    }

    void compute_q_min_max(std::vector<var>& q_) {
        for (int i=0; i<q_min.size(); ++i) {
            q_min[i] = q_[i];
        }
        for (int i=0; i<q_max.size(); ++i) {
            q_max[i] = q_[i];
        }
        for (int i=0; i<m.edgesCenters.size(); ++i) {
            const auto f0 = m.edgesFaces0[i];
            const auto f1 = m.edgesFaces1[i];
            q_min[f0] = std::min(q_[f1], q_min[f0]);
            q_max[f0] = std::max(q_[f1], q_max[f0]);
            q_min[f1] = std::min(q_[f0], q_min[f1]);
            q_max[f1] = std::max(q_[f0], q_max[f1]);
        }
        for (int i=0; i<m.boundsCenters.size(); ++i) {
            const auto f = m.boundsFaces[i];
            const auto tag = m.boundsTags[i];
            const auto n = m.boundsNormals[i];

            // Compute virtual cell boundary value
            //   Use order 1 extrapolation for wall cell value
            /*
            const auto bvar = m.boundsFuncs[i](
                q_[f],
                bcValues.at(tag), 
                varVec2(0.),
                vec2(0.), n, c
            );
            */
            const auto bvar = q_[f];// + dot(g[f], (m.boundsCenters[i] - m.facesCenters[f]));
            q_min[f] = std::min(bvar, q_min[f]);
            q_max[f] = std::max(bvar, q_max[f]);
        }
        
    }

    void compute_limiters(const std::vector<var>& q_) {
        // Minimums and maximums have been computed
        // Init limiters to really high value
        for (int i=0; i<limiters.size(); ++i) {
            limiters[i] = 1.e200;
        }
        // Compute cell limiters by looping on edges
        for (int i=0; i<m.edgesCenters.size(); ++i) {
            const auto f0 = m.edgesFaces0[i];
            const auto f1 = m.edgesFaces1[i];
            const auto l0 = limiterFunc(
                g[f0], 
                q_min[f0], 
                q_max[f0], 
                q_[f0], 
                m.edgesCenters[i] - m.facesCenters[f0], 
                m.facesAreas[f0], 
                m.edgesNormals[i]
            );
            const auto l1 = limiterFunc(
                g[f1], 
                q_min[f1], 
                q_max[f1], 
                q_[f1], 
                m.edgesCenters[i] - m.facesCenters[f1], 
                m.facesAreas[f1], 
                m.edgesNormals[i]*-1.
            );
            limiters[f0] = std::min(limiters[f0], l0);
            limiters[f1] = std::min(limiters[f1], l1);
        }
        // Also compute by looping on boundaries
        for (int i=0; i<m.boundsCenters.size(); ++i) {
            const auto f = m.boundsFaces[i];
            const auto l = limiterFunc(
                g[f], 
                q_min[f], 
                q_max[f], 
                q_[f], 
                m.boundsCenters[i] - m.facesCenters[f], 
                m.facesAreas[f], 
                m.boundsNormals[i]
            );
            limiters[f] = std::min(limiters[f], l);
        }
    }

    std::vector<var> getField(std::string field);

    void writeField(std::string field, std::string filename, inputVar q_inf_);

    void writeSavedCoefficients(std::string filename);

    std::vector<inputVar> get_solution();

    std::vector<var>& get_raw_solution() { return q; }

    std::string log();

};




#endif