#ifndef EULER_POST_H
#define EULER_POST_H


#include <euler/solver.h>



#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>




std::vector<var> solver::getField(std::string field) {
    std::vector<var> out;
    for (int i=0; i<m.boundsLengths.size(); ++i) {
        if (m.boundsTags[i] == field) {
            out.push_back(q[m.boundsFaces[i]]);
        }
    }
    return out;
}

void solver::writeField(std::string field, std::string filename, inputVar q_inf_) {
    var q_inf = convertVariable(q_inf_, c);

    real c_inf = calc_c(q_inf, c.gamma);
    real u_inf = q_inf.rhou/q_inf.rho;
    real v_inf = q_inf.rhov/q_inf.rho;
    real mach_inf = sqrt(u_inf*u_inf + v_inf*v_inf) / c_inf;
    real p_inf = calc_p(q_inf, c.gamma);

    std::string s = "x y rho U M Cp\n";
    for (int i=0; i<m.boundsLengths.size(); ++i) {
        if (m.boundsTags[i] == field) {
            uint face = m.boundsFaces[i];
            real rho = q[face].rho;
            real u = q[face].rhou/rho;
            real v = q[face].rhov/rho;
            real p = calc_p(q[face], c.gamma);
            real umag = sqrt(u*u + v*v);
            real cs = calc_c(q[face], c.gamma);
            real mach = umag/cs;
            real cp = 2./(c.gamma*mach_inf*mach_inf)*(p/p_inf - 1.);

            s += std::to_string(m.boundsCenters[i].x) + " ";
            s += std::to_string(m.boundsCenters[i].y) + " ";
            s += std::to_string(rho) + " ";
            s += std::to_string(umag) + " ";
            s += std::to_string(mach) + " ";
            s += std::to_string(cp) + "\n";
        }
    }
    std::ofstream out(filename);
    out << s;
    out.close();
}

std::vector<inputVar> solver::get_solution() {
    std::vector<inputVar> out(q.size());
    for (int i=0; i<out.size(); ++i) {
        out[i] = convertVariable(q[i], c);
    }
    return out;
}


void solver::writeResiduals(std::string filename) {
    std::string s = "Iteration rho rho*u rho*v rho*e\n";
    for (int i=0; i<residuals.size(); ++i) {
        s += std::to_string(i) + " ";
        s += formatScientific(residuals[i].rho) + " ";
        s += formatScientific(residuals[i].rhou) + " ";
        s += formatScientific(residuals[i].rhov) + " ";
        s += formatScientific(residuals[i].rhoe) + "\n";
    }
    std::ofstream out(filename);
    out << s;
    out.close();
}


vec2 solver::forceOnField(std::string field) {
    vec2 F = vec2(0., 0.);
    for (int i=0; i<m.boundsTags.size(); ++i) {
        if (field == m.boundsTags[i]) {
            const real p = calc_p(q[m.boundsFaces[i]], c.gamma);
            F += vec2(p*m.boundsNormals[i].x, p*m.boundsNormals[i].y) * m.boundsLengths[i];
        }
    }
    return F;
}


std::tuple<real, real, real> solver::computeCoefficients(std::string field, inputVar q_inf_) {
    
    var q_inf = convertVariable(q_inf_, c);
    // Calcul des coefficients aerodynamiques

    real c_inf = calc_c(q_inf, c.gamma);
    real u_inf = q_inf.rhou/q_inf.rho;
    real v_inf = q_inf.rhov/q_inf.rho;
    real mach_inf = sqrt(u_inf*u_inf + v_inf*v_inf) / c_inf;
    real p_inf = calc_p(q_inf, c.gamma);

    // Compute angle of attack
    real aoa = std::atan(v_inf/u_inf);

    // Compute chord length
    bool foundField = false;
    real minX, maxX;
    real meanY = 0.;
    int Np = 0;
    for (int i=0; i<m.boundsTags.size(); ++i) {
        if (field == m.boundsTags[i]) {
            if (!foundField) {
                foundField = true;
                minX = m.boundsCenters[i].x;
                maxX = m.boundsCenters[i].x;
            } else if (m.boundsCenters[i].x < minX) {
                minX = m.boundsCenters[i].x;
            } else if (m.boundsCenters[i].x > maxX) {
                maxX = m.boundsCenters[i].x;
            }
            meanY += m.boundsCenters[i].y;
            Np += 1;
        }
    }
    real S = maxX - minX;
    meanY /= ((real) Np);

    // Compute Forces
    vec2 M_pos = vec2(S*0.25, meanY);    // Center of moment at quarter chord
    real M = 0;
    vec2 F = vec2(0., 0.);
    vec2 Fi = vec2(0., 0.);

    for (int i=0; i<m.boundsTags.size(); ++i) {
        if (field == m.boundsTags[i]) {

            uint face = m.boundsFaces[i];
            real rho = q[face].rho;
            real u = q[face].rhou/rho;
            real v = q[face].rhov/rho;
            real p = calc_p(q[face], c.gamma);
            real umag = sqrt(u*u + v*v);
            real cs = calc_c(q[face], c.gamma);
            real mach = umag/cs;
            real cp = 2./(c.gamma*mach_inf*mach_inf)*(p/p_inf - 1.);

            Fi = m.boundsNormals[i] * cp * m.boundsLengths[i];
            F += Fi;
            M -= cross(m.boundsCenters[i] - M_pos, Fi);
        }
    }

    // Compute and return coefficients
    real CD = F.x*std::cos(aoa) + F.y*std::sin(aoa);
    real CL = -F.x*std::sin(aoa) + F.y*std::cos(aoa);
    real CM = M;
    return std::make_tuple(CD, CL, CM);
}


void solver::writeSavedCoefficients(std::string filename) {
    std::string s = "Iteration cd cl cm\n";
    for (int i=0; i<CD.size(); ++i) {
        s += std::to_string(i) + " ";
        s += formatScientific(CD[i]) + " ";
        s += formatScientific(CL[i]) + " ";
        s += formatScientific(CM[i]) + "\n";
    }
    std::ofstream out(filename);
    out << s;
    out.close();
}


std::string solver::log() {
    // Output solver log

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y %Hh-%Mm-%Ss");
    auto current_time = oss.str();

    std::string s = "";
    s += "Euler Solver\n";
    s += " - Log created at " + current_time + "\n";
    s += "\n";

    s += "Mesh metrics\n";
	s += " - Number of nodes = " + std::to_string(m.nodes.size()) + "\n";
	s +=  " - Number of edges = " + std::to_string(m.edgesNormals.size()) + "\n";
	s += " - Number of faces = " + std::to_string(m.facesAreas.size()) + "\n";
    s += " - Import time     = " + formatDefault(m.import_duration) + " s\n";
	s += "\n";

    // Solver settings
    s += "Solver Settings \n";
    s += " - CFL       = " + formatDefault(cfl) + "\n";
    s += " - Max steps = " + std::to_string(maxSteps) + "\n";
    if (isTransient()) {
        s += " - Max time  = " + formatDefault(maxTime) + "\n";
    }
    s += " - Tolerance = " + formatDefault(tolerance) + "\n";
    s += " - Solver    = " + solver_type + " ";
    if (isTransient()) {
        s += "transient\n";
    } else {
        s += "steady\n";
    }
    s += " - Order     = ";
    s += (isSecondOrder ? "2" : "1");
    s += "\n";
    if (isSecondOrder) {
        s += " - Limiter   = " + limiter_type + "\n";
    }
    s += " - Smoother\n   - Smoother Type = " + smoother_type + "\n";
    if (smoother_type != "none") {
        s += "   - Iterations    = " + formatDefault(smoother_iters) + "\n";
        s += "   - Coefficient   = " + formatDefault(smoother_coefficient) + "\n";
    } else { s += "\n"; }
    s += "\n";

    // Solution metrics
    s += "Solution metrics \n";
    s += " - Solver time   = " + formatDefault(info.solver_time) + " s\n";
    s += " - Solver steps  = " + std::to_string(info.solver_steps) + "\n";
    if (isTransient()) {
        s += " - Physical time = " + formatDefault(info.final_time) + "\n";
    }
    s += " - Final residuals\n";
    auto R = info.final_residuals;
    s += "   - R(rho)   = " + formatDefault(R.rho) + "\n";
    s += "   - R(rho*u) = " + formatDefault(R.rhou) + "\n";
    s += "   - R(rho*v) = " + formatDefault(R.rhov) + "\n";
    s += "   - R(rho*e) = " + formatDefault(R.rhoe) + "\n";

    return s;
}



#endif
