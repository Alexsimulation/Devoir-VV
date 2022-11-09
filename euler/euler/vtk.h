#ifndef EULER_VTK_H
#define EULER_VTK_H


#include <euler/classes.h>
#include <euler/mesh.h>
#include <euler/solver.h>

#include <fstream>


void solver::writeVtk(
    const std::string filename
) {
    // Write solution q to .vtk file filename
    std::string s = "";

    s += "# vtk DataFile Version 2.0\n";
    s += filename + " \n";
    s += "ASCII\n";
    s += "DATASET UNSTRUCTURED_GRID\n";
    s += "POINTS " + std::to_string(m.nodes.size()) + " float\n";
    for (auto ni : m.nodes) {
        s += std::to_string(ni.x) + " " + std::to_string(ni.y) + " 0\n";
    }

    // Figure out the total number of integer tags in the cells
    int n_face_ints = 0;
    for (auto isTri : m.faceIsTriangle) {
        if (isTri) {
            n_face_ints += 4;
        } else {
            n_face_ints += 5;
        }
    }
    s += "CELLS " + std::to_string(m.facesAreas.size()) + " ";
    s += std::to_string(n_face_ints) + "\n";

    // Add each face
    for (int fi=0; fi<m.facesAreas.size(); ++fi) {

        std::vector<uint> ns = {
            m.facesNodes0[fi],
            m.facesNodes1[fi],
            m.facesNodes2[fi]
        };
        if (!m.faceIsTriangle[fi]) {
            s += "4";
            ns.push_back(m.facesNodes3[fi]);
        } else {
            s += "3";
        }
        for (auto ni : ns) {
            s += " " + std::to_string(ni);
        }
        s += "\n";
    }

    // Write cell types
    s += "CELL_TYPES " + std::to_string(m.facesAreas.size()) + "\n";
    for (auto isTriangle : m.faceIsTriangle) {
        if (isTriangle) {
            // Tri face
            s += "5\n";
        } else {
            // Quad face
            s += "9\n";
        }
    }
    // Save data
    s += "CELL_DATA " + std::to_string(m.facesAreas.size()) + "\n";
    // Save velocity
    s += "VECTORS U float\n";
    for (int i=0; i<m.facesAreas.size(); ++i) {
        s += std::to_string(q[i].rhou/q[i].rho) + " ";
        s += std::to_string(q[i].rhov/q[i].rho) + " 0\n";
    }

    // Save density
    s += "SCALARS rho float 1\n";
    s += "LOOKUP_TABLE default\n";
    for (int i=0; i<q.size(); ++i) {
        s += std::to_string(q[i].rho) + "\n";
    }

    // Save pressure
    s += "SCALARS p float 1\n";
    s += "LOOKUP_TABLE default\n";
    for (int i=0; i<q.size(); ++i) {
        s += std::to_string(calc_p(q[i], c.gamma)) + "\n";
    }

    // Save mach number
    s += "SCALARS M float 1\n";
    s += "LOOKUP_TABLE default\n";
    for (int i=0; i<q.size(); ++i) {
        s += std::to_string(
                sqrt(q[i].rhou/q[i].rho*q[i].rhou/q[i].rho + q[i].rhov/q[i].rho*q[i].rhov/q[i].rho) /
                calc_c(q[i], c.gamma))
             + "\n";
    }

    // Save data to file
    std::ofstream out(filename);
    out << s;
    out.close();

}



#endif
