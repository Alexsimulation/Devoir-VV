#ifndef EULER_GUI_PLOT_H
#define EULER_GUI_PLOT_H



void draw_axis(std::vector<double> x, std::vector<double> y) {



}


void ExampleWindow::gen_plot_gl_arrays() {
    // Populate vertices array with plot info

    // Get residuals
    auto res = solver.residuals;
    uint n_lines_added = 0;

    if (res.size() > 4) {

        // Get the number of cells in mesh
        uint n_lines = res.size()*4;
        uint n_vertex = n_lines * 2;
        uint n_tags = n_vertex * 6; // Color tags are added

        VERTEX_NUMBER_L = n_vertex;
        if (vertices_l.size() != solver.get_max_steps()) {
            vertices_l.resize(solver.get_max_steps());
        }

        for (uint i=0; i<res.size()-1; ++i) {
            
            for (uint j=0; j<4; ++j) {

                std::array<vec2, 2> ns = {
                    vec2(((double) i)/((float) res.size()-1), log(res[i][j])/4. + 1.),
                    vec2(((double) i+1)/((float) res.size()-1), log(res[i+1][j])/4. + 1.)
                };
                
                for (uint jj=0; jj<ns.size(); ++jj) {
                    uint k = n_lines_added*6*2 + jj*6;   // This vertex start position
                    
                    // Position
                    vertices_l[k] = 2.*ns[jj].x - 1.;
                    vertices_l[k+1] = 2.*ns[jj].y - 1.;
                    vertices_l[k+2] = 0.;

                    // Color
                    vertices_l[k+3] = 1.0f;
                    vertices_l[k+4] = 1.0f;
                    vertices_l[k+5] = 1.0f;
                }

                n_lines_added += 1;
            }
        }
    }

    std::cout << "generated plot arrays, size = " << n_lines_added << std::endl;


    remake_buffers();
}



#endif