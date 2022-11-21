#ifndef EULER_GUI_RENDER_H
#define EULER_GUI_RENDER_H


#include <euler/euler.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include <gtkmm.h>

#define GLEW_STATIC
#include <GL/glew.h>

#include <gui/core.h>

// Include shaders
#include <gui/shaderFragment>
#include <gui/shaderVertex>
#include <tuple>



uint VERTEX_NUMBER = 0;
std::vector<GLfloat> vertices;

uint VERTEX_NUMBER_Q = 0;
std::vector<GLfloat> vertices_q;

uint VERTEX_NUMBER_L = 0;
std::vector<GLfloat> vertices_l;

GLfloat transform[] = {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f,
        };

bool DRAW_POLY_AS_LINES;




std::array<float, 3> cmap_soft_jet(const float& t) {
    // Jet colormap
    std::array<float, 3> c;
    if (t < 0.25) {
        c[0] = 0.;
        c[1] = 4.*t;
        c[2] = 1.;
    } else if (t < 0.5) {
        c[0] = 0.;
        c[1] = 1.;
        c[2] = 1. - (t - 0.25)/0.25;
    } else if (t < 0.75) {
        c[0] = 1. + (t - 0.75)/0.25;
        c[1] = 1.;
        c[2] = 0.;
    } else {
        c[0] = 1.;
        c[1] = 1. - (t - 0.75)/0.25;
        c[2] = 0.;
    }
    return c;
}


std::array<float, 3> cmap_rainbow(const float& t) {
    // Jet colormap
    std::array<float, 3> c;
    c[1] = 1. - 2.*(t - 0.5)*(t - 0.5);
    c[2] = 1. - t * t;
    if (t < 0.25) {
        c[0] = 0.5 - 2.*t;
    } else if (t < 0.75) {
        c[0] = 2.*(t - 0.25);
    } else {
        c[0] = 1.;
    }
    return c;
}

std::array<float, 3> cmap_inferno(const float& t) {
    // Jet colormap
    std::array<float, 3> c;
    c[0] = (2*t - t*t);
    c[1] = t*t;
    c[2] = (2*t - 2.6*t*t + 1.1*std::pow(t, 10));
    return c;
}

std::array<float, 3> cmap_viridis(const float& t) {
    // Jet colormap
    std::array<float, 3> c;
    c[0] = 0.3 - 1.3*t*t + 2.*t*t*t*t;
    c[1] = 0.9 - 0.9*(t-1.)*(t-1.);
    c[2] = 0.6 - 1.35*(t-0.35)*(t-0.35);
    return c;
}

std::vector<std::array<float, 3> (*)(const float&)> cmap_map = {
    &cmap_rainbow, &cmap_inferno, &cmap_soft_jet,
    &cmap_viridis
    };
std::map<std::string, uint> cmap_ref = {
    {"rainbow", 0}, {"inferno", 1}, {"soft jet", 2},
    {"viridis", 3}
};

float get_rho_fun(const var& q) { return q.rho; }
float get_u_fun(const var& q) { return q.rhou/q.rho; }
float get_v_fun(const var& q) { return q.rhov/q.rho; }
float get_p_fun(const var& q) { return (1.4 - 1.)*(q.rhoe - 0.5/q.rho*(q.rhou*q.rhou + q.rhov*q.rhov)); }
float get_u_norm_fun(const var& q) { 
    return sqrt(q.rhou*q.rhou + q.rhov*q.rhov)/q.rho; 
}
float get_mach_fun(const var& q) { 
    float p = (1.4 - 1.)*(q.rhoe - 0.5/q.rho*(q.rhou*q.rhou + q.rhov*q.rhov));
    float c = sqrt(1.4 * p/q.rho);
    return get_u_norm_fun(q)/c; 
}

std::vector<float (*)(const var&)> plot_field_map = {
    &get_rho_fun, &get_u_fun, &get_v_fun, &get_p_fun, 
    &get_u_norm_fun, &get_mach_fun
};
std::map<std::string, uint> plot_field_ref = {
    {"rho", 0}, {"u", 1}, {"v", 2}, {"p", 3},
    {"U", 4}, {"mach", 5}
};




void ExampleWindow::remake_buffers() {

    glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
    glBufferData(GL_ARRAY_BUFFER, VERTEX_NUMBER*6*sizeof(GLfloat), &vertices[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ARRAY_BUFFER, m_VBO_q);
    glBufferData(GL_ARRAY_BUFFER, VERTEX_NUMBER_Q*6*sizeof(GLfloat), &vertices_q[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ARRAY_BUFFER, m_VBO_l);
    glBufferData(GL_ARRAY_BUFFER, VERTEX_NUMBER_L*6*sizeof(GLfloat), &vertices_l[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void ExampleWindow::clear_vertices_array() {
    vertices.clear();
    VERTEX_NUMBER = 0;

    vertices_q.clear();
    VERTEX_NUMBER_Q = 0;

    vertices_l.clear();
    VERTEX_NUMBER_L = 0;

    remake_buffers();
}



#include <gui/plot.h>




void ExampleWindow::demo_vertices_array() {

    // Draw filled triangles
    DRAW_POLY_AS_LINES = false;

    /*
    VERTEX_NUMBER = 6;
    std::vector<GLfloat> init = {
        // POSITION		

        // Triangle 1
        // X, Y, Z position of vertex #1
        -0.5f, -0.0f, 0.0f,	1.0f, 0.0f, 0.0f,	
        
        // X, Y, Z position of vertex #2
        0.5, -0.0f, 0.0f, 0.0f, 1.0f, 0.0f,	
        
        // X, Y, Z position of vertex #3
        0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f,	

        // Triangle 2
        // X, Y, Z position of vertex #1
        0.5f, 0.0f, 0.0f,	1.0f, 0.0f, 0.0f,	
        
        // X, Y, Z position of vertex #2
        -0.5, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,	
        
        // X, Y, Z position of vertex #3
        0.0f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f,	
    };

    uint n_tags = VERTEX_NUMBER * 6;
    */
    VERTEX_NUMBER_L = 6;
    std::vector<GLfloat> init = {
        // POSITION		

        // Line 1
        // X, Y, Z position of vertex #1
        -0.5f, -0.0f, 0.0f,	1.0f, 0.0f, 0.0f,	
        
        // X, Y, Z position of vertex #2
        0.5, -0.0f, 0.0f, 0.0f, 1.0f, 0.0f,	

        // Line 2

        // X, Y, Z position of vertex #2
        0.5, -0.0f, 0.0f, 0.0f, 1.0f, 0.0f,	
        
        // X, Y, Z position of vertex #3
        0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f,	

        // Line 3

        // X, Y, Z position of vertex #2
        0.5, -0.0f, 0.0f, 0.0f, 1.0f, 0.0f,	
        
        // X, Y, Z position of vertex #3
        0.0f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f,		
    };

    uint n_tags = VERTEX_NUMBER_L * 6;

    std::cout << "demo vertices array" << std::endl;
    vertices_l.resize(n_tags);
    for (int i=0; i<init.size(); ++i) {
        vertices_l[i] = init[i];
    }
    remake_buffers();
}




void ExampleWindow::gen_mesh_gl_arrays() {
    // Populate vertices array with mesh info

    // Draw polygons as lines
    DRAW_POLY_AS_LINES = true;

    // Define vertices array for tris

    // Get the number of cells in mesh
    uint n_cells = tri_ref.size();
    uint n_vertex = n_cells * 3;
    uint n_tags = n_vertex * 6; // Color tags are added

    VERTEX_NUMBER = n_vertex;
    //vertices = (GLfloat*)realloc(&vertices[0], n_tags*sizeof(GLfloat));
    if (n_tags != vertices.size()) {
        vertices.resize(n_tags);
    }

    uint n_tri_added = 0;
    for (uint i : tri_ref) {
        std::array<vec2, 3> ns = {
            m_mesh.nodes[m_mesh.facesNodes0[i]],
            m_mesh.nodes[m_mesh.facesNodes1[i]],
            m_mesh.nodes[m_mesh.facesNodes2[i]]
        };
        for (uint j=0; j<ns.size(); ++j) {
            uint k = n_tri_added*6*3 + j*6;   // This vertex start position

            // Position
            vertices[k] = ns[j].x;
            vertices[k+1] = ns[j].y;
            vertices[k+2] = 0.;

            // Color
            vertices[k+3] = 1.0f;
            vertices[k+4] = 1.0f;
            vertices[k+5] = 1.0f;
        }
        n_tri_added += 1;
    }

    // Define vertices array for quads
    // Draw quads as lines
    n_cells = quad_ref.size();
    n_vertex = n_cells * 8;
    n_tags = n_vertex * 6; // Color tags are added

    VERTEX_NUMBER_L = n_vertex;
    if (vertices_l.size() != n_tags) {
        vertices_l.resize(n_tags);
    }

    uint n_quad_added = 0;
    for (uint i : quad_ref) {
        std::array<vec2, 8> ns = {
            m_mesh.nodes[m_mesh.facesNodes0[i]],
            m_mesh.nodes[m_mesh.facesNodes1[i]],

            m_mesh.nodes[m_mesh.facesNodes1[i]],
            m_mesh.nodes[m_mesh.facesNodes2[i]],

            m_mesh.nodes[m_mesh.facesNodes2[i]],
            m_mesh.nodes[m_mesh.facesNodes3[i]],

            m_mesh.nodes[m_mesh.facesNodes3[i]],
            m_mesh.nodes[m_mesh.facesNodes0[i]]
        };
        for (uint j=0; j<ns.size(); ++j) {
            uint k = n_quad_added*8*6 + j*6;   // This vertex start position

            // Position
            vertices_l[k] = ns[j].x;
            vertices_l[k+1] = ns[j].y;
            vertices_l[k+2] = 0.;

            // Color
            vertices_l[k+3] = 1.0f;
            vertices_l[k+4] = 1.0f;
            vertices_l[k+5] = 1.0f;
        }
        n_quad_added += 1;
    }
    

    std::cout << "generated mesh array" << std::endl;
    remake_buffers();
}


void ExampleWindow::gen_solution_gl_arrays() {
    // Populate vertices array with mesh info


    if (solution.size() == 0) {
        std::cout << "Empty solution, solve before displaying" << std::endl;
    } else {

        // Draw polygons as lines
        DRAW_POLY_AS_LINES = false;

        // Get the plot variable function pointer
        auto plot_field_fun = &get_rho_fun;

        std::string plot_field_name = (*(m_refTreeModel->children()[3]->children()[1])).get_value(m_Columns.m_col_text);
        plot_field_fun = plot_field_map[plot_field_ref[plot_field_name]];

        // Get the colormap
        auto cmap = &cmap_rainbow;

        std::string cmap_name = (*(m_refTreeModel->children()[3]->children()[0])).get_value(m_Columns.m_col_text);
        cmap = cmap_map[cmap_ref[cmap_name]];

        // Get the maximal and minimal value in solution
        float plot_field_min = plot_field_fun(solution[0]);
        float plot_field_max = plot_field_fun(solution[0]);
        for (int i=1; i<solution.size(); ++i) {
            plot_field_min = std::min(plot_field_min, plot_field_fun(solution[i]));
            plot_field_max = std::max(plot_field_max, plot_field_fun(solution[i]));
        }

        // Define vertices array for tris

        // Get the number of cells in mesh
        uint n_cells = tri_ref.size();
        uint n_vertex = n_cells * 3;
        uint n_tags = n_vertex * 6; // Color tags are added

        VERTEX_NUMBER = n_vertex;
        //vertices = (GLfloat*) realloc(&vertices[0], n_tags*sizeof(GLfloat));
        if (n_tags != vertices.size()) {
            vertices.resize(n_tags);
            std::cout << "Resized vertices" << std::endl;
        }

        uint n_tri_added = 0;
        for (uint i : tri_ref) {
            std::array<vec2, 3> ns = {
                m_mesh.nodes[m_mesh.facesNodes0[i]],
                m_mesh.nodes[m_mesh.facesNodes1[i]],
                m_mesh.nodes[m_mesh.facesNodes2[i]]
            };

            float qi = plot_field_fun(solution[i]);

            auto col = cmap(
                (qi - plot_field_min)/(plot_field_max - plot_field_min)
            );
            
            for (uint j=0; j<ns.size(); ++j) {
                uint k =n_tri_added*6*3 + j*6;   // This vertex start position
                
                // Position
                vertices[k] = ns[j].x;
                vertices[k+1] = ns[j].y;
                vertices[k+2] = 0.;

                // Color
                vertices[k+3] = col[0];
                vertices[k+4] = col[1];
                vertices[k+5] = col[2];
            }
            n_tri_added += 1;
        }

        // Define vertices array for quads

        // Get the number of cells in mesh
        n_cells = quad_ref.size();
        n_vertex = n_cells * 3 * 2;
        n_tags = n_vertex * 6; // Color tags are added

        VERTEX_NUMBER_Q = n_vertex;
        if (vertices_q.size() != n_tags) {
            vertices_q.resize(n_tags);
            std::cout << "Resized vertices_q" << std::endl;
        }

        uint n_quad_added = 0;

        for (uint i : quad_ref) {

            float qi = plot_field_fun(solution[i]);

            auto col = cmap(
                (qi - plot_field_min)/(plot_field_max - plot_field_min)
            );

            std::array<vec2, 3> ns = {
                m_mesh.nodes[m_mesh.facesNodes0[i]],
                m_mesh.nodes[m_mesh.facesNodes1[i]],
                m_mesh.nodes[m_mesh.facesNodes2[i]]
            };
            
            for (uint j=0; j<ns.size(); ++j) {
                uint k = n_quad_added*6*3 + j*6;   // This vertex start position
                
                // Position
                vertices_q[k] = ns[j].x;
                vertices_q[k+1] = ns[j].y;
                vertices_q[k+2] = 0.;

                // Color
                vertices_q[k+3] = col[0];
                vertices_q[k+4] = col[1];
                vertices_q[k+5] = col[2];
            }
            n_quad_added += 1;

            std::array<vec2, 3> ns2 = {
                m_mesh.nodes[m_mesh.facesNodes2[i]],
                m_mesh.nodes[m_mesh.facesNodes3[i]],
                m_mesh.nodes[m_mesh.facesNodes0[i]]
            };
            
            for (uint j=0; j<ns2.size(); ++j) {
                uint k = n_quad_added*6*3 + j*6;   // This vertex start position
                
                // Position
                vertices_q[k] = ns2[j].x;
                vertices_q[k+1] = ns2[j].y;
                vertices_q[k+2] = 0.;

                // Color
                vertices_q[k+3] = col[0];
                vertices_q[k+4] = col[1];
                vertices_q[k+5] = col[2];
            }
            n_quad_added += 1;
        }

        std::cout << "generated solution array" << std::endl;
    }

    remake_buffers();
}




void ExampleWindow::realize() {
    m_GLArea.make_current();
    try {
        m_GLArea.throw_if_error();

        if (GLEW_OK != glewInit())
        {
            // If the initalization is not successful, print out the message and exit the program with return value EXIT_FAILURE
            std::cout << "Failed to initialize GLEW" << std::endl;
        }

        // Triangles
        glGenVertexArrays(1, &m_VAO);
        glGenBuffers(1, &m_VBO);
        glBindVertexArray(m_VAO);

        glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
        glBufferData(GL_ARRAY_BUFFER, VERTEX_NUMBER*6*sizeof(GLfloat), &vertices[0], GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(0);
        // color attribute
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        // Quads
        glGenVertexArrays(1, &m_VAO_q);
        glGenBuffers(1, &m_VBO_q);
        glBindVertexArray(m_VAO_q);

        glBindBuffer(GL_ARRAY_BUFFER, m_VBO_q);
        glBufferData(GL_ARRAY_BUFFER, VERTEX_NUMBER_Q*6*sizeof(GLfloat), &vertices_q[0], GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(0);
        // color attribute
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        // Lines
        glGenVertexArrays(1, &m_VAO_l);
        glGenBuffers(1, &m_VBO_l);
        glBindVertexArray(m_VAO_l);

        glBindBuffer(GL_ARRAY_BUFFER, m_VBO_l);
        glBufferData(GL_ARRAY_BUFFER, VERTEX_NUMBER_L*6*sizeof(GLfloat), &vertices_l[0], GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(0);
        // color attribute
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        glUseProgram(m_Program);

        init_shaders();

    } catch(const Gdk::GLError& gle) {
        std::cout << "An error occured making the context current during realize:" << std::endl;
        std::cout << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
    }
}


void ExampleWindow::unrealize() {
    m_GLArea.make_current();
    try {
        m_GLArea.throw_if_error();

        // Triangles
        glDeleteVertexArrays(1, &m_VAO);
        glDeleteBuffers(1, &m_VBO);

        // Quads
        glDeleteVertexArrays(1, &m_VAO_q);
        glDeleteBuffers(1, &m_VBO_q);

        // Lines
        glDeleteVertexArrays(1, &m_VAO_l);
        glDeleteBuffers(1, &m_VBO_l);

        // Delete buffers and program
        glDeleteProgram(m_Program);
    } catch(const Gdk::GLError& gle) {
        std::cout << "An error occured making the context current during unrealize" << std::endl;
        std::cout << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
    }
}

bool ExampleWindow::render(const Glib::RefPtr<Gdk::GLContext>& context) {
    try {
        m_GLArea.throw_if_error();

        // Specifies the RGBA values which will be used by glClear to clear the color buffer
        glUseProgram(m_Program);

		glClearColor(clear_color[0], clear_color[1], clear_color[2], clear_color[3]);

		glClear(GL_COLOR_BUFFER_BIT);

        auto transform_location = glGetUniformLocation(m_Program, "transform");

        // Build transform matrix
        transform[0] = 400.0f * camera->scale / camera->width;
        transform[5] = 400.0f * camera->scale / camera->height;
        transform[12] = camera->x_cam;
        transform[13] = camera->y_cam;

        glUniformMatrix4fv(transform_location, 1, GL_FALSE, transform);

        // Triangles
        glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
		glBindVertexArray(m_VAO);

        if (DRAW_POLY_AS_LINES) 
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawArrays(GL_TRIANGLES, 0, VERTEX_NUMBER);
        if (DRAW_POLY_AS_LINES)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        
        
        // Quads
        glBindBuffer(GL_ARRAY_BUFFER, m_VBO_q);
		glBindVertexArray(m_VAO_q);

        if (DRAW_POLY_AS_LINES) 
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawArrays(GL_TRIANGLES, 0, VERTEX_NUMBER_Q);
        if (DRAW_POLY_AS_LINES)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        
        // Lines
        glBindBuffer(GL_ARRAY_BUFFER, m_VBO_l);
		glBindVertexArray(m_VAO_l);

		glDrawArrays(GL_LINES, 0, VERTEX_NUMBER_L);
		
        
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        glUseProgram(0);

        glFlush();

        return true;
    } catch(const Gdk::GLError& gle) {
        std::cout << "An error occurred in the render callback of the GLArea" << std::endl;
        std::cout << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
        return false;
    }
}



void ExampleWindow::init_shaders() {

    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    // check for shader compile errors
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // fragment shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    // check for shader compile errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // link shaders
    m_Program = glCreateProgram();
    glAttachShader(m_Program, vertexShader);
    glAttachShader(m_Program, fragmentShader);
    glLinkProgram(m_Program);
    // check for linking errors
    glGetProgramiv(m_Program, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(m_Program, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

}



#endif