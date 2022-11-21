#ifndef EULER_GUI_ACTIONS_H
#define EULER_GUI_ACTIONS_H


#include <euler/euler.h>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include <gtkmm.h>

#include <gui/core.h>


ExampleWindow::~ExampleWindow()
{
  clear_vertices_array();
}


void ExampleWindow::on_action_file_quit()
{
    clear_vertices_array();
    hide(); //Closes the main window to stop the app->run().
}

void ExampleWindow::on_action_others()
{
  std::cout << "A menu item was selected." << std::endl;
}


void ExampleWindow::on_action_file_new()
{
  std::cout << "No new file method implemented yet." << std::endl;
}

void ExampleWindow::on_action_file_open()
{
  std::cout << "No open file method implemented yet." << std::endl;
}

void ExampleWindow::on_action_file_save()
{
  std::cout << "No save file method implemented yet." << std::endl;
}


void ExampleWindow::on_action_mesh_generate() {
  std::cout << "Generate mesh command called." << std::endl;

  std::string geo_file_name = get_in_tree("Geometry", "Geo file");
  std::string mesh_file_name = get_in_tree("Geometry", "Mesh file");;
  
  // Run system command to generate mesh file
  std::cout << "Generate mesh (gmsh command output)\n===================================" << std::endl;
  std::string command = "gmsh " + geo_file_name + " -2 -o " + mesh_file_name;
  system(command.c_str());
  std::cout << "===================================\n" << std::endl;

}


void ExampleWindow::on_action_mesh_load() {
  std::cout << "Load mesh command called." << std::endl;

  std::string mesh_file_name = get_in_tree("Geometry", "Mesh file");

  // Load mesh
  std::cout << "===================================\n" << std::endl;
  m_mesh = mesh(mesh_file_name);
  std::cout << "===================================\n" << std::endl;

  // Fill tri ref and quad ref vectors
  tri_ref.clear();
  quad_ref.clear();
  for (uint k = 0; k<m_mesh.facesAreas.size(); ++k) {
    if (m_mesh.faceIsTriangle[k]) {
      tri_ref.push_back(k);
    } else {
      quad_ref.push_back(k);
    }
  }
  std::cout << "Tri ref size  = " << tri_ref.size() << std::endl;
  std::cout << "Quad ref size = " << quad_ref.size() << std::endl;

  // Clear boundary conditions
  std::vector<uint> row_vector;

  auto children = m_refTreeModel->children();
  for (auto iter = children.begin(), end = children.end(); iter != end; ++iter)
  {
    auto row = *iter;
    auto row_id = row[m_Columns.m_col_number];
    if (row_id == 2) {
      auto children2 = row->children();

      for (auto iter = children2.begin(), end = children2.end(); iter != end; ++iter)
      {
        auto crow = *iter;
        row_vector.push_back(crow[m_Columns.m_col_number]);
      }
    }
  }

  // Fill new boundary conditions
  for (uint row_id : row_vector) {
    for (auto iter = children.begin(), end = children.end(); iter != end; ++iter)
    {
      auto row = *iter;
      auto children2 = row->children();
      for (auto iter = children2.begin(), end = children2.end(); iter != end; ++iter)
      {
        auto crow = *iter;
        if (crow[m_Columns.m_col_number] == row_id) {
          m_refTreeModel->erase(iter);
          break;
        }
      }
    }
  }


  // Set boundary conditions
  children = m_refTreeModel->children();
  for (auto iter = children.begin(), end = children.end(); iter != end; ++iter)
  {
    auto row = *iter;
    auto row_id = row[m_Columns.m_col_number];
    if (row_id == 2) {
      uint k = 1;
      Gtk::TreeModel::Row childrow;
      for (auto s : m_mesh.uniqueBoundsTags) {
        childrow = *(m_refTreeModel->append(row.children()));
        childrow[m_Columns.m_col_number] = 20+k;
        childrow[m_Columns.m_col_name] = s;
        childrow[m_Columns.m_col_text] = "wall";
        k += 1;
      }
    }
  }

  // Display mesh
  gen_mesh_gl_arrays();

}



void ExampleWindow::on_action_solve() {
  std::cout << "Solve command called." << std::endl;

  uint n_cells = m_mesh.facesAreas.size();
  inputVar q0i;

  // Set q0
  auto q0_text = get_in_tree("Solver Settings", "Init. Conditions");
  auto line_cut = cut_str(q0_text, ' ');
  q0i = inputVar(
    std::stod(line_cut[0]), std::stod(line_cut[1]), std::stod(line_cut[2]), std::stod(line_cut[3])
  );

  // Generate initial conditions
  std::vector<inputVar> q0(n_cells);
  for (uint i=0; i<n_cells; ++i) {
    q0[i] = q0i;
  }

  constants c = constants(1.4);

  solver = steadyRk5Solver(m_mesh, q0, c);
  
  // Set CFL
  solver.set_max_steps(std::stod(get_in_tree("Solver Settings", "CFL")));

  // Set max steps
  solver.set_max_steps(std::stoi(get_in_tree("Solver Settings", "Max Steps")));
  
  // Set order
  solver.set_order(std::stoi(get_in_tree("Solver Settings", "Order")));
  
  // Set boundary conditions
  for (auto bc_name : m_mesh.uniqueBoundsTags) {
    std::string bc_text = get_in_tree("Boundary Conditions", bc_name);
    std::cout << bc_text << std::endl;
    auto line_cut = cut_str(bc_text, ' ');
    if (line_cut[0] != "wall") {
      solver.set_bc(
        bc_name, 
        line_cut[0], 
        inputVar(std::stod(line_cut[1]), std::stod(line_cut[2]), std::stod(line_cut[3]), std::stod(line_cut[4]))
      );
    } else {
      solver.set_bc(
        bc_name, 
        line_cut[0]
      );
    }
  }

  // Init simulation

  solver.init_solver(SOLVER_VERB);
  solution = solver.get_raw_solution();

  SOLVER_RUNNING = true;

  Glib::signal_idle().connect( sigc::mem_fun(*this, &ExampleWindow::run_solver_step), 0 );

}


bool ExampleWindow::run_solver_step() {
  if (SOLVER_RUNNING) {
    if (solver.convergence_check()) {
      solver.internal_solver_step(SOLVER_VERB);
    } else {
      solver.end_solver(SOLVER_VERB);
      SOLVER_RUNNING = false;
    }
    std::cout << activeGLArea << std::endl;

    if (activeGLArea == "results") {
      solution = solver.get_raw_solution();

      gen_solution_gl_arrays();

      m_GLArea.queue_draw();
    } else if (activeGLArea == "plot-residuals") {
      gen_plot_gl_arrays();

      m_GLArea.queue_draw();
    }

    return true;
  } else {
    return false;
  }
}


void ExampleWindow::on_action_mesh_show() {
  clear_vertices_array();
  gen_mesh_gl_arrays();
  activeGLArea = "mesh";
  camera = &solver_camera;
  m_GLArea.queue_render();
}

void ExampleWindow::on_action_results_show() {
  clear_vertices_array();
  solution = solver.get_raw_solution();
  gen_solution_gl_arrays();
  activeGLArea = "results";
  camera = &solver_camera;
  m_GLArea.queue_render();
}

void ExampleWindow::on_action_plot_show() {
  clear_vertices_array();
  gen_plot_gl_arrays();
  activeGLArea = "plot-residuals";
  camera = &plot_camera;
  m_GLArea.queue_render();
}


void ExampleWindow::on_action_debug_show() {
  clear_vertices_array();
  demo_vertices_array();
  activeGLArea = "debug";
  camera = &plot_camera;
  m_GLArea.queue_render();
}



void ExampleWindow::glarea_resize(GLuint width, GLuint height) {

    solver_camera.width = width;
    solver_camera.height = height;

    plot_camera.width = width;
    plot_camera.height = height;

    m_GLArea.queue_render();
    std::cout << "GLArea window resized " << width << " " << height << std::endl;
}



bool ExampleWindow::glarea_pan(GdkEventMotion* event) {
  GLint dx = camera->x_offset - (event->x - camera->x_start);
  GLint dy = camera->y_offset - (event->y - camera->y_start);

  camera->x_cam = camera->x_cam - 2.*((GLfloat) dx)/((GLfloat) camera->width);
  camera->y_cam = camera->y_cam + 2.*((GLfloat) dy)/((GLfloat) camera->height);

  camera->x_offset = event->x - camera->x_start;
  camera->y_offset = event->y - camera->y_start;

  std::cout << "Testing window pan x=" << 
    camera->x_offset << " y=" << camera->y_offset << std::endl;

  m_GLArea.queue_render();
  return true;
}


bool ExampleWindow::glarea_start_pan(GdkEventButton* event) {

  std::cout << "Testing window start pan x=" << event->x << " y=" << event->y << std::endl;
  
  camera->x_start = event->x;
  camera->y_start = event->y;

  camera->x_offset = 0;
  camera->y_offset = 0;

  return true;
}


bool ExampleWindow::glarea_scroll(GdkEventScroll* event) {

  std::cout << "Testing window scroll ";

  if (event->direction == GDK_SCROLL_UP) {
    std::cout << "Scroll up" << std::endl;
    camera->scale *= 1.1;
  } else if (event->direction == GDK_SCROLL_DOWN) {
    std::cout << "Scroll down" << std::endl;
    camera->scale *= 1./1.1;
  }

  m_GLArea.queue_render();

  return true;
}



#endif