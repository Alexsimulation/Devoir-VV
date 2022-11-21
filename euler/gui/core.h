#ifndef EULER_GUI_CORE_H
#define EULER_GUI_CORE_H


#include <euler/euler.h>
#include <gtkmm.h>

#define GLEW_STATIC
#include <GL/glew.h>

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <stdexcept>



struct camera_t {

  GLint x_offset = 0;
  GLint y_offset = 0;

  GLfloat x_cam = 0.;
  GLfloat y_cam = 0.;
  GLfloat scale = 1.;

  GLint x_start;
  GLint y_start;

  GLuint width;
  GLuint height;

};



class ExampleWindow : public Gtk::Window
{
public:
  ExampleWindow(const Glib::RefPtr<Gtk::Application>& app);
  virtual ~ExampleWindow();

private:
  //Signal handlers:
  void on_action_file_new();
  void on_action_file_open();
  void on_action_file_save();
  void on_action_file_quit();
  void on_action_others();

  void on_action_mesh_generate();
  void on_action_mesh_load();

  void on_action_solve();

  void on_action_mesh_show();
  void on_action_results_show();
  void on_action_plot_show();
  void on_action_debug_show();


  // Euler solver

  mesh m_mesh;
  std::vector<var> solution;
  steadyRk5Solver solver;

  bool SOLVER_RUNNING = false;
  bool SOLVER_VERB = true;

  bool run_solver_step();


  // Gui app 
  std::string activeFile = "new.txt";

  //Child widgets:
  Gtk::Box m_Box;
  Gtk::Paned m_Paned;
  Gtk::Grid m_Grid;

  Gtk::TextView m_TextView;

  Gtk::TreeView m_TreeView;

  Gtk::GLArea m_GLArea;


  GLfloat GLArea_scale;
  GLfloat GLArea_offset_x;
  GLfloat GLArea_offset_y;

  GLuint m_VAO {0};
  GLuint m_VBO {0};
  std::vector<uint> tri_ref;

  GLuint m_VAO_q {0};
  GLuint m_VBO_q {0};
  std::vector<uint> quad_ref;

  GLuint m_VAO_l {0};
  GLuint m_VBO_l {0};

  GLuint m_Program {0};
  GLuint m_Mvp {0};

  std::array<GLfloat, 4> clear_color = {0.1f, 0.1f, 0.1f, 1.0f};

  std::string activeGLArea = "none";

  std::string get_in_tree(std::string section, std::string sub_section);

  void realize();
  void unrealize();
  bool render(const Glib::RefPtr<Gdk::GLContext>& context);
  void init_shaders();
  void remake_buffers();
  void clear_vertices_array();
  void demo_vertices_array();

  void ScreenToWorldSpace(int mousex, int mousey);

  void gen_mesh_gl_arrays();
  void gen_solution_gl_arrays();
  void gen_plot_gl_arrays();

  void glarea_resize(GLuint width, GLuint height);

  bool glarea_start_pan(GdkEventButton* event);
  bool glarea_pan(GdkEventMotion* event);
  bool glarea_scroll(GdkEventScroll* event);

  camera_t solver_camera;
  camera_t plot_camera;
  camera_t* camera = &solver_camera;


  Glib::RefPtr<Gtk::Builder> m_refBuilder;
  Glib::RefPtr<Gio::SimpleActionGroup> m_refActionGroup;
  Glib::RefPtr<Gio::SimpleAction> m_refActionRain;
  Glib::RefPtr<Gtk::TreeStore> m_refTreeModel;

  class ModelColumns : public Gtk::TreeModelColumnRecord
  {
  public:

    ModelColumns()
      { add(m_col_number); add(m_col_name); add(m_col_text); }

    Gtk::TreeModelColumn<int> m_col_number;
    Gtk::TreeModelColumn<Glib::ustring> m_col_name;
    Gtk::TreeModelColumn<Glib::ustring> m_col_text;
  };

  ModelColumns m_Columns;

};





std::string readFileIntoString(const std::string& path) {
    std::ifstream input_file(path);
    if (!input_file.is_open()) {
        std::cerr << "Could not open the file - '"
             << path << "'" << std::endl;
        return "";
    } else {
      return std::string((std::istreambuf_iterator<char>(input_file)), std::istreambuf_iterator<char>());
    }
}





ExampleWindow::ExampleWindow(const Glib::RefPtr<Gtk::Application>& app)
: m_Box(Gtk::ORIENTATION_VERTICAL)
{
  set_title("Euler Solver");
  const int SIZE_W = 720;
  const int SIZE_H = 480;
  set_default_size(SIZE_W, SIZE_H);

  add(m_Box); //We can put a MenuBar at the top of the box and other stuff below it.

  //Define the actions:
  m_refActionGroup = Gio::SimpleActionGroup::create();

  m_refActionGroup->add_action("new",
    sigc::mem_fun(*this, &ExampleWindow::on_action_file_new) );
  m_refActionGroup->add_action("open",
    sigc::mem_fun(*this, &ExampleWindow::on_action_file_open) );
  m_refActionGroup->add_action("save",
    sigc::mem_fun(*this, &ExampleWindow::on_action_file_save) );

  m_refActionGroup->add_action("quit",
    sigc::mem_fun(*this, &ExampleWindow::on_action_file_quit) );

  m_refActionGroup->add_action("cut",
    sigc::mem_fun(*this, &ExampleWindow::on_action_others) );
  m_refActionGroup->add_action("copy",
    sigc::mem_fun(*this, &ExampleWindow::on_action_others) );
  m_refActionGroup->add_action("paste",
    sigc::mem_fun(*this, &ExampleWindow::on_action_others) );
  
  m_refActionGroup->add_action("generate_mesh",
    sigc::mem_fun(*this, &ExampleWindow::on_action_mesh_generate) );
  m_refActionGroup->add_action("load_mesh",
    sigc::mem_fun(*this, &ExampleWindow::on_action_mesh_load) );
  
  m_refActionGroup->add_action("solve",
    sigc::mem_fun(*this, &ExampleWindow::on_action_solve) );
  
  m_refActionGroup->add_action("show_mesh",
    sigc::mem_fun(*this, &ExampleWindow::on_action_mesh_show) );
  m_refActionGroup->add_action("show_results",
    sigc::mem_fun(*this, &ExampleWindow::on_action_results_show) );
  m_refActionGroup->add_action("show_plot",
    sigc::mem_fun(*this, &ExampleWindow::on_action_plot_show) );
  m_refActionGroup->add_action("show_debug",
    sigc::mem_fun(*this, &ExampleWindow::on_action_debug_show) );

  insert_action_group("example", m_refActionGroup);

  //Define how the actions are presented in the menus and toolbars:
  m_refBuilder = Gtk::Builder::create();

  //Layout the actions in a menubar and toolbar:
  const char* ui_info =
    "<interface>"
    "  <menu id='menubar'>"
    "    <submenu>"
    "      <attribute name='label' translatable='yes'>_File</attribute>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_New</attribute>"
    "        <attribute name='action'>example.new</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;n</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Open</attribute>"
    "        <attribute name='action'>example.open</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;o</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Save</attribute>"
    "        <attribute name='action'>example.save</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;s</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Quit</attribute>"
    "        <attribute name='action'>example.quit</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;q</attribute>"
    "      </item>"
    "    </submenu>"
    "    <submenu>"
    "      <attribute name='label' translatable='yes'>_Edit</attribute>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Cut</attribute>"
    "        <attribute name='action'>example.cut</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;x</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Copy</attribute>"
    "        <attribute name='action'>example.copy</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;c</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Paste</attribute>"
    "        <attribute name='action'>example.paste</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;v</attribute>"
    "      </item>"
    "    </submenu>"
    "    <submenu>"
    "      <attribute name='label' translatable='yes'>_Mesh</attribute>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Generate Mesh</attribute>"
    "        <attribute name='action'>example.generate_mesh</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;g</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Load Mesh</attribute>"
    "        <attribute name='action'>example.load_mesh</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;l</attribute>"
    "      </item>"
    "    </submenu>"
    "    <submenu>"
    "      <attribute name='label' translatable='yes'>_Solver</attribute>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Solve</attribute>"
    "        <attribute name='action'>example.solve</attribute>"
    "      </item>"
    "    </submenu>"
    "    <submenu>"
    "      <attribute name='label' translatable='yes'>_Display</attribute>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Mesh</attribute>"
    "        <attribute name='action'>example.show_mesh</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;m</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Results</attribute>"
    "        <attribute name='action'>example.show_results</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;r</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Plot</attribute>"
    "        <attribute name='action'>example.show_plot</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;p</attribute>"
    "      </item>"
    "      <item>"
    "        <attribute name='label' translatable='yes'>_Debug</attribute>"
    "        <attribute name='action'>example.show_debug</attribute>"
    "        <attribute name='accel'>&lt;Primary&gt;d</attribute>"
    "      </item>"
    "    </submenu>"
    "  </menu>"
    "</interface>";

  // When the menubar is a child of a Gtk::Window, keyboard accelerators are not
  // automatically fetched from the Gio::Menu.
  // See the examples/book/menus/main_menu example for an alternative way of
  // adding the menubar when using Gtk::ApplicationWindow.
  // Gtk::Application::set_accel_for_action() is new in gtkmm 3.11.9.
  app->set_accel_for_action("example.new", "<Primary>n");
  app->set_accel_for_action("example.open", "<Primary>o");
  app->set_accel_for_action("example.save", "<Primary>s");
  app->set_accel_for_action("example.quit", "<Primary>q");
  app->set_accel_for_action("example.cut", "<Primary>x");
  app->set_accel_for_action("example.copy", "<Primary>c");
  app->set_accel_for_action("example.paste", "<Primary>v");
  app->set_accel_for_action("example.generate_mesh", "<Primary>g");
  app->set_accel_for_action("example.load_mesh", "<Primary>l");
  app->set_accel_for_action("example.show_mesh", "<Primary>m");
  app->set_accel_for_action("example.show_results", "<Primary>r");
  app->set_accel_for_action("example.show_plot", "<Primary>p");
  app->set_accel_for_action("example.show_debug", "<Primary>d");

  try
  {
    m_refBuilder->add_from_string(ui_info);
  }
  catch(const Glib::Error& ex)
  {
    std::cerr << "Building menus and toolbar failed: " <<  ex.what();
  }

  //Get the menubar:
  auto object = m_refBuilder->get_object("menubar");
  auto gmenu = Glib::RefPtr<Gio::Menu>::cast_dynamic(object);
  if (!gmenu)
    g_warning("GMenu not found");
  else
  {
    auto pMenuBar = Gtk::make_managed<Gtk::MenuBar>(gmenu);

    //Add the MenuBar to the window:
    m_Box.pack_start(*pMenuBar, Gtk::PACK_SHRINK);
  }

  // Add textview
  //m_TextView.set_hexpand (true);
  //m_TextView.set_size_request(100, 100);
  m_TreeView.set_size_request(300, 200);

  m_refTreeModel = Gtk::TreeStore::create(m_Columns);
  m_TreeView.set_model(m_refTreeModel);
  Gtk::TreeModel::Row row;
  Gtk::TreeModel::Row childrow;

  // Add geometry row
  row = *(m_refTreeModel->append());
  row[m_Columns.m_col_number] = 1;
  row[m_Columns.m_col_name] = "Geometry";
  row[m_Columns.m_col_text] = "";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 11;
  childrow[m_Columns.m_col_name] = "Geo file";
  childrow[m_Columns.m_col_text] = "new.geo";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 12;
  childrow[m_Columns.m_col_name] = "Mesh file";
  childrow[m_Columns.m_col_text] = "new.su2";

  // Add boundary conditions row
  row = *(m_refTreeModel->append());
  row[m_Columns.m_col_number] = 2;
  row[m_Columns.m_col_name] = "Boundary Conditions";
  row[m_Columns.m_col_text] = "";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 21;
  childrow[m_Columns.m_col_name] = "airfoil";
  childrow[m_Columns.m_col_text] = "wall";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 22;
  childrow[m_Columns.m_col_name] = "test";
  childrow[m_Columns.m_col_text] = "wall";

  // Add solver settings row
  row = *(m_refTreeModel->append());
  row[m_Columns.m_col_number] = 3;
  row[m_Columns.m_col_name] = "Solver Settings";
  row[m_Columns.m_col_text] = "";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 31;
  childrow[m_Columns.m_col_name] = "Init. Conditions";
  childrow[m_Columns.m_col_text] = "1 0 0 1";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 32;
  childrow[m_Columns.m_col_name] = "Time scheme";
  childrow[m_Columns.m_col_text] = "steady";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 33;
  childrow[m_Columns.m_col_name] = "CFL";
  childrow[m_Columns.m_col_text] = "5";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 34;
  childrow[m_Columns.m_col_name] = "Order";
  childrow[m_Columns.m_col_text] = "1";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 35;
  childrow[m_Columns.m_col_name] = "Max Steps";
  childrow[m_Columns.m_col_text] = "1000";

  // Add results settings row
  row = *(m_refTreeModel->append());
  row[m_Columns.m_col_number] = 4;
  row[m_Columns.m_col_name] = "Results";
  row[m_Columns.m_col_text] = "";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 41;
  childrow[m_Columns.m_col_name] = "Colormap";
  childrow[m_Columns.m_col_text] = "rainbow";

  childrow = *(m_refTreeModel->append(row.children()));
  childrow[m_Columns.m_col_number] = 42;
  childrow[m_Columns.m_col_name] = "Plot field";
  childrow[m_Columns.m_col_text] = "rho";

  // Add columns
  m_TreeView.append_column("Name", m_Columns.m_col_name);
  m_TreeView.append_column_editable("Value", m_Columns.m_col_text);

  // Set columns to be resized
  for(guint i = 0; i < 2; i++)
  {
    auto column = m_TreeView.get_column(i);
    column->set_resizable();
  }

  m_Paned.add1(m_TreeView);

  // Add GLArea
  m_GLArea.set_hexpand (true);
  m_GLArea.set_vexpand (true);
  m_GLArea.set_size_request(200, 200);
  m_GLArea.set_auto_render(true);
  m_Paned.add2(m_GLArea);

  // Connect gl area signals
  m_GLArea.signal_realize().connect(sigc::mem_fun(*this, &ExampleWindow::realize));
  // Important that the unrealize signal calls our handler to clean up
  // GL resources _before_ the default unrealize handler is called (the "false")
  m_GLArea.signal_unrealize().connect(sigc::mem_fun(*this, &ExampleWindow::unrealize), false);
  m_GLArea.signal_render().connect(sigc::mem_fun(*this, &ExampleWindow::render), false);
  m_GLArea.signal_resize().connect(sigc::mem_fun(*this, &ExampleWindow::glarea_resize), false);

  m_GLArea.add_events(Gdk::BUTTON_PRESS_MASK);
  m_GLArea.signal_button_press_event().connect(sigc::mem_fun(*this,&ExampleWindow::glarea_start_pan));

  m_GLArea.add_events(Gdk::BUTTON1_MOTION_MASK | Gdk::BUTTON_PRESS_MASK);
  m_GLArea.signal_motion_notify_event().connect(sigc::mem_fun(*this,&ExampleWindow::glarea_pan));

  m_GLArea.add_events(Gdk::SCROLL_MASK | Gdk::SMOOTH_SCROLL_MASK);
  m_GLArea.signal_scroll_event().connect(sigc::mem_fun(*this,&ExampleWindow::glarea_scroll));


  // Add paned
  m_Paned.set_wide_handle (true);
  m_Box.add(m_Paned);


  show_all_children();
}




std::string ExampleWindow::get_in_tree(std::string section, std::string sub_section) {
  auto children = m_refTreeModel->children();
  for (auto iter = children.begin(), end = children.end(); iter != end; ++iter) {

    auto row = *iter;
    if (row[m_Columns.m_col_name] == section) {

      auto children2 = row->children();
      for (auto iter = children2.begin(), end = children2.end(); iter != end; ++iter) {
        auto crow = *iter;
        if (crow[m_Columns.m_col_name] == sub_section) {
          std::string bc_text = crow.get_value(m_Columns.m_col_text);
          return bc_text;
        }
      }
    }
    
  }

  throw std::invalid_argument( "section " + section + " sub section " + sub_section );
  
  return "";
}



#endif // GTKMM_EXAMPLE_HELLOWORLD_H