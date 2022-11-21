/*

Compile with
g++ main.cpp gui.cpp -o main `pkg-config --cflags --libs gtkmm-3.0` -std=c++17

*/


#include <gui/gui.h>
#include <gtkmm.h>
#include <iostream>




int main(int argc, char *argv[])
{

  auto app = Gtk::Application::create(argc, argv, "org.gtkmm.example");

  ExampleWindow window(app);

  //Shows the window and returns when it is closed.
  return app->run(window);
}

