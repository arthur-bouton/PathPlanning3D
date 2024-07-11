#include <igl/opengl/glfw/Viewer.h>

#include <iostream>

#include "mesh_gui_menu.h"

int main(int argc, char* argv[]) {
  igl::opengl::glfw::Viewer viewer;

  PathPlanning3D::MeshGuiMenu gui_menu(argc > 1 ? argv[1] : "");
  viewer.plugins.push_back(&gui_menu);
  viewer.launch();
}
