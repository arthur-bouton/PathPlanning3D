#include <igl/opengl/glfw/Viewer.h>

#include <iostream>

#include "mesh_gui_menu.h"

int main(int /*argc*/, char* /*argv*/[]) {
  igl::opengl::glfw::Viewer viewer;

  PathPlanning3D::MeshGuiMenu gui_menu;
  viewer.plugins.push_back(&gui_menu);
  viewer.launch();
}
