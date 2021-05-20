#pragma once

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <Eigen/Core>
#include <string>

namespace squaremind {

namespace iggui = igl::opengl::glfw::imgui;

class MeshGuiMenu : public iggui::ImGuiMenu {
   public:
    void init(igl::opengl::glfw::Viewer* viewer) override;
    bool load(std::string filename) override;

   private:
    void drawMenu();
    bool updateMesh();
    bool updatePath();
    void updateMeshEdges();
    void updateMeshColors();

    void centerCamera();
    void switchOrientation();

    // Mesh
    int mesh_id_ = -1;
    Eigen::MatrixXd mesh_vertices_;
    Eigen::MatrixXi mesh_faces_;
    Eigen::MatrixXd mesh_colors_;
    bool show_mesh_edges_ = false;
    bool use_color_distance_ = false;
    igl::ColorMapType mesh_colormap_ = igl::ColorMapType::COLOR_MAP_TYPE_JET;

    // Shortest path
    int path_id_ = -1;
    int path_source_1_ = 0;
    int path_source_2_ = 1;
    int path_target_ = 7;
    bool path_show_ = false;
};
}  // namespace squaremind
