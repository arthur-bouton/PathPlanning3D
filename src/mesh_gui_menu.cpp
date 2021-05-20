#include "mesh_gui_menu.h"

#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/opengl/MeshGL.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <imgui/imgui.h>

namespace squaremind {

void MeshGuiMenu::init(igl::opengl::glfw::Viewer* viewer) {
    ImGuiMenu::init(viewer);
    callback_draw_viewer_menu = [ this ](){ drawMenu(); };

    // Inline mesh of a cube
    mesh_vertices_ = (Eigen::MatrixXd(24, 3) <<
         0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
         0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
         0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0,
         1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0)
            .finished();

    mesh_faces_ = (Eigen::MatrixXi(12, 3) <<
        0, 1, 2, 0, 2, 3, 4, 6, 5, 4, 7, 6,
        8, 9, 10, 8, 10, 11, 12, 14, 13, 12, 15, 14,
        16, 17, 18, 16, 18, 19, 20, 22, 21, 20, 23, 22)
            .finished();

    mesh_colors_ = Eigen::MatrixXd::Constant(mesh_vertices_.rows(), 3, 1.0);
    mesh_id_ = viewer->append_mesh(true);
    path_id_ = viewer->append_mesh(true);
    updateMesh();
}

bool MeshGuiMenu::updateMesh() {
    viewer->data(mesh_id_).clear();
    viewer->data(mesh_id_).set_mesh(mesh_vertices_, mesh_faces_);
    updateMeshEdges();
    updateMeshColors();
    viewer->data(mesh_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
    viewer->data(mesh_id_).double_sided = true;
    return updatePath();
}

void MeshGuiMenu::updateMeshEdges() {
    viewer->data(mesh_id_).show_lines = show_mesh_edges_;
    viewer->data(mesh_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
}

void MeshGuiMenu::updateMeshColors() {
    // this function can be used to add color on the model w.r.t distances to
    // source points
    mesh_colors_ = Eigen::MatrixXd::Constant(1, 3, 1.0);
    viewer->data(mesh_id_).set_colors(mesh_colors_);
    viewer->data(mesh_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
}

bool MeshGuiMenu::updatePath() {
    viewer->data(path_id_).clear();
    // reset source/target if outside the mesh vertices
    if (path_source_1_ >= mesh_vertices_.rows()) path_source_1_ = 0;
    if (path_source_2_ >= mesh_vertices_.rows()) path_source_2_ = 0;
    if (path_target_ >= mesh_vertices_.rows()) path_target_ = 0;

    /** Add your path computation function calls here
     *
     *
     *
     *
     *
     **/
    // updateMeshColors();

    // Quick example to draw the path :
    // 1. Update a matrix with the set of vertices visited by your path
    // Eigen::MatrixXd path_vertices(vertices_in_my_path.size(), 3);
    // path_vertices = ...
    // 2. Define the connectivity between each vertices to create edges
    // Eigen::MatrixXi path_edges(vertices.size() - 1, 2);
    // path_edges = ...
    // 3. Set a color to your trajectory
    // Eigen::MatrixXd path_colors = ... 1.0, 0.0, 0.0 ... (red)
    // 4. Set path in viewer
    // viewer->data(path_id_).set_edges(path_vertices, path_edges,
    //                                      path_colors);
    // viewer->data(path_id_).line_width = 2.0f;
    // 5. Draw
    // viewer->data(path_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
    return true;
}

void MeshGuiMenu::centerCamera() {
    viewer->core().camera_translation = Eigen::Vector3f { 0.0f, 0.0f, 0.0f };
    viewer->core().camera_zoom = 1.0f;
    viewer->core().trackball_angle = Eigen::Quaternionf::Identity();
    viewer->core().align_camera_center(viewer->data(mesh_id_).V);
}

void MeshGuiMenu::switchOrientation() {
    for (std::size_t i = 0; i < mesh_vertices_.rows(); ++i) {
        auto v = Eigen::Vector3d { mesh_vertices_.row(i) };
        mesh_vertices_.row(i) = Eigen::Vector3d { v.y(), v.z(), v.x() };
    }
    updateMesh();
}


bool MeshGuiMenu::load(std::string filename) {
    viewer->data(mesh_id_).clear();
    igl::readPLY(filename, mesh_vertices_, mesh_faces_);
    std::cout << "[MeshGuiMenu::load] Loaded model: " << filename << std::endl;
    return updateMesh();
}

void MeshGuiMenu::drawMenu() {
    // Helper for setting viewport specific mesh options
    auto make_checkbox = [&](const char* label, unsigned int& option) {
        return ImGui::Checkbox(
            label, [&]() { return viewer->core().is_set(option); },
            [&](bool value) { return viewer->core().set(option, value); });
    };

    // Draw option
    if (ImGui::CollapsingHeader("Draw Options",
                                ImGuiTreeNodeFlags_DefaultOpen)) {
        make_checkbox("Show overlay", viewer->data().show_overlay);
    }

    // Mesh
    if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen)) {
        float w = ImGui::GetContentRegionAvailWidth();
        float p = ImGui::GetStyle().FramePadding.x;
        if (ImGui::Button("Load##Mesh", ImVec2(-1, 0))) {
            viewer->open_dialog_load_mesh();
            centerCamera();
        }
        // Draw option
        if (ImGui::Button("Center object", ImVec2(-1, 0))) {
            centerCamera();
        }
        if (ImGui::Button("Switch orientation", ImVec2(-1, 0))) {
            switchOrientation();
            centerCamera();
        }
        if (ImGui::Checkbox("Show mesh edges", &show_mesh_edges_)) {
            updateMeshEdges();
        }
        if (ImGui::Checkbox("Use distance color", &use_color_distance_)) {
            updateMeshColors();
        }
    }

    // Path
    if (ImGui::CollapsingHeader("Path", ImGuiTreeNodeFlags_DefaultOpen)) {
        // This slider allows to select the path's first source vertex
        if (ImGui::SliderInt("Source 1", &path_source_1_, 0,
                             mesh_vertices_.rows() - 1)) {
            updatePath();
        }
        // This slider allows to select the path's second source vertex
        if (ImGui::SliderInt("Source 2", &path_source_2_, 0,
                             mesh_vertices_.rows() - 1)) {
            updatePath();
        }
        // This slider allows to select the path's target vertex
        if (ImGui::SliderInt("Target", &path_target_, 0,
                             mesh_vertices_.rows() - 1)) {
            updatePath();
        }
        // Toggle visibility of computed path on mesh
        if (ImGui::Checkbox("Show path", &path_show_)) {
            viewer->data(path_id_).is_visible = path_show_;
            viewer->data(path_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
        }
        // You can add here any other elements to configure your algorihtm
    }
}

}  // namespace squaremind
