#include "mesh_gui_menu.h"

#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/opengl/MeshGL.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <imgui/imgui.h>


#define BASELINE_SOURCE_1  2934
#define BASELINE_SOURCE_2 12439
#define BASELINE_TARGET    5868


namespace PathPlanning3D {

MeshGuiMenu::MeshGuiMenu(std::string filename) {
	filename_ = filename;
}

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

	if (filename_.length() > 0)
		load(filename_);
}

bool MeshGuiMenu::updateMesh() {
    viewer->data(mesh_id_).clear();
    viewer->data(mesh_id_).set_mesh(mesh_vertices_, mesh_faces_);
    updateMeshEdges();
    viewer->data(mesh_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
    viewer->data(mesh_id_).double_sided = true;


    std::cout << "[MeshGuiMenu::updateMesh] Preprocessing mesh..." << std::endl;

	// Initialize the search algorithm with the new vertices and faces:
	path_search_ptr_ = PathSearch::ptr_t( new PathSearch( mesh_vertices_, mesh_faces_ ) );

    std::cout << "[MeshGuiMenu::updateMesh] Mesh ready" << std::endl;


    updateMeshColors();
    return updatePath();
}

void MeshGuiMenu::updateMeshEdges() {
    viewer->data(mesh_id_).show_lines = show_mesh_edges_;
    viewer->data(mesh_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
}

void MeshGuiMenu::updateMeshColors() {
    // reset source/target if outside the mesh vertices
    if (path_source_1_ >= mesh_vertices_.rows()) path_source_1_ = 0;
    if (path_source_2_ >= mesh_vertices_.rows()) path_source_2_ = 0;
    if (path_target_ >= mesh_vertices_.rows()) path_target_ = 0;


	// Define the source vertices:
	Eigen::RowVector3d vertex_source_1 = mesh_vertices_.row( path_source_1_ );
	Eigen::RowVector3d vertex_source_2 = mesh_vertices_.row( path_source_2_ );


	if ( use_color_distance_ )
	{
		// Get the distances of all vertices from the sources:
		Eigen::VectorXd distances = path_search_ptr_->dijkstra_scan( { vertex_source_1, vertex_source_2 }, mesh_vertices_ );

		float distance_max = distances.maxCoeff();

		// Fill the color matrix:
		mesh_colors_ = Eigen::MatrixXd( mesh_vertices_.rows(), 3 );
		for ( std::size_t i = 0 ; i < distances.size() ; ++i )
		{
			float normalized_distance = distances[i]/distance_max;

			// Convert the normalized distance to a rainbow colormap:
			float a = ( 1 - normalized_distance )/0.25;
			int X = floor( a );
			float Y = a - X;
			float R, G, B;
			switch( X )
			{
				case 0: R = 1;   G= Y;   B = 0; break;
				case 1: R = 1-Y; G= 1;   B = 0; break;
				case 2: R = 0;   G= 1;   B = Y; break;
				case 3: R = 0;   G= 1-Y; B = 1; break;
				case 4: R = 0;   G= 0;   B = 1; break;
			}

			mesh_colors_.row( i ) = Eigen::RowVector3d( R, G, B );
		}
	}
	else
		// Otherwise, set all vertices white:
		mesh_colors_ = Eigen::MatrixXd::Constant( 1, 3, 1.0 );


	// Apply the colors:
    viewer->data(mesh_id_).set_colors(mesh_colors_);
    viewer->data(mesh_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
}

bool MeshGuiMenu::updatePath() {
    viewer->data(path_id_).clear();
    // reset source/target if outside the mesh vertices
    if (path_source_1_ >= mesh_vertices_.rows()) path_source_1_ = 0;
    if (path_source_2_ >= mesh_vertices_.rows()) path_source_2_ = 0;
    if (path_target_ >= mesh_vertices_.rows()) path_target_ = 0;


	// Define the starting and goal vertices:
	Eigen::RowVector3d vertex_source_1 = mesh_vertices_.row( path_source_1_ );
	Eigen::RowVector3d vertex_source_2 = mesh_vertices_.row( path_source_2_ );
	Eigen::RowVector3d vertex_target = mesh_vertices_.row( path_target_ );

	// Test the path between points outside the mesh:
	if ( surface_points_ )
	{
		// Points on the surface of the cube:
		vertex_source_1 << 0.3, 0.7, 1.0;
		vertex_target   << 0.7, 0.3, 0.0;
	}


	Eigen::MatrixXd path_vertices = vertex_target;
	if ( path_show_ )
	{
		// Change the metric used to compute the cost of the trajectory:
		auto metric = [=]( Eigen::RowVector3d vertex_1, Eigen::RowVector3d vertex_2 )
		{
			Eigen::RowVector3d vec = vertex_2 - vertex_1;
			// Add a penalty based on the deviation of the trajectory relative to the horizontal origin plane:
			return vec.norm() + fabs( vertex_2[2] )*metric_factor_;
		};
		path_search_ptr_->setMetric( metric );

		// Search for the shortest path:
		if ( surface_points_ )
			path_vertices = path_search_ptr_->A_star_search( { vertex_source_1 }, vertex_target );
		else
			path_vertices = path_search_ptr_->A_star_search( { vertex_source_1, vertex_source_2 }, vertex_target );
	}


	// Fill the edge matrix:
	Eigen::MatrixXi path_edges( path_vertices.rows() - 1, 2 );
	for ( std::size_t i = 0 ; i < path_vertices.rows() - 1 ; ++i )
		path_edges.row( i ) << i, i + 1;

	// Color of the path:
	const Eigen::RowVector3d path_colors( 1.0, 0.0, 0.0 );

	// Draw the path:
	viewer->data( path_id_ ).set_edges( path_vertices, path_edges, path_colors );
	if ( surface_points_ )
		viewer->data( path_id_ ).line_width = 10.0f;
	else
		viewer->data( path_id_ ).line_width = 10.0f;
	viewer->data( path_id_ ).is_visible = path_show_;
	viewer->data( path_id_ ).dirty = igl::opengl::MeshGL::DIRTY_ALL;

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
	filename_ = filename;
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
			surface_points_ = false;
            updatePath();
			updateMeshColors();
        }
        // This slider allows to select the path's second source vertex
        if (ImGui::SliderInt("Source 2", &path_source_2_, 0,
                             mesh_vertices_.rows() - 1)) {
			surface_points_ = false;
            updatePath();
			updateMeshColors();
        }
        // This slider allows to select the path's target vertex
        if (ImGui::SliderInt("Target", &path_target_, 0,
                             mesh_vertices_.rows() - 1)) {
			surface_points_ = false;
            updatePath();
        }
        // Toggle visibility of computed path on mesh
        if (ImGui::Checkbox("Show path", &path_show_)) {
            updatePath();
            viewer->data(path_id_).is_visible = path_show_;
            viewer->data(path_id_).dirty = igl::opengl::MeshGL::DIRTY_ALL;
        }

		// Test the path between points outside the mesh:
		if ( ImGui::Button( "Surface-point test", ImVec2( -1, 0 ) ) )
		{
			surface_points_ = true;
            updatePath();
			updateMeshColors();
		}
        // Set baseline source and target points to be able to compare the running time:
		if ( ImGui::Button( "Baseline test", ImVec2( -1, 0 ) ) )
		{
			surface_points_ = false;
			path_source_1_ = BASELINE_SOURCE_1;
			path_source_2_ = BASELINE_SOURCE_2;
			path_target_   = BASELINE_TARGET;
            updatePath();
			updateMeshColors();
		}
		// Slider changing the metric used to compute the optimal path:
        if ( ImGui::SliderFloat( "Metric factor", &metric_factor_, 0, 1 ) )
		{
            updatePath();
			updateMeshColors();
		}
    }
}

}  // namespace PathPlanning3D
