<p align="center">
	<img src="snapshot.png?raw=true" height="500">
</p>

Efficient implementation of path planning algorithms in 3D in C++. The shortest path is computed in real time and displayed in red on the 3D mesh.


## Installation

First, you'll need OpenGL and Eigen:

	sudo apt install xorg-dev libglu1-mesa-dev freeglut3-dev zenity libeigen3-dev

Then, build the application:

	mkdir build && cd build && cmake .. && make


## Usage

Run the application with:

	./main

By default, a simple cube is displayed, whose vertices and faces are defined in [src/mesh_gui_menu.cpp](src/mesh_gui_menu.cpp). To display another mesh, click on `Load` in the top left menu and select a PLY file (examples can be found in [meshes](meshes)). Otherwise, you can open a mesh directly by specifying the file as an argument:

	./main ../meshes/vase-lion.ply

Move the source and target cursors in the top left menu to change the start and goal vertices of the path and see the resulting path being updated in real time.
