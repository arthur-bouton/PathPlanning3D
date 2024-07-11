option(GLFW_BUILD_EXAMPLES "Build the GLFW example programs" OFF)
option(GLFW_BUILD_TESTS "Build the GLFW test programs" OFF)
option(GLFW_BUILD_DOCS "Build the GLFW documentation" OFF)
option(GLFW_INSTALL "Generate installation target" OFF)
option(GLFW_VULKAN_STATIC
       "Assume the Vulkan loader is linked with the application" OFF)
option(GLFW_USE_OSMESA "Use OSMesa for offscreen context creation" OFF)
option(GLFW_USE_HYBRID_HPG "Force use of high-performance GPU on hybrid systems"
       OFF)
option(GLFW_USE_WAYLAND "Use Wayland for window creation" OFF)
option(USE_MSVC_RUNTIME_LIBRARY_DLL "Use MSVC runtime library DLL" ON)

add_subdirectory(glfw)
set_property(TARGET glfw PROPERTY FOLDER ${CMAKE_FOLDER})
