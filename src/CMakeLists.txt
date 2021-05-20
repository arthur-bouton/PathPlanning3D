file(GLOB
     SOURCES
     "*.cpp"
     "*.h*")
add_executable(main ${SOURCES})
target_link_libraries(main
                      PUBLIC eigen3
                      PUBLIC igl_core
                      PUBLIC igl_opengl_glfw
                      PUBLIC igl_opengl_glfw_imgui)

set_target_properties( main PROPERTIES COMPILE_FLAGS "-Wno-unused-result" )
