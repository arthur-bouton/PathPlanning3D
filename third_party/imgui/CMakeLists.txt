set(IMGUI_DIR ${CMAKE_CURRENT_LIST_DIR}/imgui)
file(GLOB
     IMGUI_SOURCES
     ${IMGUI_DIR}/*.h
     ${IMGUI_DIR}/*.cpp)
add_library(imgui STATIC ${IMGUI_SOURCES}
                         ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
                         ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp)
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD)
target_include_directories(imgui
                           PUBLIC ${IMGUI_DIR}
                           PUBLIC ${IMGUI_DIR}/..
                           PUBLIC ${IMGUI_DIR}/backends)
target_link_libraries(imgui PRIVATE glad PRIVATE glfw)
