set(IGL_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/libigl/include")

set(LIBIGL_STATIC_LIBRARY 0)

# Core of IGL
if(LIBIGL_STATIC_LIBRARY)
  file(GLOB
       SOURCES_IGL_CORE
       "${IGL_SOURCE_DIR}/igl/*.cpp"
       "${IGL_SOURCE_DIR}/igl/*.h*")
  add_library(igl_core STATIC ${SOURCES_IGL_CORE})
  target_compile_definitions(igl_core PUBLIC IGL_STATIC_LIBRARY)
  target_compile_features(igl_core PUBLIC cxx_std_11)
  target_include_directories(igl_core PUBLIC ${IGL_SOURCE_DIR})
  target_link_libraries(igl_core PUBLIC Eigen3::Eigen Threads::Threads)
else()
  add_library(igl_core INTERFACE)
  target_compile_features(igl_core INTERFACE cxx_std_11)
  target_include_directories(igl_core INTERFACE ${IGL_SOURCE_DIR})
  target_link_libraries(igl_core INTERFACE Eigen3::Eigen Threads::Threads)
endif()

# Modules
macro(add_igl_module name)
  set(target igl_${name})
  if(LIBIGL_STATIC_LIBRARY)
    string(REPLACE "_"
                   "/"
                   path
                   ${name})
    file(GLOB
         SOURCES_IGL_${name}
         "${IGL_SOURCE_DIR}/igl/${path}/*.cpp"
         "${IGL_SOURCE_DIR}/igl/${path}/*.h*")
    add_library(${target} STATIC ${SOURCES_IGL_${name}})
    target_link_libraries(${target} PUBLIC igl_core ${ARGN})
  else()
    add_library(${target} INTERFACE)
    target_link_libraries(${target} INTERFACE igl_core ${ARGN})
  endif()
endmacro()

add_igl_module(opengl OpenGL::GL glad)
add_igl_module(opengl_glfw igl_opengl glfw)
add_igl_module(opengl_glfw_imgui igl_opengl_glfw imgui)
add_igl_module(triangle triangle)
add_igl_module(png igl_stb_image)
