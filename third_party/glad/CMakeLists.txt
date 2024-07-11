add_library(glad STATIC src/glad.c include/glad/glad.h
                        include/KHR/khrplatform.h)
target_include_directories(glad PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(glad PUBLIC ${CMAKE_DL_LIBS})
