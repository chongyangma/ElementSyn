cmake_minimum_required(VERSION 2.8.12)
project(DynamicElement)

file(GLOB DET_SOURCES "*.cpp")
file(GLOB DET_HEADERS "*.h")

set(LODEPNG_DIR "../../3rdparty/lodepng")
list(APPEND DET_SOURCES ${LODEPNG_DIR}/lodepng.cpp)
list(APPEND DET_HEADERS ${LODEPNG_DIR}/lodepng.h)

link_directories(${PROJECT_BINARY_DIR}/../../lib)

add_library(DynamicElement STATIC ${DET_SOURCES} ${DET_HEADERS})

target_link_libraries(DynamicElement ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})

if(NOT WIN32)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libdet.a" DESTINATION ${CMAKE_SOURCE_DIR}/lib)
endif(NOT WIN32)
