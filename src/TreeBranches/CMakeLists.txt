cmake_minimum_required(VERSION 2.8.12)
project(TreeBranches)

file(GLOB TREE_BRANCHES_SOURCES "*.cpp")
file(GLOB TREE_BRANCHES_HEADERS "*.h")

add_executable(TreeBranches ${TREE_BRANCHES_SOURCES} ${TREE_BRANCHES_HEADERS})

target_link_libraries(TreeBranches ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} DynamicElement)

if(NOT WIN32)
    find_library(OPENGL_glu_LIBRARY GLU PATHS /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
    target_link_libraries(TreeBranches ${OPENGL_glu_LIBRARY})
endif()
