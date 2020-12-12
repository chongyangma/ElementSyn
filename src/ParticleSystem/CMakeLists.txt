cmake_minimum_required(VERSION 2.8.12)
project(ParticleSystem)

file(GLOB PARTICLE_SYSTEM_SOURCES "*.cpp")
file(GLOB PARTICLE_SYSTEM_HEADERS "*.h")

add_executable(ParticleSystem ${PARTICLE_SYSTEM_SOURCES} ${PARTICLE_SYSTEM_HEADERS})

target_link_libraries(ParticleSystem ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} DynamicElement)

if(NOT WIN32)
    find_library(OPENGL_glu_LIBRARY GLU PATHS /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
    target_link_libraries(ParticleSystem ${OPENGL_glu_LIBRARY})
endif()
