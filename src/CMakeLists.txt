add_definitions(-DGL_SILENCE_DEPRECATION=1)

find_package(OpenMP)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIR})

include_directories(../3rdparty/eigen)

add_subdirectory(DynamicElement)
add_subdirectory(ParticleSystem)
add_subdirectory(TreeBranches)
