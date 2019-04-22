cmake_minimum_required(VERSION 2.8)
project (k_amplitudes)

# Dir for my own files
include_directories(include)
include_directories(belle2018)
include_directories(evtgen)
include_directories(HepMC-2.06.08)
include_directories(tclap/include)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO)
include_directories(${ROOT_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

MESSAGE(${ROOT_INCLUDE_DIRS})
MESSAGE(${ROOT_CXX_FLAGS})

# Suppress a warning on MacOSX
set(CMAKE_MACOSX_RPATH 1)

# Compile the HepMC module (the CMakeLists have been modified so only needed lib is made)
# necessary options:
set(momentum GEV)
set(length MM)
add_subdirectory(HepMC-2.06.08)


# Files from EvtGen
# First build the needed EvtGen libraries
file(GLOB SOURCES_EVTGENBASE "evtgen/src/EvtGenBase/*.cpp")
add_library(EvtGenBaselib SHARED
    ${SOURCES_EVTGENBASE}
)
target_link_libraries(EvtGenBaselib HepMC)

# Then build main program
add_executable(k_amplitudes 
    src/main.cpp
    src/DtoKpipiAmplitude.cpp
    src/EvtGenAmplitude.cpp
    src/OutputMaker.cpp
    src/TDalitz.cpp
)

link_directories(${ROOTSYS}/lib)
target_link_libraries(k_amplitudes 
    EvtGenBaselib
    ${ROOT_LIBRARIES}
)

INSTALL(PROGRAMS build/k_amplitudes DESTINATION ${CMAKE_SOURCE_DIR})

