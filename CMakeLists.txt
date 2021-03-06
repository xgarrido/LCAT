# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(FalaiseLCATPlugin VERSION "2.1.0")

# Modules use Falaise, so we need to locate this or fail
# find_package(Falaise REQUIRED)

# Ensure our code can see the Falaise headers
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)

############################################################################################
# - LCAT library:

list(APPEND LCAT_HEADERS
  source/CAT/broken_line.h
  source/CAT/calorimeter_hit.h
  source/CAT/CAT_interface.h
  source/CAT/cell.h
  source/CAT/cell_couplet.h
  source/CAT/cell_triplet.h
  source/CAT/circle_base.h
  source/CAT/CircleRegression.h
  source/CAT/clockable.h
  source/CAT/Clock.h
  source/CAT/cluster.h
  source/CAT/clusterizer.h
  source/CAT/experimental_double.h
  source/CAT/experimental_point.h
  source/CAT/experimental_vector.h
  source/CAT/helix.h
  source/CAT/joint.h
  source/CAT/LinearRegression.h
  source/CAT/line.h
  source/CAT/node.h
  source/CAT/plane.h
  source/CAT/scenario.h
  source/CAT/sequence_base.h
  source/CAT/sequentiator.h
  source/CAT/tracked_data.h
  source/CAT/utilities.h
  )

list(APPEND LCAT_SOURCES
  source/CAT/broken_line.cpp
  source/CAT/calorimeter_hit.cpp
  source/CAT/CAT_interface.cpp
  source/CAT/cell.cc
  source/CAT/cell_couplet.cpp
  source/CAT/cell_triplet.cpp
  source/CAT/circle_base.cpp
  source/CAT/CircleRegression.cpp
  source/CAT/clockable.cpp
  source/CAT/Clock.cpp
  source/CAT/cluster.cpp
  source/CAT/clusterizer.cc
  source/CAT/experimental_double.cpp
  source/CAT/experimental_point.cpp
  source/CAT/experimental_vector.cpp
  source/CAT/joint.cpp
  source/CAT/LinearRegression.cpp
  source/CAT/line.cpp
  source/CAT/node.cpp
  source/CAT/plane.cpp
  source/CAT/scenario.cpp
  source/CAT/sequence_base.cpp
  source/CAT/sequentiator.cpp
  source/CAT/tracked_data.cc
  source/CAT/utilities.cc
  )

# Build a dynamic library from our sources
add_library(LCAT SHARED ${LCAT_HEADERS} ${LCAT_SOURCES})

target_link_libraries(LCAT Bayeux::Bayeux)

# Apple linker requires dynamic lookup of symbols, so we
# add link flags on this platform
if(APPLE)
  set_target_properties(LCAT
    PROPERTIES LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

# Install it:
install(TARGETS LCAT DESTINATION ${CMAKE_INSTALL_PLUGINDIR})


############################################################################################
# - LCAT modules:
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/falaise)

# - Headers:
list(APPEND FalaiseLCATPlugin_HEADERS
  source/falaise/snemo/reconstruction/cat_driver.h
  source/falaise/snemo/reconstruction/cat_tracker_clustering_module.h
  )

# - Sources:
list(APPEND FalaiseLCATPlugin_SOURCES
  source/falaise/snemo/reconstruction/cat_driver.cc
  source/falaise/snemo/reconstruction/cat_tracker_clustering_module.cc
  )

############################################################################################

# Build a dynamic library from our sources
add_library(Falaise_LCAT SHARED
  ${FalaiseLCATPlugin_HEADERS}
  ${FalaiseLCATPlugin_SOURCES})

target_link_libraries(Falaise_LCAT LCAT Falaise)

# Apple linker requires dynamic lookup of symbols, so we
# add link flags on this platform
if(APPLE)
  set_target_properties(Falaise_LCAT
    PROPERTIES
    LINK_FLAGS "-undefined dynamic_lookup"
    INSTALL_RPATH "@loader_path"
    )
endif()

# Install it:
install(TARGETS Falaise_LCAT DESTINATION ${CMAKE_INSTALL_PLUGINDIR})

# Test support:
option(FalaiseLCATPlugin_ENABLE_TESTING "Build unit testing system for FalaiseLCATPlugin" ON)
if(FalaiseLCATPlugin_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(testing)
endif()
