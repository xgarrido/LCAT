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
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/CAT)

############################################################################################
# - LCAT library:

list(APPEND LCAT_HEADERS
  source/CAT/CATAlgorithm/broken_line.h
  source/CAT/CATAlgorithm/calorimeter_hit.h
  source/CAT/CATAlgorithm/CAT_interface.h
  source/CAT/CATAlgorithm/cell_base.h
  source/CAT/CATAlgorithm/cell_couplet.h
  source/CAT/CATAlgorithm/cell_triplet.h
  source/CAT/CATAlgorithm/circle_base.h
  source/CAT/CATAlgorithm/CircleRegression.h
  source/CAT/CATAlgorithm/clockable.h
  source/CAT/CATAlgorithm/Clock.h
  source/CAT/CATAlgorithm/cluster.h
  source/CAT/CATAlgorithm/clusterizer.h
  source/CAT/CATAlgorithm/cylinder.h
  source/CAT/CATAlgorithm/experimental_double.h
  source/CAT/CATAlgorithm/experimental_point.h
  source/CAT/CATAlgorithm/experimental_vector.h
  source/CAT/CATAlgorithm/handle.h
  source/CAT/CATAlgorithm/helix.h
  source/CAT/CATAlgorithm/i_predicate.h
  source/CAT/CATAlgorithm/joint.h
  source/CAT/CATAlgorithm/LinearRegression.h
  source/CAT/CATAlgorithm/line.h
  source/CAT/CATAlgorithm/lt_utils.h
  source/CAT/CATAlgorithm/node.h
  source/CAT/CATAlgorithm/plane.h
  source/CAT/CATAlgorithm/scenario.h
  source/CAT/CATAlgorithm/sequence_base.h
  source/CAT/CATAlgorithm/Sequence.h
  source/CAT/CATAlgorithm/sequentiator.h
  source/CAT/CATAlgorithm/tracked_data_base.h
  source/CAT/CATAlgorithm/tracking_object.h
  source/CAT/CATAlgorithm/utilities.h
  )

list(APPEND LCAT_SOURCES
  source/CAT/CATAlgorithm/broken_line.cpp
  source/CAT/CATAlgorithm/calorimeter_hit.cpp
  source/CAT/CATAlgorithm/CAT_interface.cpp
  source/CAT/CATAlgorithm/cell_base.cpp
  source/CAT/CATAlgorithm/cell_couplet.cpp
  source/CAT/CATAlgorithm/cell_triplet.cpp
  source/CAT/CATAlgorithm/circle_base.cpp
  source/CAT/CATAlgorithm/CircleRegression.cpp
  source/CAT/CATAlgorithm/clockable.cpp
  source/CAT/CATAlgorithm/Clock.cpp
  source/CAT/CATAlgorithm/cluster.cpp
  source/CAT/CATAlgorithm/clusterizer.cpp
  source/CAT/CATAlgorithm/experimental_double.cpp
  source/CAT/CATAlgorithm/experimental_point.cpp
  source/CAT/CATAlgorithm/experimental_vector.cpp
  source/CAT/CATAlgorithm/joint.cpp
  source/CAT/CATAlgorithm/LinearRegression.cpp
  source/CAT/CATAlgorithm/line.cpp
  source/CAT/CATAlgorithm/node.cpp
  source/CAT/CATAlgorithm/plane.cpp
  source/CAT/CATAlgorithm/scenario.cpp
  source/CAT/CATAlgorithm/sequence_base.cpp
  source/CAT/CATAlgorithm/sequentiator.cpp
  source/CAT/CATAlgorithm/tracking_object.cpp
  )

# Build a dynamic library from our sources
add_library(LCAT SHARED ${LCAT_HEADERS} ${LCAT_SOURCES})

target_link_libraries(LCAT ${Bayeux_LIBRARIES})

# Apple linker requires dynamic lookup of symbols, so we
# add link flags on this platform
if(APPLE)
  set_target_properties(LCAT
    PROPERTIES LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

# Install it:
install(TARGETS LCAT DESTINATION ${CMAKE_INSTALL_LIBDIR}/${FALAISE_PLUGINLIBDIR})


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
install(TARGETS Falaise_LCAT DESTINATION ${CMAKE_INSTALL_LIBDIR}/${FALAISE_PLUGINLIBDIR})

# Test support:
# option(FalaiseLCATPlugin_ENABLE_TESTING "Build unit testing system for FalaiseLCATPlugin" ON)
# if(FalaiseLCATPlugin_ENABLE_TESTING)
#   enable_testing()
#   add_subdirectory(testing)
# endif()
