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
  source/CAT/experimental_point.h
  source/CAT/experimental_double.h
  source/CAT/CAT_interface.h
  source/CAT/clusterizer.h
  source/CAT/cell_base.h
  source/CAT/calorimeter_hit.h
  # source/TrackFit/fit_utils.h
  # source/TrackFit/gg_hit.h
  # source/TrackFit/helix_fit_mgr.h
  # source/TrackFit/i_drift_time_calibration.h
  # source/TrackFit/line_fit_mgr.h
  )

list(APPEND LCAT_SOURCES
  source/CAT/experimental_point.cc
  source/CAT/experimental_double.cc
  source/CAT/CAT_interface.cc
  source/CAT/clusterizer.cc
  source/CAT/cell_base.cc
  source/CAT/calorimeter_hit.cc
  # source/TrackFit/fit_utils.cc
  # source/TrackFit/gg_hit.cc
  # source/TrackFit/helix_fit_mgr.cc
  # source/TrackFit/i_drift_time_calibration.cc
  # source/TrackFit/line_fit_mgr.cc
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
option(FalaiseLCATPlugin_ENABLE_TESTING "Build unit testing system for FalaiseLCATPlugin" ON)
if(FalaiseLCATPlugin_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(testing)
endif()
