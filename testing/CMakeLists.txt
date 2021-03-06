# - List of test programs:
set(FalaiseLCATPlugin_TESTS
  # test_cat_driver.cxx
  # test_cat_tracker_clustering_module.cxx
  test_gg_cell.cxx
  )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(testing_SOURCES utilities.cc)

set(_cat_TEST_ENVIRONMENT "FALAISE_RESOURCE_DIR=${FALAISE_DATAROOTDIR}/resources")

if(FalaiseLCATPlugin_ENABLE_TESTING)
  foreach(_testsource ${FalaiseLCATPlugin_TESTS})
    get_filename_component(_testname "${_testsource}" NAME_WE)
    set(_testname "falaiselcatplugin-${_testname}")
    add_executable(${_testname} ${_testsource} ${testing_SOURCES})
    target_link_libraries(${_testname} Falaise_LCAT Falaise)
    # - On Apple, ensure dynamic_lookup of undefined symbols
    if(APPLE)
      set_target_properties(${_testname} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()
    add_test(NAME ${_testname} COMMAND ${_testname})
    set_property(TEST ${_testname}
      APPEND PROPERTY ENVIRONMENT ${_cat_TEST_ENVIRONMENT}
      )
    # - For now, dump them into the testing output directory
    set_target_properties(${_testname}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fltests/modules
      ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fltests/modules
      )
  endforeach()
endif()
