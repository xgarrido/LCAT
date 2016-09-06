// Standard library:
#include <iostream>
#include <exception>
#include <cstdlib>
#include <string>

// Third party:
// - Bayeux/datatools:
#include <datatools/logger.h>
#include <datatools/properties.h>
#include <datatools/utils.h>
#include <datatools/clhep_units.h>
// - Bayeux/geomtools:
#include <geomtools/manager.h>

// Falaise:
#include <falaise/falaise.h>
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/geometry/locator_plugin.h>
#include <falaise/snemo/geometry/gg_locator.h>

// This project:
#include <snemo/reconstruction/cat_driver.h>
#include <snemo/reconstruction/cat_tracker_clustering_module.h>

// Testing resources:
#include <utilities.h>

int main(int argc_, char ** argv_)
{
  falaise::initialize(argc_, argv_);
  int error_code = EXIT_SUCCESS;
  datatools::logger::priority logging = datatools::logger::PRIO_NOTICE;
  try {
    std::clog << "Test program for class 'snemo::reconstruction::cat_tracker_clustering_module'!" << std::endl;

    bool draw = false;
    int iarg = 1;
    while (iarg < argc_) {
      std::string token = argv_[iarg];
      if (token[0] == '-') {
        std::string option = token;
        if ((option == "-D") || (option == "--draw")) {
          draw = true;
        } else {
          DT_LOG_WARNING(logging, "Ignoring option '" << option << "'!");
        }
      } else {
        DT_LOG_WARNING(logging, "Ignoring argument '" << token << "'!");
      }
      iarg++;
    }

    srand48(314159);

        // Geometry manager:
    std::string GeoConfigFile = "@falaise:config/snemo/demonstrator/geometry/4.0/manager.conf";
    datatools::fetch_path_with_env(GeoConfigFile);
    datatools::properties GeoConfig;
    datatools::properties::read_config(GeoConfigFile, GeoConfig);
    geomtools::manager Geo;
    Geo.initialize(GeoConfig);

    // Extract Geiger locator:
    const snemo::geometry::gg_locator * gg_locator = 0;
    std::string locator_plugin_name = "locators_driver";
    if (Geo.has_plugin(locator_plugin_name)
        && Geo.is_plugin_a<snemo::geometry::locator_plugin>(locator_plugin_name)) {
      DT_LOG_NOTICE(logging, "Found locator plugin named '" << locator_plugin_name << "'");
      const snemo::geometry::locator_plugin & lp
        = Geo.get_plugin<snemo::geometry::locator_plugin>(locator_plugin_name);
      // Set the Geiger cell locator :
      gg_locator = &(lp.get_gg_locator());
    }

    // The CAT driver:
    // Parameters for the CAT driver:
    datatools::properties CATconfig;
    // CATconfig.store_string("CAT.level",           "normal");
    // CATconfig.store_real("CAT.max_time",          5000.0 * CLHEP::ms); // microsec ?
    // CATconfig.store_real("CAT.small_radius",      2.0 * CLHEP::mm);
    // CATconfig.store_real("CAT.probmin",           0.0);
    // CATconfig.store_integer("CAT.nofflayers",     1);
    // CATconfig.store_integer("CAT.first_event",    -1);
    // CATconfig.store_real("CAT.ratio",             10000.0);
    snemo::reconstruction::cat_driver CAT;
    CAT.set_logging_priority(logging);
    CAT.set_geometry_manager(Geo);
    CAT.initialize(CATconfig);

    // The CAT module:
    snemo::reconstruction::cat_tracker_clustering_module CATmod;
    CATmod.set_logging_priority(logging);
    CATmod.set_geometry_manager(Geo);
    CATmod.initialize_standalone(CATconfig);

    // Event loop:
    for (size_t i = 0; i < 3; i++) {
      DT_LOG_NOTICE(logging, "Processing event #" << i);
      datatools::things eventRecord;
      snemo::datamodel::calibrated_data & CD
        = eventRecord.add<snemo::datamodel::calibrated_data>(snemo::datamodel::data_info::default_calibrated_data_label());
      snemo::datamodel::calibrated_data::tracker_hit_collection_type & gghits = CD.calibrated_tracker_hits();
      generate_gg_hits(*gg_locator, gghits);
      dpp::base_module::process_status status = CATmod.process(eventRecord);
      if (status != 0) {
        DT_LOG_FATAL(logging, "CAT module failed!");
        break;
      }
      const snemo::datamodel::tracker_clustering_data & TCD
        = eventRecord.get<snemo::datamodel::tracker_clustering_data>(snemo::datamodel::data_info::default_tracker_clustering_data_label());
      TCD.tree_dump(std::clog, "Tracker clustering data: ", "[trace]: ");
      if (draw) display_event(*gg_locator, gghits, TCD);
      eventRecord.tree_dump(std::clog, "Event record: ", "[trace]: ");
    }

    // Terminate the CAT module:
    CATmod.reset();

    DT_LOG_NOTICE(logging, "The end.");
  }
  catch (std::exception & error) {
    DT_LOG_FATAL(logging, error.what());
    error_code = EXIT_FAILURE;
  }
  catch (...) {
    DT_LOG_FATAL(logging, "Unexpected error!");
    error_code = EXIT_FAILURE;
  }
  falaise::terminate();
  return error_code;
}
