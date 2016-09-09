/// \file falaise/snemo/reconstruction/clustering_comparison_module.cc

// Ourselves:
#include <snemo/reconstruction/clustering_comparison_module.h>

// Standard library:
#include <stdexcept>
#include <sstream>

// This project:
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>

namespace snemo {

  namespace reconstruction {

    // Registration instantiation macro :
    DPP_MODULE_REGISTRATION_IMPLEMENT(clustering_comparison_module,
                                      "snemo::reconstruction::clustering_comparison_module")

    void clustering_comparison_module::_set_defaults()
    {
      _TCD_labels_.clear();
      return;
    }

    void clustering_comparison_module::initialize(const datatools::properties  & setup_,
                                                  datatools::service_manager   & /* service_manager_ */,
                                                  dpp::module_handle_dict_type & /* module_dict_ */)
    {
      DT_THROW_IF(is_initialized(), std::logic_error,
                  "Module '" << get_name() << "' is already initialized ! ");

      dpp::base_module::_common_initialize(setup_);

      if (_TCD_labels_.empty()) {
        if (setup_.has_key("TCD_labels")) {
          setup_.fetch("TCD_labels", _TCD_labels_);
        }
      }
      _set_initialized(true);
      return;
    }

    void clustering_comparison_module::reset()
    {
      DT_THROW_IF(! is_initialized(), std::logic_error,
                  "Module '" << get_name() << "' is not initialized !");
      _set_initialized(false);
      _set_defaults();
      return;
    }

    // Constructor :
    clustering_comparison_module::clustering_comparison_module(datatools::logger::priority logging_priority_)
      : dpp::base_module(logging_priority_)
    {
      _set_defaults();
      return;
    }

    // Destructor :
    clustering_comparison_module::~clustering_comparison_module()
    {
      if (is_initialized()) clustering_comparison_module::reset();
      return;
    }

    // Processing :
    dpp::base_module::process_status clustering_comparison_module::process(datatools::things & data_record_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");
      DT_THROW_IF(! is_initialized(), std::logic_error,
                  "Module '" << get_name() << "' is not initialized !");

      auto exit_message = [] (auto & bank1, auto & bank2)
        {
          DT_LOG_WARNING(datatools::logger::PRIO_ALWAYS,
                         "'" << bank1 << "' and '" << bank2 << "' banks are different !");
          return dpp::base_module::PROCESS_STOP;
        };

      for (auto i = _TCD_labels_.begin(); i != _TCD_labels_.end(); ++i) {
        if (! data_record_.has(*i)) {
          DT_LOG_WARNING(get_logging_priority(), "No data bank with label '" << *i << "' found !");
          continue;
        }
        for (auto j = std::next(i); j != _TCD_labels_.end(); ++j) {
          if (! data_record_.has(*j)) {
            DT_LOG_WARNING(get_logging_priority(), "No data bank with label '" << *j << "' found !");
            continue;
          }

          DT_LOG_DEBUG(get_logging_priority(), "Comparing '" << *i << "' bank with '" << *j << "' bank");

          const auto & tcd1 = data_record_.get<snemo::datamodel::tracker_clustering_data>(*i);
          const auto & tcd2 = data_record_.get<snemo::datamodel::tracker_clustering_data>(*j);
          if (tcd1.get_number_of_solutions() != tcd2.get_number_of_solutions()) {
            exit_message(*i, *j);
          }

          for (size_t isol = 0; isol < tcd1.get_number_of_solutions(); isol++) {
            const auto & sol1 = tcd1.get_solution(isol);
            const auto & sol2 = tcd2.get_solution(isol);
            const auto & clusters1 = sol1.get_clusters();
            const auto & clusters2 = sol2.get_clusters();
            if (clusters1.size() != clusters2.size()) {
              exit_message(*i, *j);
            }
            // for (const auto & cluster1 : clusters1) {

            // }
          }




          DT_LOG_DEBUG(get_logging_priority(), "Banks '" << *i << "' and '" << *j << "' have similar content");
        }
      }
      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
      return dpp::base_module::PROCESS_SUCCESS;
    }

  } // end of namespace reconstruction

} // end of namespace snemo
