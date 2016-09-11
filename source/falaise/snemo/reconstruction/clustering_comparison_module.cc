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

      if (setup_.has_key("TCD_labels")) {
        std::vector<std::string> TCD_labels;
        setup_.fetch("TCD_labels", TCD_labels);
        // Build list of pairs of data bank label to be compared
        for (auto i = TCD_labels.begin(); i != TCD_labels.end(); ++i) {
          for (auto j = std::next(i); j != TCD_labels.end(); ++j) {
            _TCD_labels_.push_back(std::make_pair(*i, *j));
          }
        }
      }
      _set_initialized(true);
      return;
    }

    void clustering_comparison_module::reset()
    {
      DT_THROW_IF(! is_initialized(), std::logic_error,
                  "Module '" << get_name() << "' is not initialized !");

      DT_LOG_WARNING(datatools::logger::PRIO_WARNING, "Entering...");

      for (const auto & i : _nbr_event_different_) {
        const auto & a_pair = i.first;
        DT_LOG_NOTICE(get_logging_priority(), i.second << " events different between '"
                      << a_pair.first << "' and '" << a_pair.second << "' data banks.");
      }

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

      for (const auto & i : _TCD_labels_) {
        const auto & ilabel = i.first;
        const auto & jlabel = i.second;
        if (! data_record_.has(ilabel)) {
          DT_LOG_WARNING(get_logging_priority(), "No data bank with label '" << ilabel << "' found !");
          continue;
        }
        if (! data_record_.has(jlabel)) {
          DT_LOG_WARNING(get_logging_priority(), "No data bank with label '" << jlabel << "' found !");
          continue;
        }

        DT_LOG_DEBUG(get_logging_priority(), "Comparing '" << ilabel << "' bank with '" << jlabel << "' bank");

        const auto & tcd1 = data_record_.get<snemo::datamodel::tracker_clustering_data>(ilabel);
        const auto & tcd2 = data_record_.get<snemo::datamodel::tracker_clustering_data>(jlabel);

        std::ostringstream oss1, oss2;
        tcd1.tree_dump(oss1);
        tcd2.tree_dump(oss2);
        if (oss1.str() != oss2.str()) {
          DT_LOG_WARNING(datatools::logger::PRIO_ALWAYS,
                         "'" << ilabel << "' and '" << jlabel << "' banks are different !");
          DT_LOG_TRACE(get_logging_priority(), oss1.str());
          DT_LOG_TRACE(get_logging_priority(), oss2.str());
          _nbr_event_different_[std::make_pair(ilabel, jlabel)]++;
          return dpp::base_module::PROCESS_CONTINUE;
        }

        DT_LOG_DEBUG(get_logging_priority(), "Banks '" << ilabel << "' and '" << jlabel << "' have similar content");
      }
      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
      return dpp::base_module::PROCESS_SUCCESS;
    }

  } // end of namespace reconstruction

} // end of namespace snemo
