/** \file snemo/reconstruction/clustering_comparison_module.h
 *
 * Description:
 *
 *   Module for Geiger hits clustering
 *
 * History:
 *
 */

#ifndef FALAISE_CAT_PLUGIN_SNEMO_RECONSTRUCTION_CLUSTERING_COMPARISON_MODULE_H
#define FALAISE_CAT_PLUGIN_SNEMO_RECONSTRUCTION_CLUSTERING_COMPARISON_MODULE_H 1

// Third party:
// - Bayeux/dpp:
#include <bayeux/dpp/base_module.h>

namespace snemo {

  namespace reconstruction {

    /// \brief A module to compare tracker clustering algorithm
    class clustering_comparison_module : public dpp::base_module
    {

    public:

      /// Constructor
      clustering_comparison_module(datatools::logger::priority = datatools::logger::PRIO_FATAL);

      /// Destructor
      virtual ~clustering_comparison_module();

      /// Initialization
      virtual void initialize(const datatools::properties  & setup_,
                              datatools::service_manager   & service_manager_,
                              dpp::module_handle_dict_type & module_dict_);

      /// Reset
      virtual void reset();

      /// Data record processing
      virtual process_status process(datatools::things & data_);

    protected:

      /// Give default values to specific class members
      void _set_defaults();

    private:

      std::vector<std::string> _TCD_labels_; //!< A list of TCD labels

      // Macro to automate the registration of the module :
      DPP_MODULE_REGISTRATION_INTERFACE(clustering_comparison_module)

    };

  } // end of namespace reconstruction

} // end of namespace snemo

#endif // FALAISE_CAT_PLUGIN_SNEMO_RECONSTRUCTION_CLUSTERING_COMPARISON_MODULE_H
/*
** Local Variables: --
** mode: c++ --
** c-file-style: "gnu" --
** tab-width: 2 --
** End: --
*/
