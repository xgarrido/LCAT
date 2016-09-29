/// \file falaise/snemo/reconstruction/cat_driver.cc

// Ourselves:
#include <falaise/snemo/reconstruction/cat_driver.h>

// Standard library:
#include <sstream>
#include <stdexcept>

// Third party:
// - Boost :
#include <boost/foreach.hpp>
// - Bayeux/datatools :
#include <datatools/properties.h>
#include <datatools/exception.h>
#include <datatools/logger.h>
#include <datatools/utils.h>
// - Bayeux/geomtools :
#include <geomtools/manager.h>

// This project :
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/geometry/gg_locator.h>
#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/xcalo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/geometry/locator_plugin.h>

namespace snemo {

  namespace reconstruction {

    const std::string cat_driver::CAT_ID = "CAT";

    void cat_driver::set_magfield(double value_)
    {
      DT_THROW_IF(is_initialized(), std::logic_error, "CAT driver is already initialized!");
      _magfield_ = value_;
      return;
    }

    double cat_driver::get_magfield() const
    {
      return _magfield_;
    }

    cat_driver::cat_driver() :
      ::snemo::processing::base_tracker_clusterizer(cat_driver::CAT_ID)
    {
      _set_defaults();
      return;
    }

    cat_driver::~cat_driver()
    {
      if (is_initialized()) {
        this->cat_driver::reset();
      }
      return;
    }

    // Initialize the clusterizer through configuration properties
    void cat_driver::initialize(const datatools::properties & setup_)
    {
      // Invoke initialization at parent level :
      this->snemo::processing::base_tracker_clusterizer::_initialize(setup_);

      _CAT_setup_.SuperNemo = true;
      _CAT_setup_.FoilRadius = 0.0;

      double default_magfield_unit = CLHEP::tesla;
      // Forcing magnetic field (a temporary trick)
      if (! datatools::is_valid(_magfield_)) {
        if (setup_.has_key("CAT.magnetic_field")) {
          double magfield = setup_.fetch_real("CAT.magnetic_field");
          if (! setup_.has_explicit_unit("CAT.magnetic_field")) {
            magfield *= default_magfield_unit;
          }
          set_magfield(magfield);
        }
      }
      if (!datatools::is_valid(_magfield_)) {
        set_magfield(0.0025 * CLHEP::tesla);
      }

      // Verbosity level
      if (setup_.has_key("CAT.logging_priority")) {
        _CAT_setup_.level = datatools::logger::get_priority(setup_.fetch_string("CAT.level"));
        DT_THROW_IF(_CAT_setup_.level == datatools::logger::PRIO_UNDEFINED,
                    std::logic_error, "Invalid logging priority label '" << _CAT_setup_.level << "' !");
      }

      // Process calorimeter hits
      if (setup_.has_key("CAT.process_calo_hits")) {
        _process_calo_hits_ = setup_.fetch_boolean("CAT.process_calo_hits");
      }

      // Maximum processing time
      if (setup_.has_key("CAT.max_time")) {
        _CAT_setup_.MaxTime = setup_.fetch_real("CAT.max_time");
        if (! setup_.has_explicit_unit("CAT.max_time")) {
          _CAT_setup_.MaxTime *= CLHEP::ms;
        }
      }

      // Max radius of cells to be not treated as points in distance unit
      if (setup_.has_key("CAT.small_radius")) {
        _CAT_setup_.SmallRadius = setup_.fetch_real("CAT.small_radius");
        if (! setup_.has_explicit_unit("CAT.small_radius")) {
          _CAT_setup_.SmallRadius *= CLHEP::mm;
        }
      }

      // Minimal probability away from the straight line
      if (setup_.has_key("CAT.probmin")) {
        _CAT_setup_.probmin = setup_.fetch_real("CAT.probmin");
      }

      // Number of cells which can be skipped (because the cell did not
      // work) and still the cluster is continuous
      if (setup_.has_key("CAT.nofflayers")) {
        _CAT_setup_.nofflayers = setup_.fetch_integer("CAT.nofflayers");
        if (_CAT_setup_.nofflayers < 0) {
          DT_THROW_IF(true, std::logic_error,
                      "Invalid number of 'off' layers(" << _CAT_setup_.nofflayers << ") !");
        }
      }

      // Ratio of 2nd best to best chi2 which is acceptable as 2nd solution
      if (setup_.has_key("CAT.ratio")) {
        _CAT_setup_.Ratio = setup_.fetch_real("CAT.ratio");
      }

      // Sigma Z factor
      if (setup_.has_key("CAT.sigma_z_factor")) {
        _sigma_z_factor_ = setup_.fetch_real("CAT.sigma_z_factor");
        DT_THROW_IF (_sigma_z_factor_ <= 0.0 || _sigma_z_factor_ >= 100.0,
                     std::logic_error,
                     "Invalid Sigma Z factor(" << _sigma_z_factor_ << ") !");
      }

      // Get the calorimeter locators from a geometry plugin :
      const geomtools::manager & geo_mgr = get_geometry_manager();
      std::string locator_plugin_name;
      if (setup_.has_key("locator_plugin_name")) {
        locator_plugin_name = setup_.fetch_string("locator_plugin_name");
      }
      // If no locator plugin name is set, then search for the first one
      if (locator_plugin_name.empty()) {
        const geomtools::manager::plugins_dict_type & plugins = geo_mgr.get_plugins();
        for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin();
             ip != plugins.end();
             ip++) {
          const std::string & plugin_name = ip->first;
          if (geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(plugin_name)) {
            DT_LOG_DEBUG(get_logging_priority(), "Find locator plugin with name = " << plugin_name);
            locator_plugin_name = plugin_name;
            break;
          }
        }
      }
      // Access to a given plugin by name and type :
      if (geo_mgr.has_plugin(locator_plugin_name)
          && geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(locator_plugin_name)) {
        DT_LOG_NOTICE(get_logging_priority(), "Found locator plugin named '" << locator_plugin_name << "'");
        const snemo::geometry::locator_plugin & lp
          = geo_mgr.get_plugin<snemo::geometry::locator_plugin>(locator_plugin_name);
        // Set the calo cell locator :
        _calo_locator_ = &(lp.get_calo_locator());
        _xcalo_locator_ = &(lp.get_xcalo_locator());
        _gveto_locator_ = &(lp.get_gveto_locator());
      }
      if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) {
        DT_LOG_DEBUG(get_logging_priority(), "Calo locator :");
        _calo_locator_->tree_dump(std::clog, "", "[debug]: ");
        DT_LOG_DEBUG(get_logging_priority(), "X-calo locator :");
        _xcalo_locator_->tree_dump(std::clog, "", "[debug]: ");
        DT_LOG_DEBUG(get_logging_priority(), "G-veto locator :");
        _gveto_locator_->tree_dump(std::clog, "", "[debug]: ");
      }

      // Geometry description :
      _CAT_setup_.num_blocks = 1;
      _CAT_setup_.planes_per_block.clear();
      _CAT_setup_.planes_per_block.push_back(_CAT_setup_.num_blocks);
      _CAT_setup_.planes_per_block.at(0) = get_gg_locator().get_number_of_layers(0);
      _CAT_setup_.num_cells_per_plane    = get_gg_locator().get_number_of_rows(0);
      _CAT_setup_.cell_size              = get_gg_locator().get_cell_diameter();

      // Hard-coded values of bfield and chamber size
      _CAT_setup_.bfield = _magfield_ / CLHEP::tesla;
      _CAT_setup_.xsize  = 2500. * CLHEP::mm; // this is y in SnWare coordinates
      _CAT_setup_.ysize  = 1350. * CLHEP::mm; // this is z in SnWare coordinates
      _CAT_setup_.zsize  =  450. * CLHEP::mm; // this is x in SnWare coordinates

      // Check the validity of the CAT setup data :
      DT_THROW_IF(! _CAT_setup_.check(), std::logic_error,
                  "Setup data for the CAT machine is not checked : "
                  << _CAT_setup_.get_error_message() << " !");

      // Configure and initialize the CAT machine :
      // CAT::clusterizer_configure(_CAT_clusterizer_, _CAT_setup_);
      CAT::sequentiator_configure(_CAT_sequentiator_, _CAT_setup_);

      // Extract the setup of the base gamma builder :
      datatools::properties cat_setup;
      setup_.export_and_rename_starting_with(cat_setup, "CAT.", "");
      _CAT_clusterizer_.initialize(cat_setup);
      _CAT_sequentiator_.initialize();

      _set_initialized(true);
      return;
    }

    void cat_driver::_set_defaults()
    {
      _CAT_setup_.reset();
      _sigma_z_factor_ = 1.0;
      datatools::invalidate(_magfield_);
      _process_calo_hits_ = true;
      _calo_locator_  = 0;
      _xcalo_locator_ = 0;
      _gveto_locator_ = 0;
      this->base_tracker_clusterizer::_reset();
      return;
    }

    // Reset the clusterizer
    void cat_driver::reset()
    {
      DT_THROW_IF(! is_initialized(), std::logic_error,
                  "CAT driver is not initialized !");
      _set_initialized(false);
      _CAT_clusterizer_.reset();
      _CAT_sequentiator_.reset();
      _CAT_setup_.reset();
      _set_defaults();
      this->base_tracker_clusterizer::_reset();
      return;
    }


    /// Main clustering method
    int cat_driver::_process_algo(const base_tracker_clusterizer::hit_collection_type & gg_hits_,
                                  const base_tracker_clusterizer::calo_hit_collection_type & calo_hits_,
                                  snemo::datamodel::tracker_clustering_data & clustering_ )
    {
      namespace ct = CAT::topology;
      namespace sdm = snemo::datamodel;

      // CAT input data model :
      _CAT_input_.cells.clear();
      if (_CAT_input_.cells.capacity() < gg_hits_.size()){
        _CAT_input_.cells.reserve(gg_hits_.size());
      }

      _CAT_output_.tracked_data.reset();

      // Hit accounting :
      std::map<int, sdm::calibrated_data::tracker_hit_handle_type> hits_mapping;

      size_t ihit = 0;
      // GG hit loop :
      for (const auto & gg_handle : gg_hits_) {
        // Skip NULL handle :
        if (! gg_handle) continue;

        // Get a const reference on the calibrated Geiger hit :
        const sdm::calibrated_tracker_hit & snemo_gg_hit = gg_handle.get();

        // Check the geometry ID as a Geiger cell :
        const snemo::geometry::gg_locator & gg_locator = get_gg_locator();
        const geomtools::geom_id & gg_hit_gid = snemo_gg_hit.get_geom_id();
        DT_THROW_IF(! gg_locator.is_drift_cell_volume(gg_hit_gid),
                    std::logic_error,
                    "Calibrated tracker hit can not be located inside detector !");

        if (! gg_locator.is_drift_cell_volume_in_current_module(gg_hit_gid)) {
          DT_LOG_DEBUG(get_logging_priority(), "Current Geiger cell is not in the module!");
          continue;
        }

        // Extract the numbering scheme of the cell from its geom ID :
        const int side  = gg_locator.extract_side(gg_hit_gid);
        const int layer = gg_locator.extract_layer(gg_hit_gid);
        const int row   = gg_locator.extract_row(gg_hit_gid);

        // Translate into the CAT's numbering scheme :
        // -1 : negative X side; +1 : positive X side
        const int block_id = (side == 0) ? -1 : +1;
        /* Later number:
         * negative X side: layer from -8, -7, -6, -5, -4, -3, -2, -1, 0
         * positive X side: layer from 0, 1, 2 , 3, 4, 5, 6, 7, 8
         */
        const int layer_id = side == 0 ? -layer: layer;
        // Cell number :  0:113 is mapped to -56:+56
        const int cell_id  = row - (_CAT_setup_.num_cells_per_plane / 2);

        // X-Y position of the anode wire of the hit cell :
        ct::experimental_double z; // == X in sngeometry SN module frame
        ct::experimental_double x; // == Y in sngeometry SN module frame

        // Center of the cell set in CAT's own reference frame:
        z.set_value(snemo_gg_hit.get_x());
        z.set_error(0.0);
        x.set_value(snemo_gg_hit.get_y());
        x.set_error(0.0);

        // Transverse Geiger drift distance :
        ct::experimental_double y;
        // Plasma longitudinal origin along the anode wire :
        y.set_value(snemo_gg_hit.get_z());
        y.set_error(_sigma_z_factor_ * snemo_gg_hit.get_sigma_z());

        // Prompt/delayed trait of the hit :
        const bool fast = snemo_gg_hit.is_prompt();

        // Transverse Geiger drift distance :
        const double rdrift     = snemo_gg_hit.get_r();
        const double rdrift_err = snemo_gg_hit.get_sigma_r();

        // Build the Geiger hit position :
        ct::experimental_point gg_hit_position(x,y,z);

        // Add a new hit cell in the CAT input data model :
        ct::cell & c = _CAT_output_.tracked_data.add_gg_cell();
        c.set_id(ihit++);
        c.set_probmin(_CAT_setup_.probmin);
        c.set_p(gg_hit_position);
        c.set_r(fast ? rdrift : 0.25 * gg_locator.get_cell_diameter());
        c.set_er(fast ? rdrift_err : 0.25 * gg_locator.get_cell_diameter());
        c.set_layer(layer_id);
        c.set_block(block_id);
        c.set_iid(cell_id);
        c.set_fast(fast);
        c.set_small_radius(_CAT_setup_.SmallRadius);

        // Store mapping info between both data models :
        hits_mapping[c.id()] = gg_handle;

        DT_LOG_DEBUG (get_logging_priority (),
                      "Geiger cell #" << snemo_gg_hit.get_id() << " has been added "
                      << "to CAT input data with id number #" << c.id());
      } // for loop over gg_hits_

      // Take into account calo hits:
      // Calo hit accounting :
      std::map<int, sdm::calibrated_data::calorimeter_hit_handle_type> calo_hits_mapping;
      size_t jhit = 0;

      // CALO hit loop :
      for (const auto & calo_handle : calo_hits_) {
        // Skip NULL handle :
        if (! calo_handle.has_data()) continue;

        // Get a const reference on the calibrated Calo hit :
        const sdm::calibrated_calorimeter_hit & sncore_calo_hit = calo_handle.get();

        // Get calibrated calo. geom_id
        const geomtools::geom_id & a_calo_hit_gid = sncore_calo_hit.get_geom_id();
        // Extract the numbering scheme of the calo_cell from its geom ID :
        int column = -1;
        int side = -1;
        double width = datatools::invalid_real();
        double height = datatools::invalid_real();
        double thickness = datatools::invalid_real();
        ct::experimental_vector norm(0., 0., 0., 0., 0., 0.);
        geomtools::vector_3d block_position;
        // Extract the numbering scheme of the scin block from its geom ID :
        if (_calo_locator_->is_calo_block_in_current_module(a_calo_hit_gid)) {
          _calo_locator_->get_block_position(a_calo_hit_gid, block_position);
          width     = _calo_locator_->get_block_width();
          height    = _calo_locator_->get_block_height();
          thickness = _calo_locator_->get_block_thickness();
          column    = _calo_locator_->extract_column(a_calo_hit_gid);
          side      = _calo_locator_->extract_side(a_calo_hit_gid);
          const int side_number = (side == snemo::geometry::utils::SIDE_BACK) ? 1: -1;
          norm.set_x(ct::experimental_double((double) side_number, 0.));
        } else if (_xcalo_locator_->is_calo_block_in_current_module(a_calo_hit_gid)) {
          _xcalo_locator_->get_block_position(a_calo_hit_gid, block_position);
          width     = _xcalo_locator_->get_block_width();
          height    = _xcalo_locator_->get_block_height();
          thickness = _xcalo_locator_->get_block_thickness();
          column    = _xcalo_locator_->extract_column(a_calo_hit_gid);
          side      = _xcalo_locator_->extract_side(a_calo_hit_gid);
          const int side_number = (side == snemo::geometry::utils::SIDE_BACK) ? 1: -1;
          norm.set_y(ct::experimental_double((double) side_number, 0.));
        } else if (_gveto_locator_->is_calo_block_in_current_module(a_calo_hit_gid)) {
          _gveto_locator_->get_block_position(a_calo_hit_gid, block_position);
          width     = _gveto_locator_->get_block_width();
          height    = _gveto_locator_->get_block_height();
          thickness = _gveto_locator_->get_block_thickness();
          column    = _gveto_locator_->extract_column(a_calo_hit_gid);
          side      = _xcalo_locator_->extract_side(a_calo_hit_gid);
          const int side_number = (side == snemo::geometry::utils::SIDE_BACK) ? 1: -1;
          norm.set_z(ct::experimental_double((double) side_number, 0.));
        }

        ct::experimental_double energy(sncore_calo_hit.get_energy(),
                                       sncore_calo_hit.get_sigma_energy());
        ct::experimental_double time(sncore_calo_hit.get_time(),
                                     sncore_calo_hit.get_sigma_time());
        // size_t id = sncore_calo_hit.get_hit_id();
        ct::experimental_point center(block_position.x(),
                                      block_position.y(),
                                      block_position.z(),
                                      0., 0., 0.);
        ct::experimental_vector sizes(width, height, thickness,
                                      0., 0., 0.);
        ct::plane pl(center, sizes, norm);
        pl.set_probmin(_CAT_setup_.probmin);
        pl.set_type("SuperNEMO");

        // Build the Calo hit position :
        // Add a new hit calo_cell in the CAT input data model :
        ct::calorimeter_hit & c = _CAT_output_.tracked_data.add_calo_cell();
        c.set_pl(pl);
        c.set_e(energy);
        c.set_t(time);
        c.set_probmin(_CAT_setup_.probmin);
        c.set_layer(column);
        c.set_id(jhit++);

        // Store mapping info between both data models :
        calo_hits_mapping[c.id()] = calo_handle;

        DT_LOG_DEBUG(get_logging_priority(),
                     "Calo_cell #" << sncore_calo_hit.get_hit_id() << " has been added "
                     << "to CAT input data with id number #" << c.id());
      }

      // // Validate the input data :
      // if (! _CAT_input_.check()) {
      //   DT_LOG_ERROR(get_logging_priority (), "Invalid CAT input data !");
      //   return 1;
      // }

      // // Install the input data model within the algorithm object :
      // _CAT_clusterizer_.set_cells(_CAT_input_.cells);

      // // Install the input data model within the algorithm object :
      // _CAT_clusterizer_.set_calorimeter_hits(_CAT_input_.calo_cells);

      // Run the clusterizer algorithm :
      _CAT_clusterizer_.clusterize(_CAT_output_.tracked_data);

      const bool use_only_clusterizer = false;
      if (use_only_clusterizer) {
        // Store solution from clusterizer
        sdm::tracker_clustering_solution::handle_type htcs(new sdm::tracker_clustering_solution);
        clustering_.add_solution(htcs, true);
        clustering_.grab_default_solution().set_solution_id(clustering_.get_number_of_solutions() - 1);
        sdm::tracker_clustering_solution & clustering_solution = clustering_.grab_default_solution();
        clustering_solution.grab_auxiliaries().update_string(sdm::tracker_clustering_data::clusterizer_id_key(), "LCAT");

        for (const auto & a_cluster : _CAT_output_.tracked_data.get_clusters()) {
          {
            // Append a new cluster :
            sdm::tracker_cluster::handle_type tch(new sdm::tracker_cluster);
            clustering_solution.grab_clusters().push_back(tch);
          }
          sdm::tracker_cluster::handle_type & cluster_handle = clustering_solution.grab_clusters().back();
          cluster_handle.grab().set_cluster_id(clustering_solution.get_clusters().size() - 1);

          // Loop on all hits within the cluster(nodes) :
          for (const auto & a_node : a_cluster.nodes()) {
            const int hit_id = a_node.c().id();
            cluster_handle.grab().grab_hits().push_back(hits_mapping[hit_id]);
            DT_LOG_DEBUG(get_logging_priority(), "Add tracker hit with id #" << hit_id);
          }
        }
      } else {
        // Run the sequentiator algorithm :
        _CAT_sequentiator_.sequentiate(_CAT_output_.tracked_data);

        // _CAT_output_.tracked_data.dump();

        // Analyse the sequentiator output i.e. 'scenarios' made of 'sequences' of geiger cells:
        const std::vector<CAT::topology::scenario> & tss = _CAT_output_.tracked_data.get_scenarios();

        for (std::vector<CAT::topology::scenario>::const_iterator iscenario = tss.begin();
             iscenario != tss.end();
             ++iscenario) {
          DT_LOG_DEBUG(get_logging_priority(), "Number of scenarios = " << tss.size());

          sdm::tracker_clustering_solution::handle_type htcs(new sdm::tracker_clustering_solution);
          clustering_.add_solution(htcs, true);
          clustering_.grab_default_solution().set_solution_id(clustering_.get_number_of_solutions() - 1);
          sdm::tracker_clustering_solution & clustering_solution = clustering_.grab_default_solution();
          clustering_solution.grab_auxiliaries().update_string(sdm::tracker_clustering_data::clusterizer_id_key(), CAT_ID);

          // Analyse the sequentiator output :
          const std::vector<CAT::topology::sequence> & the_sequences = iscenario->sequences();
          DT_LOG_DEBUG(get_logging_priority(), "Number of sequences = " << the_sequences.size());

          for (std::vector<CAT::topology::sequence>::const_iterator isequence = the_sequences.begin();
               isequence != the_sequences.end();
               ++isequence) {
            const CAT::topology::sequence & a_sequence = *isequence;
            const size_t seqsz = a_sequence.nodes().size();
            if (seqsz > 1) {
              // A CAT cluster with more than one hit/cell(node) :
              {
                // Append a new cluster :
                sdm::tracker_cluster::handle_type tch(new sdm::tracker_cluster);
                clustering_solution.grab_clusters().push_back(tch);
              }
              sdm::tracker_cluster::handle_type & cluster_handle
                = clustering_solution.grab_clusters().back();
              cluster_handle.grab().set_cluster_id(clustering_solution.get_clusters().size() - 1);

              // Loop on all hits within the sequence(nodes) :
              for (size_t i = 0; i < seqsz; i++) {
                const CAT::topology::node & a_node = a_sequence.nodes()[i];
                const int hit_id = a_node.c().id();
                cluster_handle.grab().grab_hits().push_back(hits_mapping[hit_id]);
                DT_LOG_DEBUG(get_logging_priority(), "Add tracker hit with id #" << hit_id);
              }
            }
          } /* for sequence */
        } // finish loop on scenario
      }

      if (get_logging_priority() >= datatools::logger::PRIO_TRACE) {
        clustering_.tree_dump(std::clog, "Output clustering data : ", "[trace]: ");
      }
      return 0;
    }

    // static
    void cat_driver::init_ocd(datatools::object_configuration_description & ocd_)
    {

      // Invoke OCD support from parent class :
      ::snemo::processing::base_tracker_clusterizer::ocd_support(ocd_);

      {
        // Description of the 'CAT.magnetic_field' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.magnetic_field")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Force the magnetic field value (vertical)")
          .set_traits(datatools::TYPE_REAL)
          .set_mandatory(false)
          // .set_long_description("Default value: 25 gauss")
          .set_default_value_real(25 * CLHEP::gauss, "gauss")
          .add_example("Use no magnetic field::               \n"
                       "                                      \n"
                       "  CAT.magnetic_field : real = 0 gauss \n"
                       "                                      \n"
                       )
          ;
      }


      {
        // Description of the 'CAT.level' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.level")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Verbosity level")
          .set_traits(datatools::TYPE_STRING)
          .set_mandatory(false)
          // .set_long_description("Default value: \"normal\"")
          .set_default_value_string("normal")
          .add_example("Use normal verbosity:: \n"
                       "                                  \n"
                       "  CAT.level : string = \"normal\" \n"
                       "                                  \n"
                       )
          ;
      }

      {
        // Description of the 'CAT.max_time' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.max_time")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Maximum processing time")
          .set_traits(datatools::TYPE_REAL)
          .set_mandatory(false)
          // .set_long_description("Default value: 5000 ms")
          .set_default_value_real(5000 * CLHEP::ms, "ms")
          .add_example("Use default value::               \n"
                       "                                  \n"
                       "  CAT.max_time : real = 5000 ms   \n"
                       "                                  \n"
                       )
          ;
      }

      {
        // Description of the 'CAT.small_radius' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.small_radius")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Max radius of cells to be not treated as points in distance unit")
          .set_traits(datatools::TYPE_REAL)
          .set_mandatory(false)
          .set_long_description("Default value: 2.0 mm")
          .add_example("Use default value::                \n"
                       "                                   \n"
                       "  CAT.small_radius : real = 2.0 mm \n"
                       "                                   \n"
                       )
          ;
      }

      {
        // Description of the 'CAT.probmin' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.probmin")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Minimal probability away from the straight line")
          .set_traits(datatools::TYPE_REAL)
          .set_mandatory(false)
          .set_long_description("Default value: 0.0")
          .add_example("Use default value::               \n"
                       "                                  \n"
                       "  CAT.probmin : real = 0.0        \n"
                       "                                  \n"
                       )
          ;
      }

      {
        // Description of the 'CAT.nofflayers' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.nofflayers")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Number of cells which can be skipped (because the cell did not work) and still the cluster is continuous")
          .set_traits(datatools::TYPE_INTEGER)
          .set_mandatory(false)
          .set_long_description("Default value: 1")
          .add_example("Use default value::               \n"
                       "                                  \n"
                       "  CAT.nofflayers : integer = 1    \n"
                       "                                  \n"
                       )
          ;
      }

      {
        // Description of the 'CAT.first_event' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.first_event")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("First event to be processed")
          .set_traits(datatools::TYPE_INTEGER)
          .set_mandatory(false)
          .set_long_description("Default value: -1")
          .add_example("Do not specify any first event::  \n"
                       "                                  \n"
                       "  CAT.first_event : integer = -1  \n"
                       "                                  \n"
                       )
          ;
      }


      {
        // Description of the 'CAT.ratio' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.ratio")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Ratio of 2nd best to best chi2 which is acceptable as 2nd solution")
          .set_traits(datatools::TYPE_REAL)
          .set_mandatory(false)
          .set_long_description("Default value: 10000.0")
          .add_example("Use the default value::          \n"
                       "                                 \n"
                       "  CAT.ratio : real = 10000.0     \n"
                       "                                 \n"
                       )
          ;
      }


      {
        // Description of the 'CAT.driver.sigma_z_factor' configuration property :
        datatools::configuration_property_description & cpd
          = ocd_.add_property_info();
        cpd.set_name_pattern("CAT.sigma_z_factor")
          .set_from("snemo::reconstruction::cat_driver")
          .set_terse_description("Sigma Z factor")
          .set_traits(datatools::TYPE_REAL)
          .set_mandatory(false)
          .set_long_description("Default value: 1.0")
          .add_example("Use the default value::                  \n"
                       "                                         \n"
                       "  CAT.sigma_z_factor : real = 1.0        \n"
                       "                                         \n"
                       )
          ;
      }

      return;
    }

  }  // end of namespace reconstruction

}  // end of namespace snemo

/* OCD support */
#include <datatools/object_configuration_description.h>
DOCD_CLASS_IMPLEMENT_LOAD_BEGIN(snemo::reconstruction::cat_driver,ocd_)
{
  ocd_.set_class_name("snemo::reconstruction::cat_driver");
  ocd_.set_class_description("A driver class for the CAT algorithm");
  ocd_.set_class_library("Falaise_CAT");
  ocd_.set_class_documentation("This driver manager for the CAT clustering algorithm.");

  // Invoke specific OCD support :
  ::snemo::reconstruction::cat_driver::init_ocd(ocd_);

  ocd_.set_validation_support(true);
  ocd_.lock();
  return;
}
DOCD_CLASS_IMPLEMENT_LOAD_END() // Closing macro for implementation
// DOCD_CLASS_SYSTEM_REGISTRATION(snemo::reconstruction::cat_driver, "snemo::reconstruction::cat_driver")
