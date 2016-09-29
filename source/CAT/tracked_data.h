// -*- mode: c++ -*-

#ifndef CAT_TOPOLOGY_TRACKED_DATA_H
#define CAT_TOPOLOGY_TRACKED_DATA_H

// Standard library:
#include <iostream>
#include <vector>

// This project:
#include <CAT/cell_base.h>
#include <CAT/calorimeter_hit.h>
#include <CAT/cluster.h>
#include <CAT/scenario.h>

namespace CAT{

  namespace topology{

    /// \brief Tracked data
    ///
    /// A tracked data is composed of a list of cells
    /// a list of clusters
    /// and a list of scenarios.
    class tracked_data : public tracking_object
    {
    public:
      /// Default constructor
      tracked_data();

      /// Destructor
      virtual ~tracked_data();

      /// Reset
      void reset();

      /// Smart dump
      void tree_dump(std::ostream & out_         = std::clog,
                     const std::string & title_  = "",
                     const std::string & indent_ = "",
                     bool inherit_               = false) const;

      /// Add a Geiger hit
      topology::cell & add_gg_hit();

      /// Get a mutable reference to Geiger hits
      std::vector<cell> & grab_gg_hits();

      /// Get a non-mutable reference to Geiger hits
      const std::vector<cell> & get_gg_hits() const;

      /// Add a calorimeter cell
      topology::calorimeter_hit & add_calo_hit();

      /// Get a mutable reference to calorimeter hits
      std::vector<calorimeter_hit> & grab_calo_hits();

      /// Get a non-mutable refrence to calorimeter hits
      const std::vector<calorimeter_hit> & get_calo_hits() const;

      /// Get a mutable reference to clusters
      std::vector<cluster> & grab_clusters();

      /// Get a non-mutable reference to clusters
      const std::vector<cluster> & get_clusters() const;

      /// Get a mutable reference to scenarios
      std::vector<scenario> & grab_scenarios();

      /// Get a non-mutable reference to scenarios
      const std::vector<scenario> & get_scenarios() const;


    private:

      // list of cells
      std::vector<cell> _gg_hits_;

      // list of calos
      std::vector<calorimeter_hit> _calo_hits_;

      // list of clusters
      std::vector<cluster> _clusters_;

      // list of scenarios
      std::vector<scenario> _scenarios_;

    };

  } // namespace topology

} // namespace CAT

#endif // CAT_TOPOLOGY_TRACKED_DATA_H
