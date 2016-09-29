/* -*- mode: c++ -*- */
#ifndef CAT_CLUSTERIZER_H
#define CAT_CLUSTERIZER_H 1

#include <stdexcept>
#include <iostream>
#include <vector>

#include <CAT/line.h>
#include <CAT/cell_couplet.h>
#include <CAT/cell_triplet.h>
#include <CAT/experimental_double.h>

#include <CAT/cell_base.h>
#include <CAT/cluster.h>
#include <CAT/calorimeter_hit.h>
#include <CAT/sequence_base.h>
#include <CAT/tracked_data.h>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/logger.h>

// Forward declaration
namespace datatools {
  class properties;
}

namespace CAT{

  /// The clusterizer algorithm
  class clusterizer {

  public:

    /// Set logging priority
    void set_logging_priority(datatools::logger::priority p);

    /// Returns logging priority
    datatools::logger::priority get_logging_priority() const;

    /// Check if the clusterizer is initialized
    bool is_initialized() const;

    /// Constructor
    clusterizer();

    /// Destructor
    ~clusterizer();

    /// Initialize the clusterizer through configuration properties (not yet)
    void initialize(const datatools::properties & setup_);

    /// Reset the clusterizer
    void reset();

    /// Main algorithm
    void clusterize(topology::tracked_data & tracked_data_);

    //! get cells
    const std::vector<topology::cell>& get_cells()const;

    //! set cells
    void set_cells(const std::vector<topology::cell> & cells);

    //! get clusters
    const std::vector<topology::cluster>& get_clusters()const;

    //! set clusters
    void set_clusters(const std::vector<topology::cluster> & clusters);

    //! get calorimeter_hits
    const std::vector<topology::calorimeter_hit>& get_calorimeter_hits()const;

    //! set calorimeter_hits
    void set_calorimeter_hits(const std::vector<topology::calorimeter_hit> & calorimeter_hits);

  protected:

    /// Set the initialization flag
    void _set_initialized(bool);

    /// Set default attribute values
    void _set_defaults();

  private:

    /// Return if cells are good couplet given the neighboring
    bool _is_good_couplet_(const topology::cell & c1_,
                           const topology::cell & c2_,
                           const std::vector<topology::cell> & c1_neighbors_) const;

    /// Get level of 'closeness' of cells
    size_t _near_level_(const topology::cell & c1_, const topology::cell & c2_) const;

  private:

    bool _initialized_;           //!< Initialization status
    datatools::logger::priority _logging_;//!< Logging priority

    //limits
    double _tangent_phi_;
    double _ratio_;

    std::vector<topology::cell> _cells_;
    std::vector<topology::cluster> _clusters_;
    std::vector<topology::calorimeter_hit> _calorimeter_hits_;

  };

}

#endif
