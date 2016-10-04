/* -*- mode: c++ -*- */
#ifndef CAT_CLUSTERIZER_H
#define CAT_CLUSTERIZER_H 1

// Standard library
#include <vector>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/logger.h>

// Forward declaration
namespace datatools {
  class properties;
}

namespace CAT {

  // Forward declaration
  class tracked_data;
  class cell;

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
    void clusterize(tracked_data & tracked_data_);

  protected:

    /// Set the initialization flag
    void _set_initialized(bool);

    /// Set default attribute values
    void _set_defaults();

  private:

    /// Return if cells are good couplet given the neighboring
    bool _is_good_couplet_(const cell & c1_,
                           const cell & c2_,
                           const std::vector<cell> & c1_neighbors_) const;

    /// Get level of 'closeness' of cells
    size_t _near_level_(const cell & c1_, const cell & c2_) const;

  private:

    bool _initialized_;           //!< Initialization status
    datatools::logger::priority _logging_;//!< Logging priority

    //limits
    double _tangent_phi_;
    double _ratio_;

  };

}

#endif
