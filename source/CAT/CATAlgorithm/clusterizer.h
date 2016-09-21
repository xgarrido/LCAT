/* -*- mode: c++ -*- */
#ifndef CAT_CLUSTERIZER_H
#define CAT_CLUSTERIZER_H 1

#include <stdexcept>
#include <iostream>
#include <vector>

#include <CATAlgorithm/line.h>
#include <CATAlgorithm/cell_couplet.h>
#include <CATAlgorithm/cell_triplet.h>
#include <CATAlgorithm/experimental_double.h>

#include <CATAlgorithm/cell_base.h>
#include <CATAlgorithm/cluster.h>
#include <CATAlgorithm/calorimeter_hit.h>
#include <CATAlgorithm/sequence_base.h>
#include <CATAlgorithm/tracked_data_base.h>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/logger.h>

namespace CAT{

  /// The clusterizer algorithm
  class clusterizer{

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
    void initialize();

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

    int cell_side( const topology::cell & c);
    size_t near_level( const topology::cell & c1, const topology::cell & c2 );
    std::vector<topology::cell> get_near_cells(const topology::cell & c);

  protected:

    //limits
    double TangentPhi;
    double TangentTheta;
    double QuadrantAngle;
    double Ratio;

  private:

    bool _is_good_couplet_(const topology::cell & main_cell_,
                           const topology::cell & candidate_cell_,
                           const std::vector<topology::cell> & cells_near_main_);

  protected:

    /// Set the initialization flag
    void _set_initialized(bool);

    /// Set default attribute values
    void _set_defaults();

  public:

    void set_TangentPhi(double v);

    void set_TangentTheta(double v);

    void set_QuadrantAngle(double v);

    void set_Ratio(double v);


    //----------------------------------------


  private:

    bool _initialized_;           //!< Initialization status
    datatools::logger::priority _logging_;//!< Logging priority

    std::vector<topology::cell> cells_;
    std::vector<topology::cluster> clusters_;
    std::vector<topology::calorimeter_hit> calorimeter_hits_;
    std::vector<topology::sequence> true_sequences_;

  };

}

#endif
