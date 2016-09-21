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

    void prepare_event(topology::tracked_data & tracked_data_);
    void clusterize(topology::tracked_data & tracked_data_);
    void order_cells();

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

    int nevent;

    //geom param
    double calo_X, calo_Y, calo_Z;
    double pmax;
    double bfield;

    //limits
    double SmallRadius;
    double TangentPhi;
    double TangentTheta;
    double SmallNumber;
    double QuadrantAngle;
    double Ratio;
    double CompatibilityDistance;
    double MaxChi2;
    double probmin;
    size_t nofflayers;

    // Support numbers
    bool SuperNemo;
    bool SuperNemoChannel; /** New initialization modeof the algorithm
                            *  for SuperNEMO and usage from Channel by
                            *  Falaise and Hereward.
                            *  Use the GG_CELL_pitch as the main geoemtry parameter
                            *  of a GG cell, do not use 'rad' or 'CellDistance'
                            */
    bool NemoraOutput;
    bool N3_MC;

    int num_blocks;
    mybhep::dvector<double> planes_per_block ;
    mybhep::dvector<double> gaps_Z;
    int num_cells_per_plane;

    //  size_t dp_mode;

    //----Modification for bar-module---
  private:
    std::string  _moduleNR;
    int     _MaxBlockSize;

    bool _is_good_couplet_(const topology::cell & main_cell_,
                           const topology::cell & candidate_cell_,
                           const std::vector<topology::cell> & cells_near_main_);

  protected:

    /// Set the initialization flag
    void _set_initialized(bool);

    /// Set default attribute values
    void _set_defaults ();

  public:

    void set_lastlayer(int ll_);

    void set_num_blocks(int nb);

    void set_planes_per_block(int block, int nplanes);

    void set_num_cells_per_plane(int ncpp);

    void set_pmax(double v);

    void set_SmallRadius(double v);

    void set_TangentPhi(double v);

    void set_TangentTheta(double v);

    void set_SmallNumber(double v);

    void set_QuadrantAngle(double v);

    void set_Ratio(double v);

    void set_CompatibilityDistance(double v);

    void set_MaxChi2(double v);

    void set_probmin(double v);

    void set_nofflayers(size_t v);

    void set_SuperNemo(bool v);

    void set_SuperNemoChannel(bool v);

    void set_bfield(double v);

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
