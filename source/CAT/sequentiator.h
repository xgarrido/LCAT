// -*- mode: c++ -*-

#ifndef CAT_SEQUENTIATOR_H
#define CAT_SEQUENTIATOR_H

// Standard library:
#include <iostream>
#include <vector>

// Third party:
// - Boost:
#include <boost/cstdint.hpp>
// bayeux/datatools:
#include <bayeux/datatools/logger.h>

// This project:
#include <CAT/cell.h>
#include <CAT/line.h>
#include <CAT/cell_couplet.h>
#include <CAT/cell_triplet.h>
#include <CAT/cluster.h>
#include <CAT/calorimeter_hit.h>
#include <CAT/sequence_base.h>
#include <CAT/experimental_double.h>
#include <CAT/Clock.h>
#include <CAT/plane.h>
#include <CAT/tracked_data.h>
#include <CAT/helix.h>
#include <CAT/scenario.h>

namespace CAT {

  /// \brief The sequentiator algorithm
  class sequentiator
  {

  public:

    /// Set logging priority
    void set_logging_priority(datatools::logger::priority p);

    /// Returns logging priority
    datatools::logger::priority get_logging_priority() const;

    /// Check if the clusterizer is initialized
    bool is_initialized() const;

    /// Constructor
    sequentiator();

    /// Destructor
    ~sequentiator();

    /// Initialize the sequentiator through configuration properties (not yet)
    void initialize();

    /// Reset the sequentiator
    void reset();

    void sequentiate(tracked_data & tracked_data);
    void sequentiate_cluster(cluster & cluster);
    void make_new_sequence(node & first_node);
    void make_copy_sequence(node & first_node);
    bool evolve(sequence & sequence);
    bool good_first_node(node & node_);
    bool good_first_to_be_matched(sequence& seq);
    bool match_gaps(const std::vector<calorimeter_hit> & calos);

    //! get clusters
    const std::vector<cluster> & get_clusters() const
    {
      return clusters_;
    }

    //! set clusters
    void set_clusters(std::vector<cluster> clusters)
    {
      clusters_.clear();
      clusters_ = clusters;
    }

    //! get sequences
    const std::vector<sequence>& get_sequences() const
    {
      return sequences_;
    }

    //! set sequences
    void set_sequences(std::vector<sequence> sequences)
    {
      sequences_.clear();
      sequences_ = sequences;
    }

    bool late();

  protected:

    /// Set the initialization flag
    void _set_initialized(bool);

    /// Set default attribute values
    void _set_defaults();

  public:

    void set_num_blocks(int nb)
    {
      if (nb > 0) {
        num_blocks = nb;
        planes_per_block.assign (num_blocks, 1);
      } else {
        std::cerr << "WARNING: CAT::clusterizer::set_num_blocks: "
                  << "Invalid number of GG layer blocks !" << std::endl;
        planes_per_block.clear ();
        num_blocks = -1; // invalid value
      }
      return;
    }

    void set_planes_per_block(int block, int nplanes)
    {
      if (block< 0 || block>=(int)planes_per_block.size()) {
        throw std::range_error ("CAT::clusterizer::set_planes_per_block: Invalid GG layer block index !");
      }
      if (nplanes > 0) {
        planes_per_block.at (block) = nplanes;
      } else {
        throw std::range_error ("CAT::clusterizer::set_planes_per_block: Invalid number of GG layers in block !");
      }
      return;
    }

    void set_num_cells_per_plane(int ncpp)
    {
      if (ncpp <= 0)
        {
          num_cells_per_plane = -1; // invalid value
        }
      else
        {
          num_cells_per_plane = ncpp;
        }
      return;
    }

    void set_SOURCE_thick(double st)
    {
      if (st <= 0.0)
        {
          SOURCE_thick = std::numeric_limits<double>::quiet_NaN ();
        }
      else
        {
          SOURCE_thick = st;
        }
      return;
    }

    void set_pmax(double v){
      if ( v <= 0.0)
        {
          pmax = std::numeric_limits<double>::quiet_NaN ();
        }
      else
        {
          pmax = v;
        }
      return;
    }

    void set_MaxTime(double v)
    {
      MaxTime = v;
      return;
    }

    void set_SmallRadius(double v)
    {
      SmallRadius = v;
      return;
    }

    void set_TangentPhi(double v)
    {
      TangentPhi = v;
      return;
    }

    void set_TangentTheta(double v)
    {
      TangentTheta = v;
      return;
    }

    void set_SmallNumber(double v)
    {
      SmallNumber = v;
      return;
    }

    void set_QuadrantAngle(double v)
    {
      QuadrantAngle = v;
      return;
    }

    void set_Ratio(double v)
    {
      Ratio = v;
      return;
    }

    double get_Ratio()
    {
      return Ratio;
    }

    void set_CompatibilityDistance(double v)
    {
      CompatibilityDistance = v;
      return;
    }

    void set_MaxChi2(double v)
    {
      MaxChi2 = v;
      return;
    }

    void set_probmin(double v)
    {
      probmin = v;
      return;
    }

    void set_nofflayers(size_t v)
    {
      NOffLayers = v;
      return;
    }

    void set_len(double v)
    {
      len = v;
      return;
    }

    void set_vel(double v)
    {
      vel = v;
      return;
    }

    void set_rad(double v)
    {
      rad = v;
      return;
    }

    void set_CellDistance(double v)
    {
      CellDistance = v;
      return;
    }

    void set_SuperNemo(bool v)
    {
      SuperNemo = v;
      return;
    }
    void set_SuperNemoChannel(bool v)
    {
      if (v)
        {
          set_SuperNemo (true);
          SuperNemoChannel = true;
          set_NemoraOutput (false);
          set_N3_MC (false);
        }
      else
        {
          SuperNemoChannel = false;
        }
      return;
    }

    void set_NemoraOutput(bool no)
    {
      NemoraOutput = no;
      return;
    }

    void set_N3_MC(bool v)
    {
      N3_MC = v;
      return;
    }

    void set_FoilRadius(double v)
    {
      FoilRadius = v;
      return;
    }

    void set_xsize(double v)
    {
      xsize = v;
      return;
    }

    void set_ysize(double v)
    {
      ysize = v;
      return;
    }

    void set_zsize(double v)
    {
      zsize = v;
      return;
    }

    void set_bfield(double v)
    {
      bfield = v;
      return;
    }

    void set_gaps_Z(std::vector<double> v)
    {
      gaps_Z.clear();
      for(size_t i=0; i<v.size(); i++)
        gaps_Z.push_back(v[i]);
    }


  protected:

    bool _initialize( void );


  protected:

    Clock clock;

    // Geometry parameters:
    double vel, rad, len, CellDistance;
    double xsize,ysize,zsize; //only for plotting
    double InnerRadius, OuterRadius, FoilRadius; //only for plotting
    double pmax;

    // Limits:
    double SmallRadius;
    double TangentPhi;
    double TangentTheta;
    double SmallNumber;
    double QuadrantAngle;
    double Ratio;
    double CompatibilityDistance;
    double MaxChi2;
    double probmin;
    int NOffLayers;

    // Error parametrization
    double sigma0;
    double k0;
    double k1;
    double k2;
    double k3;

    double th0;
    double th1;
    double th2;
    double th3;

    double pnob;
    double pnot;
    double pnobt;

    double l0;
    double l1;


    // Support numbers:
    double execution_time;
    bool SuperNemo;
    bool NemoraOutput;
    bool N3_MC;
    double MaxTime;
    bool SuperNemoChannel; /** New initialization mode of the algorithm
                            *  for SuperNEMO and usage from Channel by
                            *  Falaise and Hereward.
                            *  Use the GG_CELL_pitch as the main geoemtry parameter
                            *  of a GG cell, do not use 'rad' or 'CellDistance'
                            */

    int num_blocks;
    std::vector<double> planes_per_block ;
    std::vector<double> gaps_Z;
    int num_cells_per_plane;
    double SOURCE_thick;
    double bfield;
    int cell_max_number;

  private:

    bool make_scenarios(tracked_data &td);
    void interpret_physics(const std::vector<calorimeter_hit> & calos);
    void refine_sequences_near_walls(const std::vector<calorimeter_hit> & calos);
    bool belongs_to_other_family(cell c, sequence *iseq);
    plane get_foil_plane();
    void add_pair(const sequence & sequence);
    bool clean_up_sequences();
    bool there_is_free_sequence_beginning_with(const cell &c, size_t *index);
    int gap_number(const cell &c);
    void make_name(sequence & seq);
    bool near(const cell &c, const calorimeter_hit &ch);
    double distance_from_foil(const experimental_point &ep);
    bool direct_out_of_foil(void);
    bool direct_scenarios_out_of_foil(void);
    void make_families();
    bool can_add_family(scenario &sc, size_t* jmin, size_t* nfree, double* Chi2, size_t* noverlaps, int32_t* ndof, tracked_data &td);
    size_t pick_best_scenario();
    bool can_match(sequence &s, size_t* jmin, bool& bestinvertA, bool& bestinvertB, int& with_kink, int &cells_to_delete, const std::vector<calorimeter_hit> & calos);
    bool sequence_is_within_range(node nodeA, node nodeB, sequence seq);
    size_t near_level( const cell & c1, const cell & c2 );
    void reassign_cells_based_on_helix( sequence * seq );

  private:

    datatools::logger::priority _logging_; //!< Logging priority
    bool _initialized_;                    //!< Initialization status


    int NFAMILY, NCOPY;

    cluster * local_cluster_;

    std::vector<cluster> clusters_;
    std::vector<sequence> sequences_;

    // tables from switching from true to reco sequences
    std::vector<size_t> reco_sequence_of_true_;
    std::vector<size_t> true_sequence_of_reco_;
    std::vector<size_t> n_common_hits_for_reco_track_;
    std::vector<std::vector<size_t> > families_;
    std::vector<scenario> scenarios_;


  };

} // end of namespace CAT

#endif // CAT_SEQUENTIATOR_H
