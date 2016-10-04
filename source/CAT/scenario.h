// -*- mode: c++ -*-

#ifndef FALAISE_CAT_SCENARIO_H
#define FALAISE_CAT_SCENARIO_H

// Standard library:
#include <iostream>
#include <cmath>

// Third party:
// - Boost:
#include <boost/cstdint.hpp>

// This project:
#include <CAT/experimental_point.h>
#include <CAT/experimental_vector.h>
#include <CAT/sequence_base.h>
#include <CAT/node.h>
#include <CAT/calorimeter_hit.h>

namespace CAT {

  /// \broef A reconstruction scenario is composed of
  /// a collection of tracks.
  class scenario
  {

  private:

    // chi2
    double helix_chi2_;
    double tangent_chi2_;
    int32_t ndof_;

    // n of free families
    size_t n_free_families_;

    // n of overlapping cells
    size_t n_overlaps_;

  public:

    // tracks
    std::vector<sequence> sequences_;

  public:

    //!Default constructor
    scenario();

    //!Default destructor
    virtual ~scenario();

    //! constructor
    scenario(const std::vector<sequence> & seqs);

    //! Smart print
    virtual void dump (std::ostream & a_out         = std::clog,
                       const std::string & a_title  = "",
                       const std::string & a_indent = "",
                       bool a_inherit          = false) const;


    //! set experimental_point, radius, error and id;
    void set(const std::vector<sequence> & seqs);

    //! set sequences
    void set_sequences(const std::vector<sequence> & seqs);

    //! set helix_chi2
    void set_helix_chi2(double helix_chi2);

    //! set tangent_chi2
    void set_tangent_chi2(double tangent_chi2);

    //! set n free families
    void set_n_free_families(size_t n);

    //! set n overlaps
    void set_n_overlaps(size_t n);

    //! set ndof
    void set_ndof(int32_t n);

    //! get sequences
    const std::vector<sequence> & sequences() const;

    //!get helix_chi2
    double helix_chi2() const;

    //!get tangent_chi2
    double tangent_chi2() const;

    //!get ndof
    int32_t ndof() const;

    //!get n free families
    size_t n_free_families() const;

    //!get n overlaps
    size_t n_overlaps() const;


    void calculate_n_overlaps(const std::vector<cell> & cells,
                              const std::vector<calorimeter_hit> & calos);

    void calculate_n_free_families(const std::vector<cell> &cells,
                                   const std::vector<calorimeter_hit> & calos);

    void calculate_chi2();

    double helix_Prob() const;

    double tangent_Prob() const;

    bool better_scenario_than( const scenario & s, double limit) const;

    size_t n_of_common_vertexes(double limit) const;

    size_t n_of_ends_on_wire(void) const;

  };

} // namespace CAT

#endif // FALAISE_CAT_SCENARIO_H
