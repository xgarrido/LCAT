// -*- mode: c++ -*-

#ifndef FALAISE_CAT_JOINT_H
#define FALAISE_CAT_JOINT_H

// Standard library:
#include <iostream>
#include <cmath>

// This project:
#include <CAT/experimental_point.h>
#include <CAT/experimental_vector.h>
#include <CAT/cell.h>

namespace CAT {

  /// \brief A joint is composed of three experimental points
  class joint
  {
  private:

    // first experimental point
    experimental_point epa_;

    // second experimental point
    experimental_point epb_;

    // third experimental point
    experimental_point epc_;

    // phi kink
    experimental_double kink_phi_;

    // theta kink
    experimental_double kink_theta_;

  public:

    // joint status
    bool used_;

    // chi2 of connection along a sequence
    double chi2_;
    int32_t ndof_;
    double p_;

    //!Default constructor
    joint();

    //!Default destructor
    virtual ~joint();

    //! constructor
    joint(const experimental_point &epa, const experimental_point &epb, const experimental_point &epc);

    /*** dump ***/
    virtual void dump (std::ostream & a_out         = std::clog,
                       const std::string & a_title  = "",
                       const std::string & a_indent = "",
                       bool a_inherit          = false) const;



    //! set experimental_points
    void set(const experimental_point &epa, const experimental_point &epb, const experimental_point &epc);


    //! set first experimental_points
    void set_epa(const experimental_point &epa);

    //! set second experimental_points
    void set_epb(const experimental_point &epb);

    //! set third experimental_points
    void set_epc(const experimental_point &epc);

    //! set kink phi
    void set_kink_phi(const experimental_double &phi);

    //! set kink theta
    void set_kink_theta(const experimental_double &theta);

    //! set used
    void set_used(bool used);

    //! set chi2
    void set_chi2(double chi2);
    void set_ndof(int32_t ndof);
    void set_p(double p);

    //! get experimental_point a
    const experimental_point& epa()const;

    //! get experimental_point b
    const experimental_point& epb()const;

    //! get experimental_point c
    const experimental_point& epc()const;

    //! get kink phi
    const experimental_double& kink_phi()const;

    //! get kink theta
    const experimental_double& kink_theta()const;

    //! get used
    bool used() const;

    //! get chi2
    double chi2() const;
    int32_t ndof() const;
    double p() const;

    joint invert();

    bool operator<(const joint &j) const;

    double calculate_chi2(joint j, cell A, cell B, joint * modified, bool A_is_on_gap, bool B_is_on_gap)const;

  private:

    void calculate_kinks();

  };

} // namespace CAT

#endif // FALAISE_CAT_JOINT_H
