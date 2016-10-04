/* -*- mode: c++ -*- */
#ifndef FALAISE_CAT_BROKEN_LINE_H
#define FALAISE_CAT_BROKEN_LINE_H 1

// Standard library
#include <iostream>

// This project
#include <CAT/experimental_point.h>

namespace CAT {


  /// \brief A broken_line is composed of a vector of experimental points
  class broken_line
  {
  public:

    // points
    std::vector<experimental_point> eps_;

    // chi2 of connection along a sequence
    double chi2_;
    int32_t ndof_;
    double p_;
    size_t ifirst_;
    size_t ilast_;

    //!Default constructor
    broken_line();

    //!Default destructor
    virtual ~broken_line();

    //! constructor
    broken_line(const std::vector<experimental_point> &eps);

    /*** dump ***/
    virtual void dump (std::ostream & a_out         = std::clog,
                       const std::string & a_title  = "",
                       const std::string & a_indent = "",
                       bool a_inherit          = false) const;



    //! set experimental_points
    void set(const std::vector<experimental_point> &eps);

    //! set chi2
    void set_chi2(double chi2);
    void set_ndof(int32_t ndof);
    void set_p(double p);
    void set_ifirst(size_t i);
    void set_ilast(size_t i);

    //! get experimental_points
    const std::vector<experimental_point> & eps()const;

    //! get chi2
    double chi2() const;
    int32_t ndof() const;
    double p() const;
    size_t ifirst() const;
    size_t ilast() const;

    void calculate_chi2();


  };

}
#endif // FALAISE_CAT_BROKEN_LINE_H
