// -*- mode: c++ -*-
#ifndef FALAISE_CAT_UTILITIES_H
#define FALAISE_CAT_UTILITIES_H 1

// Standard library:
#include <limits>
// Third party:
// - Boost:
#include <boost/cstdint.hpp>

namespace CAT {

  static const double plus_infinity   = std::numeric_limits<double>::max();
  static const double small_neg       = -plus_infinity;
  static const double default_min     = plus_infinity;
  static const double default_max     = -default_min;

  // \brief Singleton holding constants
  class constants {
  public:
    double get_minimal_probability() const;
    constants();
    static const constants & instance();
  private:
    double _min_probability_;
  };

  /// Invalidate integer value
  void invalidate(int32_t &);

  /// Return invalid value for integer
  int32_t invalid_integer();

  /// Shift angle values in a limited range
  void fix_angles(double & a1_, double & a2_);

  /// Return probability given chi2 value and number of degree of freedom
  double probof(double chi2_, int ndof_);

} // namespace CAT

#endif // FALAISE_CAT_UTILITIES_H
