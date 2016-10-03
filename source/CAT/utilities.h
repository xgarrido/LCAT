// -*- mode: c++ -*-

/*
 * Copyright (C) 2002 J.J. Gomez-Cadenas, J.A. Hernando
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef FALAISE_CAT_UTILITIES_H
#define FALAISE_CAT_UTILITIES_H 1

// Standard library:
#include <limits>

namespace CAT {

  static const double plus_infinity   = std::numeric_limits<double>::max();
  static const double small_neg       = -plus_infinity;
  static const int    default_integer = std::numeric_limits<int>::max();
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

  /// Shift angle values in a limited range
  void fix_angles(double & a1_, double & a2_);

  /// Return probability given chi2 value and number of degree of freedom
  double probof(double chi2_, int ndof_);

} // namespace CAT

#endif // FALAISE_CAT_UTILITIES_H
