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

#ifndef LCAT_UTILITIES
#define LCAT_UTILITIES

// Standard library:
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <limits.h>

namespace mybhep {

  static const double plus_infinity   = DBL_MAX;
  static const double small_neg       = -plus_infinity;
  static const int    default_integer = INT_MAX;
  static const double default_min     = plus_infinity;
  static const double default_max     = -default_min;

  /// Shift angle values in a limited range
  inline void fix_angles(double & a1_, double & a2_)
  {
    if (std::abs(a1_ - a2_) > M_PI) {
      if (a1_ < a2_) {
        a1_ += 2.*M_PI;
      } else {
        a2_ += 2.*M_PI;
      }
    }
    return;
  }

} // namespace mybhep

#endif // LCAT_UTILITIES
