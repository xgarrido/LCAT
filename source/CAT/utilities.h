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

  enum prlevel{MUTE,CONCISE,NORMAL,WARNING,DETAILED,VERBOSE,VVERBOSE,DUMP};

  inline prlevel get_info_level(std::string info)
  {
    if(info == "MUTE")  return MUTE;
    if(info == "CONCISE")  return CONCISE;
    if(info == "NORMAL")  return NORMAL;
    if(info == "WARNING")  return WARNING;
    if(info == "DETAILED")  return DETAILED;
    if(info == "VERBOSE")  return VERBOSE;
    if(info == "VVERBOSE")  return VVERBOSE;
    if(info == "DUMP")  return DUMP;
    else return NORMAL;
  }

  static const double plus_infinity   = DBL_MAX;
  static const double small_neg       = -plus_infinity;
  static const size_t default_integer = INT_MAX;
  static const double default_min     = plus_infinity;
  static const double default_max     = -default_min;

  /// Square function
  template <class T>
  inline T square(T a)
  {
    return a*a;
  }

  /// Cubic power function
  template <class T>
  inline T cube(T a)
  {
    return a*a*a;
  }

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

  /// Round up a double type to show n decimals
  inline double round_to(double value, int to = 1)
  {
    double places = pow(10.0, to);
    return round(value * places) / places;
  }

} // namespace mybhep

#endif // LCAT_UTILITIES
