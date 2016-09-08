/* -*- mode: c++ -*- */

#ifndef __CATAlgorithm__tracking_object_h
#define __CATAlgorithm__tracking_object_h 1

#include <string>
#include <vector>
#include <sstream>
#include <cmath>

#include <CAT/printable.h>


namespace CAT{
  namespace topology{


    // a generic tracking object
    class tracking_object : public printable
    {


    public:


      // the minimum probabilty
      double probmin_;


      void set_probmin( double probmin );

      double probmin() const;

      double probof(double chi2, int ndof)const;

      double get_probmin()const;

    };

  }
}

#endif // __CATAlgorithm__tracking_object_h
