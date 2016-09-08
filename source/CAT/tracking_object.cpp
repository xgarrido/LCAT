/* -*- mode: c++ -*- */

#include <CAT/tracking_object.h>

// Third party:
//- GSL:
#include <gsl/gsl_cdf.h>


namespace CAT{
  namespace topology{

    void tracking_object::set_probmin( double probmin )
    {
      probmin_ = probmin;
    }

    double tracking_object::probmin() const
    {
      return probmin_;
    }


    double tracking_object::probof(double chi2, int ndof) const
    {
      return gsl_cdf_chisq_Q(chi2, ndof);
    }

    double tracking_object::get_probmin()const
    {
      return probmin_;
    }

  }
}
