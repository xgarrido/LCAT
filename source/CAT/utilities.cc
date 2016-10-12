// Ourselves
#include <CAT/utilities.h>

// Standard library:
#include <cmath>

// Third party
// - Boost:
#include <boost/scoped_ptr.hpp>
// - GSL:
#include <gsl/gsl_cdf.h>

namespace CAT {

  // static
  const constants & constants::instance ()
  {
    static boost::scoped_ptr<constants> g_global_constants(0);
    if (g_global_constants.get () == 0) {
      g_global_constants.reset(new constants);
    }
    return *g_global_constants.get();
  }

  constants::constants ()
  {
    _min_probability_ = 1e-200;// need units
    return;
  }

  double constants::get_minimal_probability() const
  {
    return _min_probability_;
  }


  void invalidate(int32_t & i_)
  {
    i_ = invalid_integer();
    return;
  }

  int32_t invalid_integer()
  {
    return std::numeric_limits<int32_t>::max();
  }

  double probof(double chi2_, int ndof_)
  {
    double p = 0.;
    if (ndof_ > 0 && chi2_ > 0.0)
      p = gsl_cdf_chisq_Q(chi2_, ndof_);
    return p;
  }

  void fix_angles(double & a1_, double & a2_)
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
}
