// Ourselves:
#include <CAT/line.h>

// This project:
#include <CAT/utilities.h>

namespace CAT {

  namespace topology{

    using namespace std;

    //!Default constructor
    line::line(double probmin)
    {
      forward_axis_ = experimental_vector(mybhep::small_neg,mybhep::small_neg,mybhep::small_neg,
                                          mybhep::small_neg, mybhep::small_neg, mybhep::small_neg);
      used_ = false;
      set_probmin(probmin);
    }

    //!Default destructor
    line::~line()
    {
      return;
    }

    //! constructor
    line::line(const experimental_point & epa, const experimental_point & epb, double probmin)
    {
      set_probmin(probmin);
      epa_ = epa;
      epb_ = epb;
      set_forward_axis();
      used_ = false;
    }

    void line::dump (std::ostream & a_out       ,
                     const std::string & a_title,
                     const std::string & a_indent,
                     bool /* a_inherit */) const
    {
      std::string indent;
      if (! a_indent.empty ()) indent = a_indent;
      if (! a_title.empty ())
        {
          a_out << indent << a_title << std::endl;
        }

      a_out << indent << "line: -------------- " << std::endl;
      a_out << indent << " used: " << used() << std::endl;
      a_out << indent << " first point " << std::endl;
      this->epa().dump(a_out, "", indent + "    ");
      a_out << indent << " second point " << std::endl;
      this->epb().dump(a_out, "", indent + "    ");
      a_out << indent << " forward axis " << std::endl;
      this->forward_axis().dump(a_out, "", indent + "    ");
      a_out << indent << " -------------- " << std::endl;

      return;
    }

    //! set experimental_points
    void line::set(const experimental_point & epa, const experimental_point & epb)
    {
      epa_ = epa;
      epb_ = epb;
      set_forward_axis();
    }

    //! set first experimental_point
    void line::set_epa(const experimental_point & epa)
    {
      epa_ = epa;
    }

    //! set second experimental_point
    void line::set_epb(const experimental_point & epb)
    {
      epb_ = epb;
    }

    //! set used
    void line::set_used(bool used)
    {
      used_ = used;
    }

    //! get experimental_point a
    const experimental_point& line::epa()const
    {
      return epa_;
    }

    //! get experimental_point b
    const experimental_point& line::epb()const
    {
      return epb_;
    }

    //! get forward axis
    const experimental_vector& line::forward_axis()const
    {
      return forward_axis_;
    }

    //! get used
    bool line::used() const
    {
      return used_;
    }

    experimental_double line::phi(){
      return forward_axis_.phi();
    }

    experimental_double line::theta(){
      return forward_axis_.theta();
    }

    experimental_double line::kink_phi(const line & l){
      return forward_axis_.kink_phi(l.forward_axis());
    }

    experimental_double line::kink_theta(const line & l)
    {
      return forward_axis_.kink_theta(l.forward_axis());
    }

    double line::chi2(const line & l, bool use_theta_kink, double *chi2_just_phi)
    {

      experimental_double phi_kink = kink_phi(l);
      experimental_double theta_kink = kink_theta(l);

      double result = std::pow(phi_kink.value()/phi_kink.error(),2);

      *chi2_just_phi = result;

      if( use_theta_kink ){
        result += std::pow(theta_kink.value()/theta_kink.error(),2);
      }

      // if( print_level() > mybhep::VERBOSE ){
      //   std::clog << appname_ << " calculate chi2: phi kink : "; (phi_kink*180./M_PI).dump();
      //   if( use_theta_kink ){
      //     std::clog << " theta kink : "; (theta_kink*180./M_PI).dump();
      //   }
      //   std::clog << " chi2 " << result << std::endl;
      // }

      return result;
    }

    void line::set_a_forward_axis(const experimental_vector &  v)
    {
      forward_axis_ = v;
      return;
    }

    line line::invert()
    {
      line inverted(probmin());
      inverted.set_epa(epb());
      inverted.set_epb(epa());
      inverted.set_used(used());
      experimental_vector O(0.,0.,0.,0.,0.,0.);
      inverted.set_a_forward_axis(O-forward_axis());
      return inverted;
    }


    void line::set_forward_axis()
    {
      forward_axis_.set(epa(), epb());
      return;
    }

  } // namespace topology

} // namespace CAT
