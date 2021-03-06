// Ourselves:
#include <CAT/experimental_double.h>

// Standard library:
#include <cmath>
#include <limits>
#include <iomanip>

namespace CAT {

  bool experimental_double::is_valid () const
  {
    return is_value_valid () && is_error_valid ();
  }

  bool experimental_double::is_value_valid () const
  {
    return v_ == v_;
  }

  bool experimental_double::is_error_valid () const
  {
    return e_ == e_;
  }

  //!Default constructor
  experimental_double::experimental_double()
  {
    v_ = std::numeric_limits<double>::quiet_NaN ();
    e_ = std::numeric_limits<double>::quiet_NaN ();
    return;
  }

  //!Default destructor
  experimental_double::~experimental_double()
  {
    return;
  }

  //! constructor
  experimental_double::experimental_double(const double &v, const double &e)
  {
    v_ = v;
    e_ = e;
    return;
  }

  void experimental_double::dump (std::ostream & a_out,
                                  const std::string & /* a_title */,
                                  const std::string & /* a_indent */,
                                  bool /* a_inherit */) const
  {
    a_out << std::setprecision(15) << value() << " +- " << std::setprecision(15)  << error();
    return;
  }

  //! set value and error
  void experimental_double::set(const experimental_double &v)
  {
    v_ = v.value();
    e_ = v.error();
  }

  //! set value and error
  void experimental_double::set(const double &val, const double &err)
  {
    v_ = val;
    e_ = err;
  }

  //! set value
  void experimental_double::set_value(const double &v)
  {
    v_ = v;
  }

  //! set error
  void experimental_double::set_error(const double &e)
  {
    e_ = e;
  }

  //! get value
  const double& experimental_double::value() const
  {
    return v_;
  }

  //! get error
  const double& experimental_double::error() const
  {
    return e_;
  }

  // Operators
  //! operador +=
  experimental_double& experimental_double::operator += (const experimental_double& p2)
  {
    experimental_double& p1= *this;
    double val = p1.value() + p2.value();
    double err = std::hypot(p1.error(), p2.error());
    p1.set_value(val);
    p1.set_error(err);
    return p1;
  }

  //! operador -=
  experimental_double& experimental_double::operator -= (const experimental_double& p2)
  {
    experimental_double& p1= *this;
    double val = p1.value() - p2.value();
    double err = std::hypot(p1.error(), p2.error());
    p1.set_value(val);
    p1.set_error(err);

    return p1;
  }

  //! operador *=
  experimental_double& experimental_double::operator *= (experimental_double a)
  {
    experimental_double& p1= *this;
    double val = p1.value()*a.value();
    double err = std::hypot(a.value()*p1.error(), p1.value()*a.error());

    p1.set_value(val);
    p1.set_error(err);

    return p1;
  }

  //! operador *=
  experimental_double& experimental_double::operator *= (double a)
  {
    experimental_double& p1= *this;
    double val = p1.value()*a;
    double err = p1.error()*std::abs(a);

    p1.set_value(val);
    p1.set_error(err);

    return p1;
  }

  //! operador /=
  experimental_double& experimental_double::operator /= (experimental_double a)
  {
    experimental_double& p1= *this;

    if( a.value() == 0 ){
      std::clog << "CAT::experimental_double::operator/=: problem: division by experimental_double with value " << a.value() << std::endl;
    }

    double val = p1.value()/a.value();
    double err = std::hypot(p1.error()/a.value(),p1.value()*a.error()/std::pow(a.value(),2));
    p1.set_value(val);
    p1.set_error(err);
    return p1;
  }

  //! operador /=
  experimental_double& experimental_double::operator /= (double a)
  {
    experimental_double& p1= *this;

    if( a == 0 ){
      std::clog << "CAT::experimental_double::operator/=: problem: division by double with value " << a << std::endl;
    }

    double val = p1.value()/a;
    double err = p1.error()/std::abs(a);
    p1.set_value(val);
    p1.set_error(err);
    return p1;
  }

  // Operations with experimental_points
  // -v
  // sin(v)
  experimental_double experimental_sin (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(sin(v1.value()));
    v.set_error(std::abs(cos(v1.value()))*v1.error());
    return v;
  }

  // cos(v)
  experimental_double experimental_cos (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(cos(v1.value()));
    v.set_error(std::abs(sin(v1.value()))*v1.error());
    return v;
  }

  // tan(v)
  experimental_double experimental_tan (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(tan(v1.value()));
    v.set_error((1. + std::pow(v.value(),2))*v1.error());
    return v;
  }

  // asin(v)
  experimental_double experimental_asin (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(asin(v1.value()));
    v.set_error(v1.error()/std::sqrt(1 - std::pow(v1.value(),2)));
    return v;
  }

  // acos(v)
  experimental_double experimental_acos (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(acos(v1.value()));
    v.set_error(v1.error()/std::sqrt(1 - std::pow(v1.value(),2)));
    return v;
  }

  // atan(v)
  experimental_double experimental_atan2 (const experimental_double& v1,
                                          const experimental_double& v2)
  {
    experimental_double v;
    v.set_value(atan2(v1.value(), v2.value()));

    if( v2.value() == 0. ){  // if angle = 90 degrees,
      // obtain error from angle = 180 - (other angle)
      double den = 1 + std::pow(v2.value()/v1.value(),2);
      double num = std::pow(v2.error()/v1.value(),2) + std::pow(v2.value()*v1.error()/std::pow(v1.value(),2),2);
      v.set_error(std::sqrt(num)/den);
    }
    else{

#if 0
      if( std::abs(v1.error()) > std::abs(v1.value()) ){
        // large error propagation for the arctan !
        // when "y" has large error, standard error propagation under-estimates the error
        double r1 = atan2(v1.value() - v1.error(), v2.value());
        double r2 = atan2(v1.value() + v1.error(), v2.value());
        v.set_error(std::abs(r1 - r2)/2.);
        return v;
      }
#endif

      double den = 1 + std::pow(v1.value()/v2.value(),2);
      double num = std::pow(v1.error()/v2.value(),2) + std::pow(v1.value()*v2.error()/std::pow(v2.value(),2),2);
      v.set_error(std::sqrt(num)/den);

    }
    return v;
  }

  // square(v)
  experimental_double experimental_square (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(std::pow(v1.value(),2));
    v.set_error(2*std::abs(v1.value())*v1.error());
    return v;
  }

  // sqrt(v)
  experimental_double experimental_sqrt (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(std::sqrt(v1.value()));
    v.set_error(v1.error()/(2*v.value()));
    return v;
  }

  // cube(v)
  experimental_double experimental_cube (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(std::pow(v1.value(),3));
    v.set_error(3*std::pow(v1.value(),2)*v1.error());
    return v;
  }

  // std::abs(v)
  experimental_double experimental_fabs (const experimental_double& v1)
  {
    experimental_double v;
    v.set_value(std::abs(v1.value()));
    v.set_error(v1.error());
    return v;
  }

  experimental_double operator - (const experimental_double& v1)
  {
    experimental_double v = v1;
    v.set_value(-v1.value());
    v.set_error(v1.error());

    return v;
  }

  // v1+v2
  experimental_double operator + (const experimental_double& v1, const experimental_double& v2)
  {
    experimental_double v = v1;
    v+=v2;
    return v;
  }

  //! v1-v2
  experimental_double operator - (const experimental_double& v1, const experimental_double& v2)
  {
    experimental_double v = v1;
    v-=v2;
    return v;
  }

  // v*d
  experimental_double operator * (const experimental_double& v1, const experimental_double& d)
  {
    experimental_double v = v1;
    v*=d;
    return v;
  }

  // v/d
  experimental_double operator / (const experimental_double& v1, const experimental_double& d)
  {
    experimental_double v = v1;
    v/=d;
    return v;
  }

  // v*d
  experimental_double operator * (const experimental_double& v1, double d)
  {
    experimental_double v = v1;
    experimental_double dd(d,0.);
    v*=dd;
    return v;
  }

  // v/d
  experimental_double operator / (const experimental_double& v1, double d)
  {
    experimental_double v = v1;
    experimental_double dd(d,0.);
    v/=dd;
    return v;
  }

  // v/d
  experimental_double operator / (double d, const experimental_double& v1)
  {
    experimental_double dd(d,0.);
    experimental_double v = v1;
    dd/=v;
    return dd;
  }

  // average
  experimental_double average(const std::vector<experimental_double> & vs_)
  {
    if (vs_.empty()){
      experimental_double bad;
      std::clog << "CAT::average: problem: averaging over an empty vector !" << std::endl;
      return bad;
    }

    double mean = 0.;
    double err = 0.;

    for(std::vector<experimental_double>::const_iterator iv=vs_.begin(); iv!=vs_.end(); ++iv){
      mean += iv->value();
      err += std::pow(iv->error(),2);
    }

    return experimental_double(mean/vs_.size(), std::sqrt(err)/vs_.size());
  }

  // weighted average
  experimental_double weighted_average (const std::vector<experimental_double> & vs_)
  {
    double mean = 0.;
    double inverr = 0.;
    double newerr = 0.;
    for(std::vector<experimental_double>::const_iterator iv=vs_.begin(); iv!=vs_.end(); ++iv){
      if( iv->error() ){
        mean += iv->value()/std::pow(iv->error(),2);
        inverr += 1/std::pow(iv->error(),2);
        newerr += std::pow(iv->error(),2);
      }else{
        std::clog << "CAT::weighted_average: problem: double has error zero; switch to normal average " << std::endl;
        return average(vs_);
      }
    }
    return experimental_double(mean/inverr, sqrt(newerr));
  }

} // namespace CAT
