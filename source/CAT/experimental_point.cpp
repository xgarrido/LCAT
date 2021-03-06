// Ourselves:
#include <CAT/experimental_point.h>
// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/clhep_units.h>

// This project:
#include <CAT/utilities.h>

namespace CAT {

  //!Default constructor
  experimental_point::experimental_point()
  {
    x_ = experimental_double();
    y_ = experimental_double();
    z_ = experimental_double();
    radius_ = experimental_double();
  }

  //!Default destructor
  experimental_point::~experimental_point()
  {
  }

  //! constructor
  experimental_point::experimental_point(const experimental_double &x, const experimental_double &y, const experimental_double &z){
    x_ = x;
    y_ = y;
    z_ = z;
    set_radius();
  }

  //! constructor from coordinates with error
  experimental_point::experimental_point(double x, double y, double z, double ex, double ey, double ez){
    x_.set_value(x);
    y_.set_value(y);
    z_.set_value(z);
    x_.set_error(ex);
    y_.set_error(ey);
    z_.set_error(ez);
    set_radius();
  }

  void experimental_point::dump (std::ostream & a_out,
                                 const std::string & a_title  ,
                                 const std::string & a_indent ,
                                 bool /* a_inherit */  ) const
  {
    std::string indent;
    if (! a_indent.empty ()) indent = a_indent;
    if (! a_title.empty ())
      {
        a_out << indent << a_title << std::endl;
      }

    a_out << indent << " x : "; (x()/CLHEP::mm).dump(); a_out << " [mm] " << std::endl;
    a_out << indent << " y : "; (y()/CLHEP::mm).dump(); a_out << " [mm] " << std::endl;
    a_out << indent << " z : "; (z()/CLHEP::mm).dump(); a_out << " [mm] " << std::endl;

    return;
  }


  //! set x
  void experimental_point::set_x(const experimental_double &x)
  {
    x_.set_value(x.value());
    x_.set_error(x.error());
    set_radius();
  }

  //! set ex
  void experimental_point::set_ex(double ex)
  {
    x_.set_error(ex);
    set_radius();
  }

  //! set y
  void experimental_point::set_y(const experimental_double &y)
  {
    y_.set_value(y.value());
    y_.set_error(y.error());
  }

  //! set ey
  void experimental_point::set_ey(double ey)
  {
    y_.set_error(ey);
    set_radius();
  }

  //! set z
  void experimental_point::set_z(const experimental_double &z)
  {
    z_.set_value(z.value());
    z_.set_error(z.error());
    set_radius();
  }

  //! set ez
  void experimental_point::set_ez(double ez)
  {
    z_.set_error(ez);
    set_radius();
  }

  //!get experimental x
  const experimental_double& experimental_point::x() const
  {
    return x_;
  }

  //!get experimental y
  const experimental_double& experimental_point::y() const
  {
    return y_;
  }

  //!get experimental z
  const experimental_double& experimental_point::z() const
  {
    return z_;
  }

  //!get experimental radius
  const experimental_double& experimental_point::radius() const
  {
    return radius_;
  }

  //! distance
  experimental_double experimental_point::distance(const experimental_point& p2) const
  {
    experimental_double result;

    result.set_value(std::sqrt(std::pow(x_.value()-p2.x().value(),2) +
                               std::pow(y_.value()-p2.y().value(),2) +
                               std::pow(z_.value()-p2.z().value(),2)));
    result.set_error(std::sqrt(std::pow(x_.value()*x_.error(),2) + std::pow(p2.x().value()*p2.x().error(),2) +
                               std::pow(y_.value()*y_.error(),2) + std::pow(p2.y().value()*p2.y().error(),2) +
                               std::pow(z_.value()*z_.error(),2) + std::pow(p2.z().value()*p2.z().error(),2))/result.value());

    return result;
  }


  //! distance
  experimental_double experimental_point::hor_distance(const experimental_point& p2) const
  {
    experimental_double result;

    result.set_value(std::hypot(x_.value()-p2.x().value(),z_.value()-p2.z().value()));
    result.set_error(std::sqrt(std::pow(x_.value()*x_.error(),2) + std::pow(p2.x().value()*p2.x().error(),2) +
                               std::pow(z_.value()*z_.error(),2) + std::pow(p2.z().value()*p2.z().error(),2))/result.value());

    return result;
  }


  void experimental_point::set_radius(){
    // propagate radius error:
    //  r = sqrt(x^2 + z^2)
    //  dr/dx = x/r,
    //  dr/dz = z/r

    double rr = std::hypot(x_.value(), z_.value());
    if( std::isnan(rr) ) rr = small_neg;
    double err = std::hypot(x_.value()*x_.error(), z_.value()*z_.error())/rr;
    if( std::isnan(err) ) err = small_neg;

    radius_.set_value(rr);
    radius_.set_error(err);

    return;
  }

  experimental_point average (const std::vector<experimental_point> &vs)
  {
    std::vector<experimental_double> xs;
    std::vector<experimental_double> ys;
    std::vector<experimental_double> zs;

    for(std::vector<experimental_point>::const_iterator iv=vs.begin(); iv!=vs.end(); ++iv){
      xs.push_back(iv->x());
      ys.push_back(iv->y());
      zs.push_back(iv->z());
    }

    return experimental_point(average(xs), average(ys), average(zs));
  }

} // namespace CAT
