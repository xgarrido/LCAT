#include <CAT/cell_base.h>
#include <cmath>

namespace CAT{
namespace topology{

  using namespace mybhep;

  experimental_point cell::build_from_cell(const experimental_vector & forward_,
                                           const experimental_vector & transverse_,
                                           const experimental_double & cos_,
                                           int sign_,
                                           bool replace_r_,
                                           double max_r_) const
  {
    const experimental_double sin = experimental_sin(experimental_acos(cos_))*sign_;
    experimental_double radius = r();
    if (replace_r_) {
      radius.set_value(max_r_);
    }

    const experimental_point p = (experimental_vector(ep()) + radius*(forward_*cos_ + transverse_*sin)).point_from_vector();
    return p;
  }


  experimental_point cell::angular_average(const experimental_point & epa_, const experimental_point & epb_, experimental_double & angle_)
  {
    if (small()){
      angle_.set_value(0.);
      angle_.set_error(0.1); // fictitious value to avoid divergences
      return ep();
    }

    const experimental_vector v1(ep(), epa_);
    const experimental_vector v2(ep(), epb_);

    experimental_double phi1 = v1.phi();
    experimental_double phi2 = v2.phi();

    double rephi1 = phi1.value();
    double rephi2 = phi2.value();

    fix_angles(rephi1, rephi2);

    phi1.set_value(rephi1);
    phi2.set_value(rephi2);
    angle_ = phi1 - phi2;

    std::vector<experimental_double> phis;
    phis.push_back(phi1);
    phis.push_back(phi2);

    const experimental_double ave_phi = weighted_average(phis);

    // if( print_level() >= mybhep::VVERBOSE ){
    //   std::clog << "CAT::cell::angular_average: averaging phi1: "; (phi1*180./M_PI).dump();
    //   std::clog << " and phi2: "; (phi2*180./M_PI).dump();
    //   std::clog << " to phi_ave: "; (ave_phi*180./M_PI).dump();
    //   std::clog << " " << std::endl;
    // }

    const experimental_double cos_ave_phi = experimental_cos(ave_phi);

    const experimental_vector x(1.,0.,0.,0.,0.,0.);
    const experimental_vector z(0.,0.,1.,0.,0.,0.);
    int sign = 1;
    if( ave_phi.value() < 0. ) sign = -1;  // ave_phi in [180,360], i.e. p to the left of cell center


    // distance of each point from center of cell
    // not necessarily the same if one of the 2 points results from intersecting cells
    // in such case, keep largest of 2 values to locate average
    // otherwise it will not make sense with the other point from intersecting cells
    const double r1 = v1.hor().length().value();
    const double r2 = v2.hor().length().value();
    bool replace_r = false;
    double maxr = 0.;
    if( r1 != r2 ){
      replace_r = true;
      maxr = std::max(r1,r2);
    }

    const double errx = std::max(std::abs(epb_.x().value() - epa_.x().value()), (epb_+epa_).x().error());
    const double errz = std::max(std::abs(epb_.z().value() - epa_.z().value()), (epb_+epa_).z().error());

    experimental_point p = build_from_cell(x, z, cos_ave_phi, sign, replace_r, maxr);
    p.set_ex(errx);
    p.set_ez(errz);
    return p;
  }


  bool cell::same_quadrant(const experimental_point & epa_, const experimental_point & epb_) const
  {
    // check if the angular position of points epa and epb
    // is less than 90 degrees apart around the cell

    experimental_double initial_phi1 = experimental_vector(ep(), epa_).phi();
    experimental_double initial_phi2 = experimental_vector(ep(), epb_).phi();

    double re_initial_phi1 = initial_phi1.value();
    double re_initial_phi2 = initial_phi2.value();
    fix_angles(re_initial_phi1, re_initial_phi2);

    if (std::abs(re_initial_phi1 - re_initial_phi2) > M_PI/2.) return false;
    return true;
  }

  bool cell::same_cell(topology::cell c)const{

    if( block() == c.block() &&
        layer() == c.layer() &&
        iid() == c.iid())
      return true;

    return false;

  }

  bool cell::intersect(topology::cell c) const{

    double fraction_limit = 0.9; /// fraction of radius after which cells intersect

    double dist = experimental_vector(ep(), c.ep()).hor().length().value();
    experimental_double rsum = r() + c.r();

    if( rsum.value() > dist*fraction_limit ){
      // if( print_level() >= mybhep::VVERBOSE ){
      //   std::clog << "CAT::cell::intersect: cells " << id() << " and " << c.id() << " intersect: dist " << dist << " radii " << r().value() << " and " << c.r().value() << " rsum " << rsum.value() << std::endl;
      // }
      return true;
    }

    bool test = false;

    /*
    double chi2 = square((rsum.value() - dist)/rsum.error());

    bool test = probof(chi2, 1) > probmin();
    if( test )
      if( print_level() >= mybhep::VVERBOSE ){
        std::clog <<" cells " << id() << " and " <<c.id() << " intersect: dist " << dist << " radii " << r().value() << " and " <<c.r().value()<< " rsum " << rsum.value() << " +- " << rsum.error() << " prob " << probof(chi2, 1) << " chi2 " << chi2 << " problimit " << probmin() << std::endl;
      }
    */

    return test;

  }


}
}
