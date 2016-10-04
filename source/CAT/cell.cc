// Ourselves
#include <CAT/cell.h>
#include <CAT/utilities.h>

// Standard library
#include <cmath>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/i_tree_dump.h>
#include <bayeux/datatools/clhep_units.h>

namespace CAT {

  cell::cell()
  {
    //ep_ = experimental_point();
    //r0_= experimental_double();
    //r_= experimental_double();
    _id_ = default_integer;
    _side_ = default_integer;
    _layer_ = default_integer;
    _row_ = default_integer;
    _prompt_ = true;
    _free_ = false;
    _begun_ = false;
    _small_radius_= 0.;
    return;
  }

  cell::~cell()
  {
    reset();
    return;
  }

  void cell::reset()
  {
    return;
  }

  void cell::tree_dump(std::ostream & out_,
                       const std::string & title_,
                       const std::string & indent_,
                       bool inherit_) const
  {
    std::string indent;
    if (! indent_.empty()) indent = indent_;
    if (! title_.empty()) {
      out_ << indent << title_ << std::endl;
    }

    out_ << indent << datatools::i_tree_dumpable::tag << "ID     : " << get_id() << std::endl;
    out_ << indent << datatools::i_tree_dumpable::tag << "Side   : " << get_side() << std::endl;
    out_ << indent << datatools::i_tree_dumpable::tag << "Row    : " << get_row() << std::endl;
    out_ << indent << datatools::i_tree_dumpable::tag << "Layer  : " << get_layer() << std::endl;
    out_ << indent << datatools::i_tree_dumpable::tag << "Prompt : " << is_prompt() << std::endl;
    out_ << indent << datatools::i_tree_dumpable::tag << "Small  : " << is_small() << std::endl;
    const experimental_point & ep = get_position();
    out_ << indent << datatools::i_tree_dumpable::tag << "Position (x,y,z) : ("
         << ep.x().value()/CLHEP::mm << ", " << ep.y().value()/CLHEP::mm << ", " << ep.z().value()/CLHEP::mm << ") mm" << std::endl;
    if (is_small() && is_prompt()) {
      out_ << indent << datatools::i_tree_dumpable::tag << "Original radius : " << get_original_radius().value()/CLHEP::mm << " mm" << std::endl;
    }
    out_ << indent << datatools::i_tree_dumpable::inherit_tag(inherit_) << "Radius : " << get_radius().value()/CLHEP::mm << " mm" << std::endl;
    return;
  }

  void cell::set_position(const experimental_point & ep_)
  {
    _position_ = ep_;
    return;
  }

  void cell::set_radius(double r_)
  {
    _r0_.set_value(r_);
    _set_radius_();
    return;
  }

  void cell::set_radius_error(double er_)
  {
    _r0_.set_error(er_);
    _set_radius_();
    return;
  }

  void cell::set_small_radius(double sr_)
  {
    _small_radius_ = sr_;
    return;
  }

  void cell::set_id(int id_)
  {
    _id_ = id_;
    return;
  }

  void cell::set_layer(int layer_)
  {
    _layer_ = layer_;
    return;
  }

  void cell::set_side(int side_)
  {
    _side_ = side_;
    return;
  }

  void cell::set_row(int row_)
  {
    _row_ = row_;
    return;
  }

  void cell::set_prompt(bool p_)
  {
    _prompt_ = p_;
    return;
  }

  void cell::set_free(bool free_)
  {
    _free_ = free_;
    return;
  }

  void cell::set_begun(bool begun_)
  {
    _begun_ = begun_;
  }

  bool cell::is_small() const
  {
    if (_r0_.value() <= _small_radius_) return true;
    return false;
  }

  const experimental_point & cell::get_position() const
  {
    return _position_;
  }

  const experimental_double & cell::get_radius() const
  {
    return _r_;
  }

  const experimental_double & cell::get_original_radius() const
  {
    return _r0_;
  }

  int cell::get_id() const
  {
    return _id_;
  }

  int cell::get_layer() const
  {
    return _layer_;
  }

  int cell::get_side() const
  {
    return _side_;
  }

  int cell::get_row() const
  {
    return _row_;
  }

  bool cell::is_prompt() const
  {
    return _prompt_;
  }

  bool cell::is_free() const
  {
    return _free_;
  }

  bool cell::begun() const
  {
    return _begun_;
  }

  void cell::_set_radius_()
  {
    _r_ = _r0_;
    return;
  }


  experimental_point cell::build_from_cell(const experimental_vector & forward_,
                                           const experimental_vector & transverse_,
                                           const experimental_double & cos_,
                                           int sign_,
                                           bool replace_r_,
                                           double max_r_) const
  {
    const experimental_double sin = experimental_sin(experimental_acos(cos_))*sign_;
    experimental_double radius = get_radius();
    if (replace_r_) {
      radius.set_value(max_r_);
    }

    const experimental_point p = (experimental_vector(get_position()) + radius*(forward_*cos_ + transverse_*sin)).point_from_vector();
    return p;
  }


  experimental_point cell::angular_average(const experimental_point & epa_, const experimental_point & epb_, experimental_double & angle_)
  {
    if (is_small()){
      angle_.set_value(0.);
      angle_.set_error(0.1); // fictitious value to avoid divergences
      return get_position();
    }

    const experimental_vector v1(get_position(), epa_);
    const experimental_vector v2(get_position(), epb_);

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

    // if( print_level() >= VVERBOSE ){
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

    experimental_double initial_phi1 = experimental_vector(get_position(), epa_).phi();
    experimental_double initial_phi2 = experimental_vector(get_position(), epb_).phi();

    double re_initial_phi1 = initial_phi1.value();
    double re_initial_phi2 = initial_phi2.value();
    fix_angles(re_initial_phi1, re_initial_phi2);

    if (std::abs(re_initial_phi1 - re_initial_phi2) > M_PI/2.) return false;
    return true;
  }

  bool cell::intersect(const cell & c_) const
  {
    const double dist = experimental_vector(get_position(), c_.get_position()).hor().length().value();
    const experimental_double rsum = get_radius() + c_.get_radius();

    const double fraction_limit = 0.9; /// fraction of radius after which cells intersect
    if( rsum.value() > dist*fraction_limit ){
      // if( print_level() >= VVERBOSE ){
      //   std::clog << "CAT::cell::intersect: cells " << id() << " and " << c.id() << " intersect: dist " << dist << " radii " << r().value() << " and " << c.r().value() << " rsum " << rsum.value() << std::endl;
      // }
      return true;
    }

    return false;
  }
}
