/* -*- mode: c++ -*- */

#include <CATAlgorithm/cell_triplet.h>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/clhep_units.h>

namespace CAT{
  namespace topology{

    //!Default constructor
    cell_triplet::cell_triplet()
    {
      appname_= "cell_triplet: ";
      free_ = false;
      begun_ = false;
    }

    //!Default destructor
    cell_triplet::~cell_triplet(){}

    //! constructor
    cell_triplet::cell_triplet(const cell &ca, const cell &cb, const cell &cc, double probmin){
      set_probmin(probmin);
      appname_= "cell_triplet: ";
      ca_ = ca;
      cb_ = cb;
      cc_ = cc;
      free_ = false;
      begun_ = false;
      chi2s_.clear();
      probs_.clear();
    }

    /*** dump ***/
    void cell_triplet::dump (std::ostream & a_out,
                             const std::string & a_title,
                             const std::string & a_indent,
                             bool /* a_inherit */) const{
      std::string indent;
      if (! a_indent.empty ()) indent = a_indent;
      if (! a_title.empty ())
        {
          a_out << indent << a_title << std::endl;
        }

      a_out << indent << appname_ << " -------------- " << std::endl;
      a_out << indent  << " free: " << free() << " begun: " << begun()  << std::endl;
      a_out << indent  << " first cell " << std::endl;
      ca().dump(a_out,"", indent + "   ");
      a_out << indent << " second cell " << std::endl;
      cb().dump(a_out, "",indent + "   ");
      a_out << indent << " third cell " << std::endl;
      cc().dump(a_out, "",indent + "   ");
      for(std::vector<joint>::const_iterator ijoint=joints_.begin(); ijoint!=joints_.end(); ++ijoint)
        ijoint->dump(a_out,"",indent + "   ");
      a_out << indent  << " -------------- " << std::endl;

      return;
    }

    /*** dump ***/
    void cell_triplet::dump_joint (joint j,
                                   std::ostream & a_out,
                                   const std::string & a_title,
                                   const std::string & a_indent,
                                   bool /* a_inherit*/) const{
      std::string indent;
      if (! a_indent.empty ()) indent = a_indent;
      if (! a_title.empty ())
        {
          a_out << indent << a_title << std::endl;
        }

      a_out << indent << appname_ << " -------------- " << std::endl;
      a_out << indent << " first point " << std::endl;
      this->ca().dump_point(j.epa());
      a_out << indent << " second point " << std::endl;
      this->cb().dump_point(j.epb());
      a_out << indent << " third point " << std::endl;
      this->cc().dump_point(j.epc());
      a_out << indent << " -------------- " << std::endl;

      return;
    }

    //! set cells
    void cell_triplet::set(const cell_couplet &cca, const cell_couplet &ccb){
      cb_ = cca.ca();
      ca_ = cca.cb();
      cc_ = ccb.cb();
      if( cca.ca().id() != ccb.ca().id() ){
        std::clog << " problem: trying to form a triplet of cell with cells " << cca.ca().id() << " "
             << cca.cb().id() << " "
             << ccb.ca().id() << " "
             << ccb.cb().id() << std::endl;
      }
    }

    //! set cells
    void cell_triplet::set(const cell &ca, const cell &cb, const cell &cc){
      ca_ = ca;
      cb_ = cb;
      cc_ = cc;
    }


    //! set free level
    void cell_triplet::set_free(bool free){
      free_ = free;
    }

    //! set begun level
    void cell_triplet::set_begun(bool begun){
      begun_ = begun;
    }

    //! set joints
    void cell_triplet::set_joints(const std::vector<joint> & joints){
      joints_ = joints;
    }

    //! set chi2 list
    void cell_triplet::set_chi2s(const std::vector<double> & chi2s){
      chi2s_ = chi2s;
    }

    //! set prob list
    void cell_triplet::set_probs(const std::vector<double> & probs){
      probs_ = probs;
    }

    //! get first cell couplet
    cell_couplet cell_triplet::cca()
    {
      cell_couplet cc1(cb_, ca_, probmin());
      return cc1;
    }

    //! get second cell couplet
    cell_couplet cell_triplet::ccb()
    {
      cell_couplet cc2(cb_, cc_, probmin());
      return cc2;
    }

    //! get joints
    const std::vector<joint>& cell_triplet::joints() const
    {
      return joints_;
    }

    //! get first cell
    const cell& cell_triplet::ca()const
    {
      return ca_;
    }

    //! get second cell
    const cell& cell_triplet::cb()const
    {
      return cb_;
    }

    //! get third cell
    const cell& cell_triplet::cc()const
    {
      return cc_;
    }

    //! get list of chi2
    const std::vector<double>& cell_triplet::chi2s() const
    {
      return chi2s_;
    }

    //! get list of prob
    const std::vector<double>& cell_triplet::probs() const
    {
      return probs_;
    }


    //! get free level
    bool cell_triplet::free()const{
      return free_;
    }

    //! get begun level
    bool cell_triplet::begun()const{
      return begun_;
    }

    void cell_triplet::calculate_joints(double ratio_, double phi_limit_)
    {
      datatools::logger::priority local_priority = datatools::logger::PRIO_WARNING;

      DT_LOG_DEBUG(local_priority, "Calculate joints for cells: " << ca_.id() << " " << cb_.id() << " " << cc_.id());

      joints_.clear();
      std::vector<line> t1 = cca().tangents(); // note: this tangent goes from cell B to cell A
      std::vector<line> t2 = ccb().tangents(); // this goes from B to C
      bool intersect_ab = ca_.intersect(cb_);
      bool intersect_bc = cb_.intersect(cc_);
      bool intersect_ca = cc_.intersect(ca_);

      const bool is_fast = ca_.fast();
      if (! is_fast) {
	phi_limit_ = std::max(phi_limit_, 90. * CLHEP::degree);
        DT_LOG_DEBUG(local_priority, "Cells are delayed, reset phi_limit " << phi_limit_/CLHEP::degree);
      }

      if (local_priority >= datatools::logger::PRIO_DEBUG) {
        DT_LOG_DEBUG(local_priority, "Angles of tangents " << ca_.id() << " -> " << cb_.id() << " :");
        for (std::vector<line>::const_iterator i1 = t1.begin(); i1!=t1.end(); ++i1){
          std::clog << i1 - t1.begin() << ":  phi "; ca_.dump_point_phi(i1->epb()); std::clog << " -> "; cb_.dump_point_phi(i1->epa()); std::clog << " " << std::endl;
        }
        std::clog << appname_ << " angles of tangents " << cb_.id() << " -> " << cc_.id() << " :" << std::endl;
        for(std::vector<line>::const_iterator i2=t2.begin(); i2!=t2.end(); ++i2){
          std::clog << i2 - t2.begin() << ":  phi ";  cb_.dump_point_phi(i2->epa()); std::clog << " -> " ; cc_.dump_point_phi(i2->epb()); std::clog << " " << std::endl;
        }
        if( ca_.small() ) std::clog << " cell " << ca_.id() << " is small " << std::endl;
        if( cb_.small() ) std::clog << " cell " << cb_.id() << " is small " << std::endl;
        if( cc_.small() ) std::clog << " cell " << cc_.id() << " is small " << std::endl;
        if( intersect_ab ) std::clog << " cells " << ca_.id() << " and " << cb_.id() << " intersect " << std::endl;
        if( intersect_bc ) std::clog << " cells " << cb_.id() << " and " << cc_.id() << " intersect " << std::endl;
        if( intersect_ca ) std::clog << " cells " << cc_.id() << " and " << ca_.id() << " intersect " << std::endl;
      }

      size_t idx1 = 0;
      size_t idx2 = 0;
      for (const auto & i1 : t1) {
        for (const auto & i2 : t2) {
          DT_LOG_DEBUG(local_priority, "Tangents " << idx1++ << " and " << idx2++);

          bool shall_include_separation = true;
          experimental_vector a0 = i1.forward_axis();
          experimental_vector a = a0.hor();
          experimental_vector d0 = i2.forward_axis();
          experimental_vector d = d0.hor();

          // keep only the connections that don't invert foward sense
          const double psc = std::abs(a.kink_phi(d).value());
          if (psc < 60.*CLHEP::degree) {
            DT_LOG_DEBUG(local_priority, "Rejected because direction is reversed: psc = " << psc/CLHEP::degree << "°");
            continue;
          }

          // Middle cell has small radius (< 2 mm)
          if (cb_.small()) {
            DT_LOG_DEBUG(local_priority, "No separation: middle cells is small");
            shall_include_separation = false;
          } else if (intersect_ab) {
            experimental_vector b0 = cca().forward_axis();
            experimental_vector b = b0.hor();
            const double psc1 = a.kink_phi(b).value();

            if (std::abs(psc1 - 90.*CLHEP::degree) < 30.*CLHEP::degree ||
                std::abs(psc1 + 90.*CLHEP::degree) < 30.*CLHEP::degree ||
                std::abs(psc1 - 270.*CLHEP::degree) < 30.*CLHEP::degree) {
              // connection along the intersection
              // keep only the connection with consistent ordering of cells
              experimental_vector c0 = ccb().forward_axis();
              experimental_vector c = c0.hor();
              const double psc2 = std::abs(b.kink_phi(c).value());
              if (psc2 < 60.*CLHEP::degree) {
                DT_LOG_DEBUG(local_priority, "Rejected because first 2 cells intersect and the ordering is wrong: psc = " << psc2);
                continue;
              }
            }
          } else if (intersect_bc) {
            experimental_vector b0 = ccb().forward_axis();
            experimental_vector b = b0.hor();
            const double psc1 = d.kink_phi(b).value();
            if (std::abs(psc1 - 90.*CLHEP::degree) < 30.*CLHEP::degree ||
                std::abs(psc1 + 90.*CLHEP::degree) < 30.*CLHEP::degree ||
                std::abs(psc1 - 270.*CLHEP::degree) < 30.*CLHEP::degree) {
              // connection along the intersection
              // keep only the connection with consistent ordering of cells
              experimental_vector c0 = cca().forward_axis();
              experimental_vector c = c0.hor();
              const double psc2 = std::abs(b.kink_phi(c).value());
              if (psc2 < 60.*CLHEP::degree) {
                DT_LOG_DEBUG(local_priority, "Rejected because last 2 cells intersect and the ordering is wrong: psc = " << psc2);
                continue;
              }
            }
          }

          size_t ndof = 2;  // 2 kink angles, 0 or 1 one separation angle
          if (shall_include_separation) ndof++;

          experimental_double local_separation;
          experimental_point p;
          experimental_double newxa, newza;
          if (cb_.small()) {
            p = cb_.ep();
            newxa = p.x();
            newxa.set_error(cb_.r().error());
            newza = p.z();
            newza.set_error(cb_.r().error());
            p.set_x(newxa);
            p.set_z(newza);
          } else {
            p = cb_.angular_average(i1.epa(), i2.epa(), &local_separation);
          }

          line newt1(i1.epb(), p);
          line newt2(p, i2.epb());
          const experimental_double phi_kink = newt1.kink_phi(newt2);

          if (local_priority >= datatools::logger::PRIO_DEBUG) {
            std::clog << " p1: phi = "; ca_.dump_point(i1.epb()); std::clog << " " << std::endl;
            std::clog << " p2 average point: "; cb_.dump_point(p); std::clog << " " << std::endl;
            std::clog << " p3: "; cc_.dump_point(i2.epb()); std::clog << " " << std::endl;
            std::clog << "    separation: "; (local_separation*180/M_PI).dump(); std::clog << " " << std::endl;
            std::clog << "    phi_kink: " << (phi_kink.value()*180/M_PI) << " " << std::endl;
          }

          bool ok = false;

          auto unknown_vertical = [] (const topology::cell & cell_) -> bool
            {
              if (cell_.ep().y().value() == 0. &&
                  cell_.ep().y().error() > 1000.) return true;
              return false;
            };

          const bool use_theta_kink = !(unknown_vertical(ca_) || unknown_vertical(cb_) || unknown_vertical(cc_));
          if (! use_theta_kink) ndof--;

          double chi2_just_phi = 0.0;
          double chi2 = newt1.chi2(newt2, use_theta_kink, &chi2_just_phi);

          // also useul to check just the phi value (adding the theta information reduces the strength of the phi test)
          const double prob_just_phi = probof(chi2_just_phi,1);

          if (shall_include_separation)
            chi2 += std::pow(local_separation.value()/local_separation.error(), 2);
          chi2s_.push_back(chi2);

          const double local_prob = probof(chi2, ndof);
          probs_.push_back(local_prob);
          if (local_prob > probmin() && prob_just_phi > probmin() && std::abs(phi_kink.value()) <= phi_limit_) {
            ok = true;
          }

          DT_LOG_DEBUG(local_priority, "chi2 " << chi2 << " prob " << local_prob
                       << " prob_just_phi " << prob_just_phi << " phi_kink "
                       << phi_kink.value()/CLHEP::degree << "° limit " << phi_limit_/CLHEP::degree << "° accepted: " << ok);

          if (ok) {
            joint j(newt1.epa(),p,newt2.epb(), get_probmin());
            j.set_chi2(chi2);
            j.set_ndof(ndof);
            j.set_p(probof(chi2, ndof));
            joints_.push_back(j);
          }
        } // end of lines t2
      } // end of lines t1

      joints_ = refine(joints_, ratio_);
      return;
    }

    std::vector<joint> cell_triplet::refine(const std::vector<joint> & joints, double ratio_, size_t max_njoints)
    {
      std::vector<joint> _joints;

      // if( print_level() > mybhep::VERBOSE ){
      //   std::clog << " refining " << joints.size() << " joints " << std::endl;
      // }
      experimental_double delta_phi;
      bool found;
      for(std::vector<joint>::const_iterator ijoint=joints.begin(); ijoint != joints.end(); ++ijoint){
        found = false;
        for(std::vector<joint>::const_iterator jjoint=joints.begin(); jjoint != joints.end(); ++jjoint){
          if( jjoint == ijoint ) continue;

          if( ca_.same_quadrant(ijoint->epa(), jjoint->epa() ) &&
              cc_.same_quadrant(ijoint->epc(), jjoint->epc() ) &&
              ijoint->p() < jjoint->p() ){
            // if( print_level() > mybhep::VERBOSE ){
            //   std::clog << " ... removing joint " << ijoint - joints.begin()  << " with prob " << ijoint->p() << " because joint " << jjoint - joints.begin()  <<  " with prob " << jjoint->p() << " has the same initial and final quadrant " << std::endl;
            // }
            found = true;
            break;
          }
        }

        if( !found )
          _joints.push_back(*ijoint);

      }


      if( _joints.size() > 1 ){
        // order joints in order of increasing chi2 (best joint comes first)
        std::sort(_joints.begin(), _joints.end());
      }

      // only keep best joints
      if( _joints.size() > max_njoints ){
        _joints.erase(_joints.begin() + max_njoints, _joints.end());
      }

      if( _joints.size() >= 2 && !(ca_.intersect(cb_) || ca_.intersect(cc_) || cb_.intersect(cc_) ) ){
        std::vector<joint>::iterator ijoint = _joints.begin();
        ijoint ++;
        while( ijoint != _joints.end() ){
          if( (size_t)(ijoint - _joints.begin() + 1) > _joints.size() ) break;
          if( _joints[0].p() / ijoint->p() > ratio_ ){
            // if( print_level() > mybhep::VERBOSE )
            //   std::clog << " remove joint with p " << ijoint->p() << " in favor of 1st joint with p " << _joints[0].p() << std::endl;
            _joints.erase(ijoint);
            ijoint = _joints.begin() + (ijoint - _joints.begin());
          }else{
            ijoint++;
          }
        }
      }

      // if( print_level() > mybhep::VERBOSE ){
      //   std::clog << " after refining there are " << _joints.size() << " joints " << std::endl;
      //   for(std::vector<joint>::const_iterator ij=_joints.begin(); ij!=_joints.end(); ++ij){
      //     std::clog << " joint " << ij - _joints.begin() << " : "; dump_joint(*ij);
      //   }
      // }

      return _joints;
    }

    size_t cell_triplet::iteration()const{
      for(std::vector<joint>::const_iterator i=joints_.begin(); i!=joints_.end(); ++i)
        if( ! i->used() )
          return (size_t)(i - joints_.begin());

      return joints().size();
    }

    cell_triplet cell_triplet::invert(){
      cell_triplet inverted;
      inverted.set(ccb(),cca());
      inverted.set_free(free());
      inverted.set_begun(begun());
      inverted.set_chi2s(chi2s());
      inverted.set_probs(probs());

      std::vector<joint> inverted_joints;
      for(std::vector<joint>::iterator i=joints_.begin(); i!=joints_.end(); ++i){
        inverted_joints.push_back(i->invert());
      }
      inverted.set_joints( inverted_joints );
      return inverted;
    }


    void cell_triplet::set_all_used(){
      for(std::vector<joint>::iterator i=joints_.begin(); i!=joints_.end(); ++i)
        i->set_used(true);
      set_begun(true);
      return;
    }

    bool operator==(const cell_triplet& left,
                    const cell_triplet& right)
    {

      return ((left.ca().id() == right.ca().id()) && (left.cc().id() == right.cc().id())) ||
        ((left.ca().id() == right.cc().id()) && (left.cc().id() == right.ca().id()));

    }


    bool cell_triplet::same_last_cell(cell c)const{
      return ((this->ca().id() == c.id()) ||
              (this->cc().id() == c.id()) );

    }


  }

}
