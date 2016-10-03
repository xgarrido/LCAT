/* -*- mode: c++ -*- */

#include <CAT/node.h>
#include <CAT/utilities.h>

namespace CAT {

  node::node()
  {
    free_ = false;
    is_kink_ = false;
    chi2_ = 0.;
    ndof_=0;
    circle_phi_ = small_neg;
  }

  node::~node()
  {
    return;
  }

  node::node(const cell & c, const std::vector<cell_couplet> & cc, const std::vector<cell_triplet> & ccc)
  {
    c_ = c;
    cc_ = cc;
    ccc_ = ccc;
    free_ = false;
    is_kink_ = false;
    chi2_ = 0.;
    ndof_ = 0;
    setup_cc_maps();
    setup_ccc_maps();
    circle_phi_ = small_neg;
  }

  node::node(const cell &c, double probmin)
  {
    set_probmin(probmin);
    c_ = c;
    free_ = false;
    is_kink_ = false;
    chi2_ = 0.;
    ndof_=0;
    circle_phi_ = small_neg;
  }

  void node::dump(std::ostream & a_out,
                  const std::string & a_title,
                  const std::string & a_indent,
                  bool /*a_inherit*/) const
  {
    std::string indent;
    if (! a_indent.empty ()) indent = a_indent;
    if (! a_title.empty ()) {
      a_out << indent << a_title << std::endl;
    }
    a_out << indent << " --------------------- " << std::endl;
    a_out << indent  << " main cell " << " free " << free() << " chi2 " << chi2() << " ndof " << ndof() << std::endl;
    this->c().tree_dump(a_out,"",indent + "   ");
    a_out << indent << " fitted point: "; ep().dump();
    a_out << indent << " cell couplets: " << cc().size() << std::endl;
    for(std::vector<cell_couplet>::const_iterator icc=cc_.begin(); icc!=cc_.end(); ++icc)
      icc->dump(a_out, "",indent + "     ");
    a_out << indent << " cell triplets: " << ccc().size() << std::endl;
    for(std::vector<cell_triplet>::const_iterator iccc=ccc_.begin(); iccc!=ccc_.end(); ++iccc)
      iccc->dump(a_out, "",indent + "     ");
    a_out << indent  << " --------------------- " << std::endl;
    return;
  }

  void node::setup_cc_maps()
  {
    cc_index_.clear();
    for(std::vector<cell_couplet>::const_iterator icc=cc_.begin(); icc!=cc_.end(); ++icc){
      cc_index_[icc->cb().get_id()] = icc-cc_.begin();
    }
  }

  void node::setup_ccc_maps()
  {
    ccc_ca_index_.clear();
    ccc_cc_index_.clear();
    for(std::vector<cell_triplet>::const_iterator iccc=ccc_.begin(); iccc!=ccc_.end(); ++iccc){
      ccc_ca_index_[iccc->cb().get_id()] = iccc-ccc_.begin();
      ccc_cc_index_[iccc->cc().get_id()] = iccc-ccc_.begin();
    }
  }

  void node::set(const cell & c,
                 const std::vector<cell_couplet> & cc,
                 const std::vector<cell_triplet> & ccc)
  {
    c_ = c;
    cc_ = cc;
    ccc_ = ccc;
    setup_cc_maps();
    setup_ccc_maps();
  }

  void node::set_c(const cell& c)
  {
    c_ = c;
  }

  void node::set_cc(const std::vector<cell_couplet> &cc)
  {
    cc_ = cc;
    setup_cc_maps();
    setup_ccc_maps();
  }

  void node::set_ccc(const std::vector<cell_triplet>  &ccc)
  {
    ccc_ = ccc;
  }

  void node::set_links(const std::vector<cell>  &links)
  {
    links_ = links;
  }

  void node::set_free(bool free)
  {
    free_ = free;
  }

  void node::set_is_kink(bool is_kink)
  {
    is_kink_ = is_kink;
  }

  void node::set_chi2(double chi2)
  {
    chi2_ = chi2;
  }

  void node::set_ndof(int32_t ndof)
  {
    ndof_ = ndof;
  }

  void node::set_ep( const experimental_point &ep )
  {
    ep_ = ep;
  }

  const cell & node::c() const
  {
    return c_;
  }

  const std::vector<cell_couplet> & node::cc() const
  {
    return cc_;
  }

  const std::vector<cell_triplet> & node::ccc() const
  {
    return ccc_;
  }

  const std::vector<cell> & node::links() const
  {
    return links_;
  }

  bool node::free() const
  {
    return free_;
  }

  //! get is_kink
  bool node::is_kink() const
  {
    return is_kink_;
  }

  //! get chi2
  double node::chi2()const
  {
    return chi2_;
  }

  //! get ndof
  int32_t node::ndof()const
  {
    return ndof_;
  }

  //! get prob
  double node::Prob()const
  {
    return probof(chi2(),ndof());
  }

  //! get fitted experimental_point
  const experimental_point& node::ep()const
  {
    return ep_;
  }

  void node::calculate_triplets(double ratio_, double phi_limit_)
  {
    if( cc_.size() < 2 ) return;
    for(std::vector<cell_couplet>::const_iterator icc=cc_.begin(); icc!=cc_.end(); ++icc){
      cell c1 = icc->cb();
      for(std::vector<cell_couplet>::const_iterator jcc=cc_.begin() + (size_t)(icc - cc_.begin()); jcc!=cc_.end(); ++jcc){
        cell c2 = jcc->cb();
        if( c1.get_id() == c2.get_id() ) continue;
        cell_triplet ccc(c1,c_,c2, probmin());
        // if( print_level() >= VVERBOSE ){
        //   std::clog << appname_ << " calculate triplets for three cells: " << ccc.ca().get_id() << "  " << ccc.cb().get_id() << "  " << ccc.cc().get_id() << std::endl;
        // }
        ccc.calculate_joints(ratio_, phi_limit_);
        if( ccc.joints().size() > 0 ){
          // if( print_level() >= VVERBOSE ){
          //   std::clog << appname_ << " adding joints " << std::endl;
          //   for(std::vector<joint>::iterator ijoint = ccc.joints_.begin(); ijoint != ccc.joints_.end(); ++ ijoint )
          //     std::clog << " joint " << ijoint - ccc.joints_.begin() << " phia: " << experimental_vector(ccc.ca().ep(), ijoint->epa()).phi().value()*180./M_PI
          //               << " phib: " << experimental_vector(ccc.cb().ep(), ijoint->epb()).phi().value()*180./M_PI
          //               << " phic: " << experimental_vector(ccc.cc().ep(), ijoint->epc()).phi().value()*180./M_PI << " chi2 " << ijoint->chi2() << std::endl;
          // }
          add_triplet(ccc);
        }
      }
    }
  }

  void node::add_triplet(const cell_triplet &ccc)
  {
    ccc_.push_back(ccc);
    ccc_ca_index_[ccc.ca().get_id()] = ccc_.size() - 1;
    ccc_cc_index_[ccc.cc().get_id()] = ccc_.size() - 1;
    return;
  }

  void node::remove_couplet(size_t index)
  {
    cc_.erase(cc_.begin() + index);
    setup_cc_maps();
    return;
  }

  void node::remove_triplet(size_t index)
  {
    ccc_.erase(ccc_.begin() + index);
    setup_ccc_maps();
    return;
  }

  void node::remove_link(size_t index)
  {
    links_.erase(links_.begin() + index);
    return;
  }

  node node::invert()
  {
    node inverted;
    inverted.set_probmin(probmin());
    inverted.set_c(c());
    inverted.set_cc(cc());
    inverted.set_ccc(ccc());
    inverted.set_free(free());
    inverted.set_chi2(chi2());
    inverted.set_ndof(ndof());
    inverted.set_ep(ep());
    return inverted;
  }


  std::string node::topological_type() const
  {
    if( cc().empty() ) {
      // no cell couplets
      return "ISOLATED";
    }
    if( ccc().empty() ) {
      // no cell triplets
      if( cc().size() == 1 ) {
        return "VERTEX";
      }
      return "MULTI_VERTEX";
    }
    if( ccc().size() == 1 ) {
      // 1 cell triplet
      return "BRIDGE";
    }
    return "OTHER";
  }

  const std::map<size_t,size_t> & node::cc_index()const
  {
    return cc_index_;
  }

  const std::map<size_t,size_t> & node::ccc_ca_index()const
  {
    return ccc_ca_index_;
  }

  const std::map<size_t,size_t> & node::ccc_cc_index()const
  {
    return ccc_cc_index_;
  }

  bool node::has_couplet(const cell & a, cell_couplet * ct)const
  {
    if( !cc_index_.count(a.get_id()) ) return false;
    size_t index=cc_index().find(a.get_id())->second;
    if( index >= cc_.size()) {
      std::clog << " problem: cc index " << index << " for cell a of id " << a.get_id() << " is larger than cc size " << cc_.size() << std::endl;
      dump();
      return false;
    }
    *ct = cc_.at(index);
    return true;
  }

  bool node::has_couplet(const cell& a, size_t* index)const
  {
    if( !cc_index_.count(a.get_id()) ) return false;
    *index= cc_index().find(a.get_id())->second;
    if( *index >= cc_.size() ){
      std::clog << " problem: cc index " << *index << " for cell of id " << a.get_id() << " is larger than cc size " << cc_.size() << std::endl;
      dump();
      return false;
    }
    return true;
  }

  bool node::has_couplet(size_t idd, size_t* index)const
  {
    if( !cc_index_.count(idd) ) return false;
    *index = cc_index().find(idd)->second;
    if( *index >= cc_.size() ){
      std::clog << " problem: cc index " << *index << " for id " << idd << " is larger than cc size " << cc_.size() << std::endl;
      dump();
      return false;
    }
    return true;
  }

  bool node::has_triplet(const cell & a, const cell & c, size_t *index)const
  {
    cell null;
    const std::vector<cell_triplet>::const_iterator & found_triplet
      = std::find(ccc().begin(), ccc().end(), cell_triplet(a, null, c));
    if (found_triplet != ccc().end()) {
      *index = found_triplet - ccc().begin();
      return true;
    }
    return false;
  }

  bool node::has_triplet(const cell &a, const cell &c)const
  {
    cell null;
    if( std::find(ccc().begin(), ccc().end(), cell_triplet(a, null, c) ) != ccc().end() )
      return true;
    return false;
  }

  bool node::has_triplet(const cell &a)const
  {
    for(std::vector<cell_triplet>::const_iterator iccc=ccc_.begin(); iccc!=ccc_.end(); ++iccc){
      int ida = iccc->ca().get_id();
      int idc = iccc->cc().get_id();
      if( ( ida == a.get_id() || idc == a.get_id() ) ){
        return true;
      }
    }
    return false;
  }


  bool operator==(const node& left,
                  const node& right)
  {
    return left.c().get_id() == right.c().get_id();
  }

}
