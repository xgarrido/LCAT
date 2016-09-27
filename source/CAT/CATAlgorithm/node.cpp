/* -*- mode: c++ -*- */

#include <CATAlgorithm/node.h>

namespace CAT {
  namespace topology{

    using namespace std;
    using namespace mybhep;

      //!Default constructor
    node::node()
      {
        appname_= "node: ";
        free_ = false;
        is_kink_ = false;
        chi2_ = 0.;
        ndof_=0;
        circle_phi_ = mybhep::small_neg;
      }

      //!Default destructor
      node::~node()
      {
        return;
      }

      //! constructor
      node::node(const cell & c, const std::vector<cell_couplet> & cc, const std::vector<cell_triplet> & ccc){
        appname_= "node: ";
        c_ = c;
        cc_ = cc;
        ccc_ = ccc;
        free_ = false;
        is_kink_ = false;
        chi2_ = 0.;
        ndof_=0;
        setup_cc_maps();
        setup_ccc_maps();
        circle_phi_ = mybhep::small_neg;
      }

      //! constructor
      node::node(const cell &c, double probmin){
        set_probmin(probmin);
        appname_= "node: ";
        c_ = c;
        free_ = false;
        is_kink_ = false;
        chi2_ = 0.;
        ndof_=0;
        circle_phi_ = mybhep::small_neg;
      }

      /*** dump ***/
      void node::dump (ostream & a_out,
                       const std::string & a_title,
                       const std::string & a_indent,
                       bool /*a_inherit*/) const{
        std::string indent;
        if (! a_indent.empty ()) indent = a_indent;
        if (! a_title.empty ())
          {
            a_out << indent << a_title << std::endl;
          }

        a_out << indent << appname_ << " --------------------- " << std::endl;
        a_out << indent  << " main cell " << " free " << free() << " chi2 " << chi2() << " ndof " << ndof() << std::endl;
        this->c().dump(a_out,"",indent + "   ");
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

    void node::setup_cc_maps(){
      cc_index_.clear();
      for(std::vector<cell_couplet>::const_iterator icc=cc_.begin(); icc!=cc_.end(); ++icc){
        cc_index_[icc->cb().id()] = icc-cc_.begin();
      }
    }

    void node::setup_ccc_maps(){
      ccc_ca_index_.clear();
      ccc_cc_index_.clear();
      for(std::vector<cell_triplet>::const_iterator iccc=ccc_.begin(); iccc!=ccc_.end(); ++iccc){
        ccc_ca_index_[iccc->cb().id()] = iccc-ccc_.begin();
        ccc_cc_index_[iccc->cc().id()] = iccc-ccc_.begin();
      }
    }


      //! set cells
      void node::set(const cell &c, const std::vector<cell_couplet> &cc, const std::vector<cell_triplet> & ccc){
        c_ = c;
        cc_ = cc;
        ccc_ = ccc;
        setup_cc_maps();
        setup_ccc_maps();
      }

      //! set main cell
      void node::set_c(const cell& c){
        c_ = c;
      }

      //! set cell couplets
      void node::set_cc(const std::vector<cell_couplet> &cc){
        cc_ = cc;
        setup_cc_maps();
        setup_ccc_maps();
      }

      //! set cell triplets
      void node::set_ccc(const std::vector<cell_triplet>  &ccc){
        ccc_ = ccc;
      }

      //! set links
      void node::set_links(const std::vector<cell>  &links){
        links_ = links;
      }

      //! set free level
      void node::set_free(bool free){
        free_ = free;
      }

      //! set is_kink
      void node::set_is_kink(bool is_kink){
        is_kink_ = is_kink;
      }

      //! set chi2
      void node::set_chi2(double chi2){
        chi2_ = chi2;
      }

      //! set ndof
      void node::set_ndof(int32_t ndof){
        ndof_ = ndof;
      }

      //! set fitted experimental_point
      void node::set_ep( const experimental_point &ep )
      {
        ep_ = ep;
      }

      //! get main cell
      const cell& node::c()const
      {
        return c_;
      }

      //! get cell couplets
      const std::vector<cell_couplet> &node::cc()const{
        return cc_;
      }

      //! get cell triplets
      const std::vector<cell_triplet> &node::ccc()const{
        return ccc_;
      }

      //! get links
      const std::vector<cell> &node::links()const{
        return links_;
      }

      //! get free level
      bool node::free()const{
        return free_;
      }

      //! get is_kink
      bool node::is_kink()const{
        return is_kink_;
      }

      //! get chi2
      double node::chi2()const{
        return chi2_;
      }

      //! get ndof
      int32_t node::ndof()const{
        return ndof_;
      }

      //! get prob
      double node::Prob()const{
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
          if( c1.id() == c2.id() ) continue;
          cell_triplet ccc(c1,c_,c2, probmin());
          // if( print_level() >= mybhep::VVERBOSE ){
          //   std::clog << appname_ << " calculate triplets for three cells: " << ccc.ca().id() << "  " << ccc.cb().id() << "  " << ccc.cc().id() << std::endl;
          // }
          ccc.calculate_joints(ratio_, phi_limit_);
          if( ccc.joints().size() > 0 ){
            // if( print_level() >= mybhep::VVERBOSE ){
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

    void node::add_triplet(const cell_triplet &ccc){
      ccc_.push_back(ccc);
      ccc_ca_index_[ccc.ca().id()] = ccc_.size() - 1;
      ccc_cc_index_[ccc.cc().id()] = ccc_.size() - 1;
      return;
    }

    void node::remove_couplet(size_t index){
      cc_.erase(cc_.begin() + index);
      setup_cc_maps();
      return;
    }

    void node::remove_triplet(size_t index){
      ccc_.erase(ccc_.begin() + index);
      setup_ccc_maps();
      return;
    }

    void node::remove_link(size_t index){
      links_.erase(links_.begin() + index);
      return;
    }

      node node::invert(){
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

      std::string node::topological_type() const{

        if( cc().empty() ) // no cell couplets
          return "ISOLATED";

        if( ccc().empty() ){ // no cell triplets
          if( cc().size() == 1 )
            return "VERTEX";
          return "MULTI_VERTEX";
        }

        if( ccc().size() == 1 ) // 1 cell triplet
          return "BRIDGE";

        return "OTHER";
      }

    //! grab cc index map
    std::map<size_t,size_t> node::cc_index()const{
      return cc_index_;
    }

    //! grab ccc cb index map
    std::map<size_t,size_t> node::ccc_ca_index()const{
      return ccc_ca_index_;
    }

    //! grab ccc cc index map
    std::map<size_t,size_t> node::ccc_cc_index()const{
      return ccc_cc_index_;
    }


      bool node::has_couplet(const cell & a, cell_couplet* ct)const {

        if( !cc_index_.count(a.id()) ) return false;
        size_t index=cc_index()[a.id()];
        if( index >= cc_.size() ){
          std::clog << " problem: cc index " << index << " for cell a of id " << a.id() << " is larger than cc size " << cc_.size() << endl;
          dump();
          return false;
        }
        *ct = cc_.at(index);
        return true;

      }

      bool node::has_couplet(const cell& a, size_t* index)const{

        if( !cc_index_.count(a.id()) ) return false;
        *index= cc_index()[a.id()];
        if( *index >= cc_.size() ){
          std::clog << " problem: cc index " << *index << " for cell of id " << a.id() << " is larger than cc size " << cc_.size() << endl;
          dump();
          return false;
        }
        return true;
      }

      bool node::has_couplet(size_t idd, size_t* index)const{

        if( !cc_index_.count(idd) ) return false;
        *index= cc_index()[idd];
        if( *index >= cc_.size() ){
          std::clog << " problem: cc index " << *index << " for id " << idd << " is larger than cc size " << cc_.size() << endl;
          dump();
          return false;
        }
        return true;

      }

    bool node::has_triplet(const cell &a, const cell &c, size_t *index)const{

      cell null;
      std::vector<cell_triplet>::const_iterator ftriplet = std::find(ccc().begin(), ccc().end(),cell_triplet(a,null,c) );

      if( ftriplet != ccc().end() ){
        *index = ftriplet - ccc().begin();
        return true;
      }
      return false;
    }

      bool node::has_triplet(const cell &a, const cell &c)const{

        cell null;
        if( std::find(ccc().begin(), ccc().end(),cell_triplet(a,null,c) ) != ccc().end() )
          return true;
        return false;
      }

      bool node::has_triplet(const cell &a)const{
        for(std::vector<cell_triplet>::const_iterator iccc=ccc_.begin(); iccc!=ccc_.end(); ++iccc){
          size_t ida = iccc->ca().id();
          size_t idc = iccc->cc().id();
          if( ( ida == a.id() || idc == a.id() ) ){
            return true;
          }
        }

        return false;

      }


      bool operator==(const node& left,
                      const node& right)
      {

        return left.c().id() == right.c().id();

      }

  }
}
