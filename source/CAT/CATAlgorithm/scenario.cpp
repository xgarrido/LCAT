/* -*- mode: c++ -*- */
#include <CATAlgorithm/scenario.h>

namespace CAT {
  namespace topology{

    //!Default constructor
    scenario::scenario()
    {
      appname_= "scenario: ";
      set_probmin(10.);
      //sequences_.clear();
      helix_chi2_ = mybhep::small_neg;
      tangent_chi2_ = mybhep::small_neg;
      ndof_ = mybhep::default_integer;
      n_free_families_ = mybhep::default_integer;
      n_overlaps_ = mybhep::default_integer;
    }

    //!Default destructor
    scenario::~scenario()
    {
    }

    //! constructor
    scenario::scenario(const std::vector<sequence> & seqs, double probmin){
      appname_= "scenario: ";
      set_probmin(probmin);
      sequences_ = seqs;
      helix_chi2_ = mybhep::small_neg;
      tangent_chi2_ = mybhep::small_neg;
      ndof_ = mybhep::default_integer;
      n_free_families_ = mybhep::default_integer;
      n_overlaps_ = mybhep::default_integer;
    }

    /*** dump ***/
    void scenario::dump (std::ostream & a_out      ,
                         const std::string & a_title ,
                         const std::string & a_indent,
                         bool /* a_inherit */          )const{
      {
        std::string indent;
        if (! a_indent.empty ()) indent = a_indent;
        if (! a_title.empty ())
          {
            a_out << indent << a_title << std::endl;
          }

        a_out << indent << appname_ << " -------------- " << std::endl;
        a_out << indent << "helix_chi2 : " << helix_chi2() << "tangent_chi2 : " << tangent_chi2() << " ndof " << ndof() << " helix_prob " << helix_Prob() << " tangent_prob " << tangent_Prob() << std::endl;
        a_out << indent << "n free families : " << n_free_families() << std::endl;
        a_out << indent << "n overlaps : " << n_overlaps() << std::endl;
        for(std::vector<sequence>::const_iterator iseq = sequences_.begin(); iseq != sequences_.end(); ++iseq)
          iseq->dump();
        a_out << indent << " -------------- " << std::endl;

        return;
      }
    }



    //! set experimental_point, radius, error and id;
    void scenario::set(const std::vector<sequence> & seqs){
      appname_= "scenario: ";
      sequences_ = seqs;
      helix_chi2_ = mybhep::small_neg;
      tangent_chi2_ = mybhep::small_neg;
      n_free_families_ = mybhep::default_integer;
      n_overlaps_ = mybhep::default_integer;
    }

    //! set sequences
    void scenario::set_sequences(const std::vector<sequence> & seqs){
      sequences_ = seqs;
    }

    //! set helix_chi2
    void scenario::set_helix_chi2(double helix_chi2){
      helix_chi2_ = helix_chi2;
    }

    //! set tangent_chi2
    void scenario::set_tangent_chi2(double tangent_chi2){
      tangent_chi2_ = tangent_chi2;
    }

    //! set n free families
    void scenario::set_n_free_families(size_t n){
      n_free_families_ = n;
    }

    //! set n overlaps
    void scenario::set_n_overlaps(size_t n){
      n_overlaps_ = n;
    }

    //! set ndof
    void scenario::set_ndof(int32_t n){
      ndof_ = n;
    }

    //! get sequences
    const std::vector<sequence> & scenario::sequences()const
    {
      return sequences_;
    }

    //!get helix_chi2
    double scenario::helix_chi2() const
    {
      return helix_chi2_;
    }

    //!get tangent_chi2
    double scenario::tangent_chi2() const
    {
      return tangent_chi2_;
    }

    //!get ndof
    int32_t scenario::ndof() const
    {
      return ndof_;
    }

    //!get n free families
    size_t scenario::n_free_families() const {return n_free_families_;}

    //!get n overlaps
    size_t scenario::n_overlaps() const {return n_overlaps_;}


    void scenario::calculate_n_overlaps(const std::vector<topology::cell> & cells,
                                        const std::vector<topology::calorimeter_hit> & calos){

      std::vector<int> freecells(cells.size());
      fill(freecells.begin(), freecells.end(), 1);

      std::vector<int> freecalos(calos.size());
      fill(freecalos.begin(), freecalos.end(), 1);

      size_t counter = 0;

      for(std::vector<sequence>::iterator iseq = sequences_.begin(); iseq != sequences_.end(); ++iseq){

        for(std::vector<node>::iterator in = iseq->nodes_.begin(); in != iseq->nodes_.end(); ++in){
          if( in->c().id() >= cells.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: cell " << in->c().id() << " has larger id than n of cells " << cells.size() << std::endl;
            continue;
          }

          if( freecells[in->c().id()] )
            freecells[in->c().id()] = 0;
          else
            counter ++;

        }

        if( iseq->has_decay_helix_vertex() && iseq->decay_helix_vertex_type() == "calo" ){
          if( iseq->calo_helix_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: helix calo " << iseq->calo_helix_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }

          if( freecalos[iseq->calo_helix_id()] )
            freecalos[iseq->calo_helix_id()] = 0;
          else
            counter ++;
        }

        if( iseq->has_helix_vertex() && iseq->helix_vertex_type() == "calo" ){
          if( iseq->helix_vertex_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: helix calo-vertex " << iseq->helix_vertex_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }

          if( iseq->helix_vertex_id() != iseq->calo_helix_id() ){ // avoid double counting if both extrapolations point to the same calo
            if( freecalos[iseq->helix_vertex_id()] )
              freecalos[iseq->helix_vertex_id()] = 0;
            else
              counter ++;
          }
        }

        if( iseq->has_decay_tangent_vertex() && iseq->decay_tangent_vertex_type() == "calo" ){
          if( iseq->calo_tangent_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: tangent calo " << iseq->calo_tangent_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }

          if( iseq->calo_tangent_id() != iseq->calo_helix_id() &&
              iseq->calo_tangent_id() != iseq->helix_vertex_id() ){ // avoid double counting if both extrapolations point to the same calo
            if( freecalos[iseq->calo_tangent_id()] )
              freecalos[iseq->calo_tangent_id()] = 0;
            else
              counter ++;
          }
        }

        if( iseq->has_tangent_vertex() && iseq->tangent_vertex_type() == "calo" ){
          if( iseq->tangent_vertex_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: tangent calo-vertex " << iseq->tangent_vertex_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }

          if( iseq->tangent_vertex_id() != iseq->calo_helix_id() &&
              iseq->tangent_vertex_id() != iseq->calo_tangent_id() &&
              iseq->tangent_vertex_id() != iseq->helix_vertex_id() ){ // avoid double counting if both extrapolations point to the same calo
            if( freecalos[iseq->tangent_vertex_id()] )
              freecalos[iseq->tangent_vertex_id()] = 0;
            else
              counter ++;
          }
        }



      }


      n_overlaps_ = counter;

      return;

    }


    size_t scenario::n_of_common_vertexes(double limit)const{
      double local_distance = 0.;
      double min_local_distance = 0.;
      bool found;
      size_t counter = 0;

      for(std::vector<sequence>::const_iterator iseq = sequences_.begin(); iseq != sequences_.end(); ++iseq){
        min_local_distance = limit;
        found = false;
        for(std::vector<sequence>::const_iterator jseq = iseq; jseq != sequences_.end(); ++jseq)
          {
            if( iseq == jseq ) continue;
            local_distance = 0.;
            if( iseq->common_vertex_on_foil(&(*jseq), &local_distance) ){
              if( local_distance < min_local_distance ){
                min_local_distance = local_distance;
                found = true;
              }
            }
          }
        if( found ) counter ++;
      }
      return counter;
    }


    size_t scenario::n_of_ends_on_wire(void)const{

      size_t counter = 0;

      for(std::vector<sequence>::const_iterator iseq = sequences_.begin(); iseq != sequences_.end(); ++iseq){

        if( !iseq->has_helix_vertex() && !iseq->has_tangent_vertex() )
          counter ++;

        if( !iseq->has_decay_helix_vertex() && !iseq->has_decay_tangent_vertex() )
          counter ++;

      }

      return counter;
    }



    void scenario::calculate_n_free_families(const std::vector<topology::cell> &cells,
                                             const std::vector<topology::calorimeter_hit> & calos){

      std::vector<int> freecells(cells.size());
      fill(freecells.begin(), freecells.end(), 1);

      std::vector<int> freecalos(calos.size());
      fill(freecalos.begin(), freecalos.end(), 1);


      for(std::vector<sequence>::iterator iseq = sequences_.begin(); iseq != sequences_.end(); ++iseq){

        for(std::vector<node>::iterator in = iseq->nodes_.begin(); in != iseq->nodes_.end(); ++in){
          if( in->c().id() >= cells.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: cell " << in->c().id() << " has larger id than n of cells " << cells.size() << std::endl;
            continue;
          }
          else{
            freecells[in->c().id()] = 0;
          }
        }

        if( iseq->has_decay_helix_vertex() && iseq->decay_helix_vertex_type() == "calo" ){

          if( iseq->calo_helix_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: helix calo " << iseq->calo_helix_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }
          freecalos[iseq->calo_helix_id()] = 0;
        }

        if( iseq->has_helix_vertex() && iseq->helix_vertex_type() == "calo" ){

          if( iseq->helix_vertex_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: helix calo-vertex " << iseq->helix_vertex_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }
          freecalos[iseq->helix_vertex_id()] = 0;
        }

        if( iseq->has_decay_tangent_vertex() && iseq->decay_tangent_vertex_type() == "calo" ){

          if( iseq->calo_tangent_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: tangent calo " << iseq->calo_tangent_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }
          freecalos[iseq->calo_tangent_id()] = 0;
        }

        if( iseq->has_tangent_vertex() && iseq->tangent_vertex_type() == "calo" ){

          if( iseq->tangent_vertex_id() >= calos.size() ){
            // if( print_level() >= mybhep::VVERBOSE )
            //   std::clog << " problem: tangent calo-vertex " << iseq->tangent_vertex_id() << " has larger id than n of calos " << calos.size() << std::endl;
            continue;
          }
          freecalos[iseq->tangent_vertex_id()] = 0;
        }



      }

      size_t counter = 0;
      for(std::vector<int>::iterator i=freecells.begin(); i!= freecells.end(); ++i)
        if( *i )
          counter ++;

      for(std::vector<int>::iterator i=freecalos.begin(); i!= freecalos.end(); ++i)
        if( *i )
          counter ++;

      n_free_families_ = counter;

      return;
    }


    void scenario::calculate_chi2(){

      double helix_chi2 = 0.;
      double tangent_chi2 = 0.;
      int32_t ndof = 0;
      for(std::vector<sequence>::iterator iseq = sequences_.begin(); iseq != sequences_.end(); ++iseq){
        helix_chi2 += iseq->helix_chi2();
        tangent_chi2 += iseq->chi2();
        ndof += iseq->ndof();
      }

      helix_chi2_ = helix_chi2;
      tangent_chi2_ = tangent_chi2;
      ndof_ = ndof;

      return;
    }


    double scenario::helix_Prob()const{
      return probof(helix_chi2(), ndof());
    }

    double scenario::tangent_Prob()const{
      return probof(helix_chi2(), ndof());
    }

    bool scenario::better_scenario_than( const scenario & s, double limit)const{

      // - n of recovered cells
      int deltanfree = n_free_families() - s.n_free_families();

      // n of new overlaps
      int deltanoverls = n_overlaps() - s.n_overlaps();

      double deltaprob_helix = helix_Prob() - s.helix_Prob();
      double deltachi_helix = helix_chi2() - s.helix_chi2();
      double deltaprob_tangent = tangent_Prob() - s.tangent_Prob();
      double deltachi_tangent = tangent_chi2() - s.tangent_chi2();

      // if( print_level() >= mybhep::VVERBOSE ){
      //   std::clog << " delta n_free_families = (" << n_free_families()  << " - " << s.n_free_families() << ")= " << deltanfree
      //             << " dela n_overlaps = (" << n_overlaps() << " - " << s.n_overlaps() << ")= " << deltanoverls
      //   	  << " delta prob_helix = (" << helix_Prob()  << " - " << s.helix_Prob() << ") = " << deltaprob_helix
      //   	  << " delta prob_tangent = (" << tangent_Prob()  << " - " << s.tangent_Prob() << ") = " << deltaprob_tangent
      //             << std::endl;
      // }

      if( deltanoverls < - 2*deltanfree )
        return true;

      if( deltanoverls == - 2*deltanfree ){

        int delta_n_common_vertexes = n_of_common_vertexes(limit) - s.n_of_common_vertexes(limit);
        // if( print_level() >= mybhep::VVERBOSE )
        //   std::clog << " delta n common vertex = " << delta_n_common_vertexes << std::endl;
        if( delta_n_common_vertexes > 0 ) return true;
        if( delta_n_common_vertexes < 0 ) return false;

        int delta_n_of_ends_on_wire = n_of_ends_on_wire() - s.n_of_ends_on_wire();
        // if( print_level() >= mybhep::VVERBOSE )
        //   std::clog << " delta n ends on wire = " << delta_n_of_ends_on_wire << std::endl;
        if( delta_n_of_ends_on_wire < 0 ) return true;
        if( delta_n_of_ends_on_wire > 0 ) return false;

        if( deltaprob_helix > 0. )
          return true;

        if( deltaprob_tangent > 0. )
          return true;

        if( deltaprob_helix == 0. && deltachi_helix < 0. )
          return true;

        if( deltaprob_tangent == 0. && deltachi_tangent < 0. )
          return true;

      }

      return false;
    }

  }
}
