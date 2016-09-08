/* -*- mode: c++ -*- */
#ifndef __CATAlgorithm__ICELL
#define __CATAlgorithm__ICELL
#include <iostream>
#include <cmath>


#include <CAT/experimental_double.h>
#include <CAT/experimental_point.h>

namespace CAT {

  class cell {

    // a cell is composed of an experimental point
    // and an experimental radius

  private:
    // the minimum probabilty
    double probmin_;

    std::string appname_;

    // experimental point
    experimental_point ep_;

    // radius (original value)
    experimental_double r0_;

    // radius (modified value if the cell is "small")
    experimental_double r_;

    // id
    size_t id_;
    int user_id_; // use this value to store user id, never inside tracking

    // characterize fast and delayed cells
    bool fast_;

    // layer number
    int layer_;

    // block number
    int block_;

    // iid number
    int iid_;

    // N3 or SN
    std::string type_;

    // radius below which a cell is small
    double small_radius_;

  public:

    void set_probmin( double probmin );

    double probmin() const;

    double probof(double chi2, int ndof)const;

    double get_probmin()const;

    // status of cell couplet
    bool free_;

    // begun cell couplet
    bool begun_;

    //!Default constructor
    cell()
    {
      set_probmin(10.);
      //ep_ = experimental_point();
      //r0_= experimental_double();
      //r_= experimental_double();
      // id_ = mybhep::default_integer;
      // user_id_ = mybhep::default_integer;
      // layer_ = mybhep::default_integer;
      // block_ = mybhep::default_integer;
      // iid_ = mybhep::default_integer;
      // n3id_ = mybhep::default_integer;
      fast_ = true;
      free_ = false;
      begun_ = false;
      type_ ="SN";
      small_radius_= 0.;
    }

    //!Default destructor
    virtual ~cell(){};

    /*** dump ***/
    // virtual void dump (std::ostream & a_out         = std::clog,
    //                    const std::string & a_title  = "",
    //                    const std::string & a_indent = "",
    //                    bool /*a_inherit */         = false) const{
    //   {
    //     std::string indent;
    //     if (! a_indent.empty ()) indent = a_indent;
    //     if (! a_title.empty ())
    //       {
    //         a_out << indent << a_title << std::endl;
    //       }

    //     a_out << indent << appname_ << " -------------- " << std::endl;
    //     a_out << indent << "id : " << this->id() << " layer " << this->layer() << " block " << this->block() << " iid " << this->iid() << " n3id " << this->n3id() << " fast " << this->fast() << " small " << this->small() << " unknown vertical " << this->unknown_vertical() << std::endl;
    //     a_out << indent << " point " << std::endl;
    //     this->ep().dump(a_out,"", indent + "   ");
    //     a_out << indent << "radius : "; (r()/mybhep::mm).dump(); a_out << " [mm ] " << std::endl;
    //     if( small() && fast() ){
    //       a_out << indent << "original radius : "; (r0()/mybhep::mm).dump(); a_out << " [mm ] " << std::endl;
    //     }
    //     a_out << indent << " -------------- " << std::endl;

    //     return;
    //   }
    // }



    // //! set experimental_point, radius, error and id;
    // void set(experimental_point p, double r,double er, size_t id, bool fast)
    // {
    //   ep_ = p;
    //   r0_.set_value(r);
    //   r0_.set_error(er);
    //   id_ = id;
    //   fast_ = fast;
    //   set_radius();
    // }

    //! set experimental_point
    void set_position(const experimental_point & p)
    {
      ep_ = p;
    }

    //! set radius
    void set_r(double r)
    {
      r0_.set_value(r);
      set_radius();
    }

    //! set radius error
    void set_er(double er)
    {
      r0_.set_error(er);
      set_radius();
    }

    //! set id
    void set_id(size_t id)
    {
      id_ = id;
    }
    void set_user_id(size_t user_id)
    {
      user_id_ = user_id;
    }

    //! set small radius
    void set_small_radius(double small_radius)
    {
      small_radius_ = small_radius;
    }


    //! set layer
    void set_layer(size_t layer)
    {
      layer_ = layer;
    }

    //! set block
    void set_block(int block)
    {
      block_ = block;
    }

    //! set iid
    void set_iid(size_t iid)
    {
      iid_ = iid;
    }

    //! set fast flag
    void set_fast(bool fast)
    {
      fast_ = fast;
    }


    //! set free level
    void set_free(bool free){
      free_ = free;
    }

    //! set begun level
    void set_begun(bool begun){
      begun_ = begun;
    }

    //! set type
    void set_type(std::string type){
      type_ = type;
    }


    bool small() const
    {
      bool sm = false;

      //        if (r0_.value() <= probmin()*r0_.error() ) sm = true;
      if (r0_.value() <= small_radius() ) sm = true;

      return sm;
    }

    bool unknown_vertical() const
    {
      bool uv = false;

      if( ep().y().value() == 0. &&
          ep().y().error() > 1000. ) uv = true;

      return uv;
    }

    //! get experimental_point
    const experimental_point& ep()const
    {
      return ep_;
    }

    //!get experimental r
    const experimental_double& r() const
    {
      return r_;
    }

    //!get original experimental r
    const experimental_double& r0() const
    {
      return r0_;
    }

    //!get id
    size_t id() const {return id_;}
    const int& user_id() const {return user_id_;}

    //!get small_radius
    double small_radius() const {return small_radius_;}

    //!get layer
    int layer() const {return layer_;}

    //!get block
    int block() const {return block_;}

    //!get iid
    int iid() const {return iid_;}

    //!get fast flag
    bool fast() const {return fast_;}

    //! get free level
    bool free()const{
      return free_;
    }

    //! get begun level
    bool begun()const{
      return begun_;
    }

    //! get type
    std::string type(){
      return type_;
    }

    //! get cell number
    int cell_number() const{

      return iid();
    }


  public:

    // experimental_double distance(cell c) const;
    // experimental_point angular_average(experimental_point epa, experimental_point epb, experimental_double* angle);
    // experimental_point build_from_cell(experimental_vector forward, experimental_vector transverse, experimental_double cos, int sign, bool replace_r, double maxr) const;
    // void dump_point(experimental_point ep) const;
    // void dump_point_phi(experimental_point ep) const;
    // bool same_quadrant(experimental_point epa, experimental_point epb) const;
    // bool same_cell(topology::cell c) const;
    // bool intersect(topology::cell c) const;



    bool operator<(const CAT::cell& c) const{

      // if( id_ > mybhep::default_integer || c.id() > mybhep::default_integer ){
      //   if( print_level() >= mybhep::NORMAL )
      //     std::clog << " problem: trying to compare cells with ids " << id_ << " and " << c.id() << " just returning false " << std::endl;
      //   return false;
      // }

      if( this->id() == c.id() ) return false;

      // side of foil
      if( this->block() < 0 && c.block() > 0 ){
        return false;
      }
      if( this->block() > 0 && c.block() < 0 )
        return true;


      // layer
      if(std::abs(this->layer()) < std::abs(c.layer())){
        return false;
      }
      if(std::abs(this->layer()) > std::abs(c.layer()))
        return true;

      // iid
      if(this->iid() < c.iid()){
        return false;
      }
      if(this->iid() > c.iid()){
        return true;
      }


      return true;

    }

    static bool compare(const CAT::cell& c1, const CAT::cell& c) {

      // side of foil
      if( c1.block() < 0 && c.block() > 0 ){
        return false;
      }

      // layer
      if(std::abs(c1.layer()) < std::abs(c.layer())){
        return false;
      }
      if(std::abs(c1.layer()) > std::abs(c.layer()))
        return true;

      // iid
      if(c1.iid() < c.iid()){
        return false;
      }
      if(c1.iid() > c.iid())
        return true;



      return true;

    };

  private:
    void set_radius(){
      r_ = r0_;
      /*
        if( small() && fast() )
        r_.set_error(std::max(r0_.value(), r0_.error()));
      */
    }


  };
}

#endif
