/* -*- mode: c++ -*- */
#ifndef __CATAlgorithm__ICELL
#define __CATAlgorithm__ICELL

#include <iostream>
#include <cmath>
#include <CAT/utilities.h>
#include <CAT/tracking_object.h>
#include <CAT/experimental_point.h>
#include <CAT/experimental_vector.h>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/clhep_units.h>

namespace CAT {

  namespace topology {


    /// \brief Cell class
    // A cell is composed of an experimental point and an experimental radius
    class cell : public tracking_object {

    public:

      // status of cell couplet
      bool free_;

      // begun cell couplet
      bool begun_;

      /// Default constructor
      cell();

      /// Destructor
      virtual ~cell();

      /// Reset
      void reset();

      /// Smart dump
      void tree_dump(std::ostream & out_         = std::clog,
                     const std::string & title_  = "",
                     const std::string & indent_ = "",
                     bool inherit_               = false) const;

      //! set experimental_point
      void set_p(const experimental_point & p)
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
      void set_free(bool free)
      {
        free_ = free;
      }

      //! set begun level
      void set_begun(bool begun)
      {
        begun_ = begun;
      }

      bool small() const
      {
        if (r0_.value() <= small_radius_) return true;
        return false;
      }

      //! get experimental_point
      const experimental_point & ep() const
      {
        return ep_;
      }

      //!get experimental r
      const experimental_double & r() const
      {
        return r_;
      }

      //!get original experimental r
      const experimental_double & r0() const
      {
        return r0_;
      }

      //!get id
      size_t id() const {return id_;}

      // //!get small_radius
      //  double small_radius() const {return small_radius_;}

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

      //! get cell number
      int cell_number() const{
        return iid();
      }


    public:

      experimental_point angular_average(const experimental_point & epa_, const experimental_point & epb_, experimental_double & angle_);

      experimental_point build_from_cell(const experimental_vector & forward_,
                                         const experimental_vector & transverse_,
                                         const experimental_double & cos_,
                                         int sign_, bool replace_r_, double max_r_) const;

      bool same_quadrant(const experimental_point & epa_, const experimental_point & epb_) const;

      bool intersect(const topology::cell & c_) const;


    private:
      void set_radius(){
        r_ = r0_;
        /*
          if( small() && fast() )
          r_.set_error(std::max(r0_.value(), r0_.error()));
        */
      }

      private:

      // experimental point
      experimental_point ep_;

      // radius (original value)
      experimental_double r0_;

      // radius (modified value if the cell is "small")
      experimental_double r_;

      // id
      size_t id_;

      // characterize fast and delayed cells
      bool fast_;

      // layer number
      int layer_;

      // block number
      int block_;

      // iid number
      int iid_;

      // radius below which a cell is small
      double small_radius_;

    };
  }
}

#endif
