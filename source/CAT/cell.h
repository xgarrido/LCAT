/* -*- mode: c++ -*- */
#ifndef FALAISE_CAT_CELL_H
#define FALAISE_CAT_CELL_H 1

// Standard library
#include <iostream>

// This project
#include <CAT/experimental_point.h>
#include <CAT/experimental_vector.h>

namespace CAT {

  /// \brief Cell class
  // A cell is composed of an experimental point and an experimental radius
  class cell {

  public:

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

    /// Set Geiger cell position
    void set_position(const experimental_point &);

    /// Set Geiger cell radius
    void set_radius(double);

    /// Set Geiger cell radius error
    void set_radius_error(double);

    /// Set small radius
    void set_small_radius(double);

    /// Set Geiger cell id
    void set_id(int);

    /// Set Geiger cell side
    void set_side(int);

    /// Set Geiger cell row
    void set_row(int);

    /// Set Geiger cell layer
    void set_layer(int);

    /// Set Geiger as prompt
    void set_prompt(bool);

    //! set free level
    void set_free(bool free);

    //! set begun level
    void set_begun(bool begun);

    //! get experimental_point
    const experimental_point & get_position() const;

    //!get experimental r
    const experimental_double & get_radius() const;

    //!get original experimental r
    const experimental_double & get_original_radius() const;

    //!get id
    int get_id() const;

    //!get layer
    int get_layer() const;

    //!get side
    int get_side() const;

    //!get row
    int get_row() const;

    //!get fast flag
    bool is_prompt() const;

    //! get small flag
    bool is_small() const;

    //! get free level
    bool is_free() const;

    //! get begun level
    bool begun() const;

  public:

    experimental_point angular_average(const experimental_point & epa_, const experimental_point & epb_, experimental_double & angle_);

    experimental_point build_from_cell(const experimental_vector & forward_,
                                       const experimental_vector & transverse_,
                                       const experimental_double & cos_,
                                       int sign_, bool replace_r_, double max_r_) const;

    bool same_quadrant(const experimental_point & epa_, const experimental_point & epb_) const;

    bool intersect(const cell & c_) const;


  private:
    void _set_radius_();

  private:

    // experimental point
    experimental_point _position_;

    // radius (original value)
    experimental_double _r0_;

    // radius (modified value if the cell is "small")
    experimental_double _r_;

    // id
    int _id_;

    // layer number
    int _layer_;

    // side number
    int _side_;

    // row number
    int _row_;

    // characterize fast and delayed cells
    bool _prompt_;

    // radius below which a cell is small
    double _small_radius_;

    // status of cell couplet
    bool _free_;

    // begun cell couplet
    bool _begun_;

  };
}

#endif // FALAISE_CAT_CELL_BASE_H
