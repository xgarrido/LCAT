// -*- mode: c++ -*-

#ifndef FALAISE_CAT_PLANE_H
#define FALAISE_CAT_PLANE_H 1

// Standard library:
#include <iostream>
#include <string>

// This project:
#include <CAT/experimental_vector.h>
#include <CAT/experimental_point.h>

namespace CAT {

  /// \brief A plane (calorimeter block) is identified by
  /// a type (Nemo3 or SuperNemo)
  /// one point (center) in middle of block
  /// the plane sizes (sizex, sizey, sizez for SuperNemo
  ///                  for Nemo3 sizex is towards the foil, sizez is tangent to the calo circle)
  /// a normal to the plane (nx, ny, nz) looking towards the foil
  /// a view = ("x" or "y" or "z" for SuperNemo,
  ///           "inner", "outer", "bottom" or "top" for Nemo3)
  /// "face" returns the center of the face towards the foil
  class plane
  {

  private:

    // experimental point
    experimental_point center_;

    // sizes
    experimental_vector sizes_;

    // normal looking towards the foil
    experimental_vector norm_;

    // Nemo3 or SuperNEMO
    std::string type_;

  public:

    //!Default constructor
    plane();

    //!Default destructor
    virtual ~plane();

    //! constructor
    plane(const experimental_point &center,
          const experimental_vector &sizes,
          const experimental_vector &norm);

    /// Smart print
    void dump (std::ostream & a_out         = std::clog,
               const std::string & a_title  = "",
               const std::string & a_indent = "",
               bool a_inherit               = false) const;


    //! set center
    void set_center(const experimental_point &center);

    //! set sizes
    void set_sizes(const experimental_vector &sizes);


    //! set norm
    void set_norm(const experimental_vector &norm);

    //! set type
    void set_type(const std::string & type);

    //! get center
    const experimental_point& center() const;

    //! get sizes
    const experimental_vector& sizes() const;

    // returns the normal looking towards the foil
    const experimental_vector& norm() const;

    // get type
    const std::string& type() const;

    std::string view() const;

    //! get point of the face of the plane
    experimental_point face() const;

    bool intersect(const experimental_point &ep) const;

    bool intersect (const experimental_point &start,
                    const experimental_vector &direction,
                    experimental_point* ep) const;

    // std::vector from the face of the plane to the point
    experimental_vector norm_to_point(const experimental_point &ep) const;

  };

} // namespace CAT

#endif // FALAISE_CAT_PLANE_H
