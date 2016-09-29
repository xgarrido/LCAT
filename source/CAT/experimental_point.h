// -*- mode: c++ -*-

#ifndef CAT_TOPOLOGY_EXPERIMENTAL_POINT_H
#define CAT_TOPOLOGY_EXPERIMENTAL_POINT_H

// Standard library:
#include <iostream>
#include <cmath>

// This project:
#include <CAT/experimental_double.h>

namespace CAT {

  namespace topology {

    /// \brief An experimental point is composed of an ordinary point (x, y, z)
    /// with corresponding error (ex, ey, ez)
    class experimental_point
    {

    private:

      // x coordinate
      experimental_double x_;

      // y coordinate
      experimental_double y_;

      // z coordinate
      experimental_double z_;

      // radius in horizontal view (= x^2 + z^2)
      experimental_double radius_;

    public:
      //!Default constructor
      experimental_point();

      //!Default destructor
      virtual ~experimental_point();

      //! constructor
      experimental_point(const experimental_double &x, const experimental_double &y, const experimental_double &z);

      //! constructor from coordinates with error
      experimental_point(double x, double y, double z, double ex, double ey, double ez);

      /*** dump ***/
      virtual void dump (std::ostream & a_out         = std::clog,
                         const std::string & a_title  = "",
                         const std::string & a_indent = "",
                         bool a_inherit               = false) const;

      //! set x
      void set_x(const experimental_double &x);

      //! set ex
      void set_ex(double ex);

      //! set y
      void set_y(const experimental_double &y);

      //! set ey
      void set_ey(double ey);

      //! set z
      void set_z(const experimental_double &z);

      //! set ez
      void set_ez(double ez);

      //!get experimental x
      const experimental_double& x() const;

      //!get experimental y
      const experimental_double& y() const;

      //!get experimental z
      const experimental_double& z() const;

      //!get experimental radius
      const experimental_double& radius() const;

      //! distance
      experimental_double distance(const experimental_point& p2) const;

      //! distance
      experimental_double hor_distance(const experimental_point& p2) const;

    private:

      void set_radius();
    };

    // average
    experimental_point average(const std::vector<experimental_point> &vs);

  } // namespace topology

} // namespace CAT

#endif // CAT_TOPOLOGY_EXPERIMENTAL_POINT_H
