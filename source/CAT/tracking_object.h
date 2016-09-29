// -*- mode: c++ -*-

#ifndef CAT_TOPOLOGY_TRACKING_OBJECT_H
#define CAT_TOPOLOGY_TRACKING_OBJECT_H

namespace CAT {

  namespace topology {

    /// \brief A generic tracking object
    class tracking_object
    {
    public:

      void set_probmin(double probmin);

      double probmin() const;

      double probof(double chi2, int ndof) const;

      double get_probmin() const;

    private:

      double probmin_; ///< the minimum probability

    };

  } // namespace topology

} // namespace CAT

#endif // CAT_TOPOLOGY_TRACKING_OBJECT_H
