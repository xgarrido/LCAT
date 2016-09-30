// -*- mode: c++ -*-

#ifndef FALAISE_CAT_TRACKING_OBJECT_H
#define FALAISE_CAT_TRACKING_OBJECT_H 1

namespace CAT {

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

} // namespace CAT

#endif // FALAISE_CAT_TRACKING_OBJECT_H
