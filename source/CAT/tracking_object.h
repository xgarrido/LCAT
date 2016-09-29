/* -*- mode: c++ -*- */

#ifndef __CATAlgorithm__tracking_object_h
#define __CATAlgorithm__tracking_object_h 1

namespace CAT{
  namespace topology{


    // a generic tracking object
    class tracking_object //: public printable
    {
    public:

      void set_probmin( double probmin );

      double probmin() const;

      double probof(double chi2, int ndof)const;

      double get_probmin()const;

    private:
      // the minimum probabilty
      double probmin_;
    };

  }
}

#endif // __CATAlgorithm__tracking_object_h
