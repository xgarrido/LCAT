/* -*- mode: c++ -*- */

#ifndef __CATAlgorithm__printable_h
#define __CATAlgorithm__printable_h 1


#include <iostream>

namespace CAT {
  namespace topology{


    //! Base utility class for conditional print
    /*!
      - derived class implement the info method to print information
      \ingroup base
    */

    class printable {

    public:

      printable();

      virtual ~printable();

      //! print interface
      virtual void info(std::ostream& os = std::clog) const ;

      //! print the information using the << operator
      friend std::ostream& operator<<(std::ostream& os, const printable& v) ;

      //! print the information using the << operator
      friend std::ostream& operator<<(std::ostream& os, const printable* v) ;


    };
  }
}

#endif // __CATAlgorithm__printable_h
