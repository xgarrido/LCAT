/* -*- mode: c++ -*- */

#include <iostream>
#include <CAT/printable.h>

namespace CAT {
  namespace topology{

    printable::printable()
    {
    }

    printable::~printable()
    {
    }

    //! print interface
    void printable::info(std::ostream& os) const
    {
      os << " no information available" << std::endl;
    }

    //! print the information using the << operator
    std::ostream& operator<<(std::ostream& os, const printable& v)
    {
      v.info(os);
      return os;
    }

    //! print the information using the << operator
    std::ostream& operator<<(std::ostream& os, const printable* v)
    {
      v->info(os);
      return os;
    }


  }
}
