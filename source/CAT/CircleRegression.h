/* -*- mode: c++ -*- */
#ifndef FALAISE_CAT_CIRCLEREGRESSION_H
#define FALAISE_CAT_CIRCLEREGRESSION_H 1

#include <iostream>
#include <cmath>
#include <CAT/experimental_double.h>
#include <CAT/circle_base.h>


namespace CAT {

  // Operate weighted circle regression on points (xi, yi) to find center and
  // radius of best fitting circle
  class CircleRegression
  {

    private:

      circle c_;

      std::vector<experimental_double> yi_;
      std::vector<experimental_double> xi_;

    public:

      //!Default constructor
      CircleRegression()
      {
        c_ = circle();
        xi_.clear();
        yi_.clear();
      }


      //!Default destructor
      virtual ~CircleRegression(){}

      //! constructor
      CircleRegression(std::vector<experimental_double> &xi, std::vector<experimental_double> &yi){
        xi_ = xi;
        yi_ = yi;
      }



      /*** dump ***/
      virtual void dump (std::ostream & a_out         = std::clog,
                         const std::string & a_title  = "",
                         const std::string & a_indent = "",
                         bool /*a_inherit*/          = false){
        {
          std::string indent;
          if (! a_indent.empty ()) indent = a_indent;
          if (! a_title.empty ())
            {
              a_out << indent << a_title << std::endl;
            }

          a_out << indent << " points: " << std::endl;
          experimental_double phi(0.,0.);
          double phi_ref = 0.;
          for(std::vector<experimental_double>::iterator it=xi_.begin(); it != xi_.end(); ++it){
            experimental_double y = yi_[it - xi_.begin()];
            phi_ref = phi.value();
            phi = c_.phi_of_point(experimental_point(*it, experimental_double(0.,0.), y),phi_ref);
            a_out << indent << " .. x "; it->dump(); a_out << " y "; yi_[it - xi_.begin()].dump(); a_out << " predicted x "; position(phi).x().dump(); a_out << " y "; position(phi).z().dump(); a_out << " phi "; phi.dump(); a_out << " "  << std::endl;
          }
          a_out << indent << " circle: " << std::endl;
          circle().dump();
          a_out << indent << " -------------- " << std::endl;

          return;
        }
      }


      //! set
      void set(std::vector<experimental_double> &xi, std::vector<experimental_double> &yi);


      //! set xi
      void set_xi(std::vector<experimental_double> &xi);

      //! set yi
      void set_yi(std::vector<experimental_double> &yi);

      //! set circle
      void set_circle(circle &c);

      //! get xi
      const std::vector<experimental_double>& xi()const;

      //! get yi
      const std::vector<experimental_double>& yi()const;

      //! get circle
      const circle& c()const;


      bool fit(void);

      experimental_point position(experimental_double &phi);

      bool points_in_good_order();


    };

}

#endif
