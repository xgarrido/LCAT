/* -*- mode: c++ -*- */
#include <CATAlgorithm/CircleRegression.h>
#if CAT_WITH_DEVEL_ROOT == 1
#include <TROOT.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#endif // CAT_WITH_DEVEL_ROOT == 1

namespace CAT{
  namespace topology{

#if CAT_WITH_DEVEL_ROOT == 1

    double CircleRegression::expression_to_be_minimized(const double *par)
      {

        const Double_t xc =par[0];
        const Double_t yc =par[1];
        const Double_t rad =par[2];


        double expression = 0.;
        for (size_t i=0;i<xi_.size(); i++) {
          expression += std::abs(pow(xi_[i].value() - xc,2) + pow(yi_[i].value() - yc,2) - pow(rad,2));
          /*
          std::clog << " point " << i << " x " << xi_[i].value() << " par " << xc << " diff " << xi_[i].value() - xc
                    << " y " << yi_[i].value() << " par " << yc << " diff " << yi_[i].value() - yc
                    << " r " <<  rad << " exp " << pow(xi_[i].value() - xc,2) + pow(yi_[i].value() - yc,2) - pow(rad,2) << " cumul " << expression << std::endl;
          */
        }

        //      std::clog << " xc " << par[0] << " yc " << par[1] << " r " << par[2] << " val " << expression << std::endl;

        return expression;
      }
#endif // CAT_WITH_DEVEL_ROOT == 1


      //! set
      void CircleRegression::set(std::vector<experimental_double> &xi, std::vector<experimental_double> &yi)
      {
        xi_ = xi;
        yi_ = yi;
      }


      //! set xi
      void CircleRegression::set_xi(std::vector<experimental_double> &xi)
      {
        xi_ = xi;
      }

      //! set yi
      void CircleRegression::set_yi(std::vector<experimental_double> &yi)
      {
        yi_ = yi;
      }

      //! set circle
      void CircleRegression::set_circle(circle &c)
      {
        c_ = c;
      }

      //! get xi
      const std::vector<experimental_double>& CircleRegression::xi()const
      {
        return xi_;
      }

      //! get yi
      const std::vector<experimental_double>& CircleRegression::yi()const
      {
        return yi_;
      }

      //! get circle
      const circle& CircleRegression::c()const
      {
        return c_;
      }



      bool CircleRegression::fit(void){


        if( xi_.size() != yi_.size() ){
          // if( print_level() >= mybhep::NORMAL ){
          //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, sizes x " << xi_.size() << " y " << yi_.size() << std::endl;
          // }
          return false;
        }

        experimental_double xc, yc, r;


        bool method1 = false;
        // method 1: R. Bullock, http://www.dtcenter.org/met/users/docs/write_ups/circle_fit.pdf
        // method 2: D. Umbach, K. N. Jones, http://www.cs.bsu.edu/homepages/kerryj/kjones/circles.pdf

        std::vector<double> ui;
        std::vector<double> vi;
        double Sw = 0.;
        double Swx = 0.;
        double Swy = 0.;
        double Swxx = 0.;
        double Swxy = 0.;
        double Swyy = 0.;
        double Swxxx = 0.;
        double Swxyy = 0.;
        double Swxxy = 0.;
        double Swyyy = 0.;
        double Swuu = 0.;
        double Swuv = 0.;
        double Swvv = 0.;
        double Swuuu = 0.;
        double Swuuv = 0.;
        double Swuvv = 0.;
        double Swvvv = 0.;

        double delta, xave, yave;

        if( method1 ){
          xave = average(xi_).value();
          yave = average(yi_).value();
          std::clog << "CAT::CircleRegression::fit: xave " << xave << " yave " << yave << std::endl;

        }

	// calculate parameters for circle
        for(std::vector<experimental_double>::iterator it=xi_.begin(); it != xi_.end(); ++it)
          {
            double y = yi_[it - xi_.begin()].value();
            double yerr = yi_[it - xi_.begin()].error();
            double w = 1./(std::pow(it->error(),2)) + 1./(std::pow(yerr,2));
            if( std::isnan(w) || std::isinf(w) )
              w = 1.;
            Sw += w;

            if( method1 ){
              double u = it->value() - xave;
              double v = y - yave;
              Swuu += w*u*u;
              Swuv += w*u*v;
              Swvv += w*v*v;
              Swuuu += w*u*u*u;
              Swuuv += w*u*u*v;
              Swuvv += w*u*v*v;
              Swvvv += w*v*v*v;
            } else {
              double x = it->value();
              Swx += w*x;
              Swy += w*y;
              Swxx += w*x*x;
              Swxy += w*x*y;
              Swyy += w*y*y;
              Swxyy += w*x*y*y;
              Swxxy += w*x*x*y;
              Swyyy += w*y*y*y;
              Swxxx += w*x*x*x;
            }
          }

        if( method1 ){
          delta = Swuu*Swvv - Swuv*Swuv;

          if( delta == 0.){
            // if( print_level() >= mybhep::NORMAL ){
            //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, delta " << delta << " Swuu " << Swuu << " Swvv " << Swvv << " Swuv " << Swuv << std::endl;
            // }
            return false;
          }

          double uc = (Swuuu + Swuvv)/(2.*delta);
          double vc = (Swuuv + Swvvv)/(2.*delta);
          double erruc = 0.;
          double errvc = 0.;
          double alpha = uc*uc + vc*vc + (Swuu + Swvv)/Sw;
          double erralpha = 0.;

          std::clog << "CAT::CircleRegression::fit: uc " << uc << " vc " << vc << std::endl;

          xc.set(uc + xave, erruc);
          yc.set(vc + yave, errvc);
          r.set(std::sqrt(alpha), erralpha/(2.*std::sqrt(alpha)));
        }
        else{
          double A = Sw*Swxx - Swx*Swx;
          double B = Sw*Swxy - Swx*Swy;
          double C = Sw*Swyy - Swy*Swy;
          double D = (Sw*Swxyy - Swx*Swyy + Sw*Swxxx - Swx*Swxx)/2.;
          double E = (Sw*Swxxy - Swy*Swxx + Sw*Swyyy - Swy*Swyy)/2.;
          delta = A*C - B*B;

          if( std::isnan(A) || std::isinf(A) ){
            // if( print_level() >= mybhep::NORMAL ){
            //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, A " << A << " Sw " << Sw << " Swxx " << Swxx << " Swx " << Swx << std::endl;
            // }
            return false;
          }
          if( std::isnan(B) || std::isinf(B) ){
            // if( print_level() >= mybhep::NORMAL ){
            //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, B " << B << " Sw " << Sw << " Swxy " << Swxy << " Swx " << Swx << " Swy " << Swy << std::endl;
            // }
            return false;

          }
          if( std::isnan(C) || std::isinf(C) ){
            // if( print_level() >= mybhep::NORMAL ){
            //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, C " << C << " Sw " << Sw << " Swyy " << Swyy << " Swy " << Swy << std::endl;
            // }
            return false;

          }
          if( std::isnan(D) || std::isinf(D) ){
            // if( print_level() >= mybhep::NORMAL ){
            //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, D " << D << " Sw " << Sw << " Swxyy " << Swxyy << " Swx " << Swx << " Swyy " << Swyy << " Swxxx " << Swxxx << " Swxx " << Swxx <<std::endl;
            // }
            return false;

          }
          if( std::isnan(E) || std::isinf(E) ){
            // if( print_level() >= mybhep::NORMAL ){
            //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, E " << E << " Sw " << Sw << " Swxxy " << Swxxy << " Swy " << Swy << " Swxx " << Swxx << " Swyyy " << Swyyy << " Swyy " << Swyy <<std::endl;
            // }
            return false;

          }


          if( delta == 0.){
            // if( print_level() >= mybhep::NORMAL ){
            //   std::clog << "CAT::CircleRegression::fit: problem: in circle regression, delta " << delta << " A " << A << " C " << C << " B " << B << std::endl;
            // }
            return false;
          }

	  double XC = (D*C - B*E)/delta;
	  double YC = (A*E - B*D)/delta;
          double rsum = 0.;
          for(std::vector<experimental_double>::iterator it=xi_.begin(); it != xi_.end(); ++it)
            {
              double u = it->value() - XC;
              double y = yi_[it - xi_.begin()].value();
              double v = y - YC;
              rsum += std::hypot(u, v);
            }
	  rsum /= xi_.size();

          xc.set(XC, 0.);
          yc.set(YC, 0.);
          r.set(rsum , 0. );
	  c_ = circle(experimental_point(xc, experimental_double(0.,0.), yc), r, probmin());

	  double mean_error_x = 0.;
	  double mean_error_y = 0.;
	  double mean_error_r = 0.;
	  for(std::vector<experimental_double>::iterator it=xi_.begin(); it != xi_.end(); ++it)
	    {
	      topology::experimental_point local_point(*it, experimental_double(0.,0.), yi_[it - xi_.begin()]);
	      topology::experimental_point circle_point = c_.position(local_point);
	      topology::experimental_vector distance(local_point, circle_point);
	      mean_error_x += pow(distance.x().value(),2);
	      mean_error_y += pow(distance.z().value(),2);
	      mean_error_r += pow(distance.x().value(),2) + pow(distance.z().value(),2);
	    }

          xc.set_error(std::sqrt(mean_error_x/xi_.size()));
          yc.set_error(std::sqrt(mean_error_y/xi_.size()));
	  r.set_error(std::sqrt(mean_error_r/xi_.size()));

          // if( print_level() >= mybhep::VVERBOSE ){
          //   std::clog << "CAT::CircleRegression::fit: fitted circle through " << xi_.size() << " points: xc: "; xc.dump(); std::clog << " yc: "; yc.dump(); std::clog << " r: "; r.dump(); std::clog << " " << std::endl;
          // }
        }

        c_ = circle(experimental_point(xc, experimental_double(0.,0.), yc), r, probmin());

        return points_in_good_order();

      }

    experimental_point CircleRegression::position(experimental_double &phi){

      return c_.position(phi);
    }

    bool CircleRegression::points_in_good_order(void){
      // check the points are monotonically moving along the circle

      size_t s = xi_.size();
      if( s < 3 ) return true;

      experimental_double deltaphiAB, deltaphi_product;
      size_t index;

      experimental_double phi_second =  c_.phi_of_point(experimental_point(xi_[1], experimental_double(0.,0.), yi_[1]));
      experimental_double phi_initial =  c_.phi_of_point(experimental_point(xi_[0], experimental_double(0.,0.), yi_[0]), phi_second.value());
      experimental_double phi_final_minus_one =  c_.phi_of_point(experimental_point(xi_[s-2], experimental_double(0.,0.), yi_[s-2]), phi_initial.value());
      experimental_double phi_final =  c_.phi_of_point(experimental_point(xi_.back(), experimental_double(0.,0.), yi_.back()), phi_final_minus_one.value());
      experimental_double deltaphi_overall =  phi_final - phi_initial;
      experimental_double phiA = phi_initial;
      experimental_double phiB = phi_initial;

      for(std::vector<experimental_double>::iterator it=xi_.begin(); it != xi_.end(); ++it)
        {
          index = it - xi_.begin();
          if( index >= 1 ){
            phiA = c_.phi_of_point(experimental_point(xi_[index-1], experimental_double(0.,0.), yi_[index-1]), phiA.value());
            phiB = c_.phi_of_point(experimental_point(*it, experimental_double(0.,0.), yi_[index]), phiA.value());
            deltaphiAB = phiB-phiA;
            deltaphi_product = deltaphiAB*deltaphi_overall;
            if( deltaphi_product.value() < - deltaphi_product.error() ){
              // if( print_level() >= mybhep::VVERBOSE ){
              //   std::clog << "CAT::CircleRegression::points_in_good_order: points are not in good order: phi[ " << index-1 << "] = " << phiA.value() << ", phi[" << index << "] = " << phiB.value() << ", deltaphi_overall = " << deltaphi_overall.value() << std::endl;
              // }
              return false;
            }
          }
        }
      return true;
    }

  }

}
