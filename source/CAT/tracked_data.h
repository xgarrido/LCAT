/* -*- mode: c++ -*- */
#ifndef __CATAlgorithm__ITRACKED_DATA
#define __CATAlgorithm__ITRACKED_DATA
#include <iostream>
#include <cmath>
#include <vector>
#include <CAT/utilities.h>
#include <CAT/experimental_point.h>
#include <CAT/experimental_vector.h>
#include <CAT/cell_base.h>
#include <CAT/line.h>
#include <CAT/cell_couplet.h>
#include <CAT/cell_triplet.h>
#include <CAT/node.h>
#include <CAT/calorimeter_hit.h>
#include <CAT/scenario.h>
// #include <CAT/lt_ utils.h>


namespace CAT{
  namespace topology{


    class tracked_data : public tracking_object{

      // a tracked_data is composed of a list of cells
      // a list of clusters
      // and a list of scenarios

    public:

      // list of cells
      std::vector<cell> cells_;

      // list of calos
      std::vector<calorimeter_hit> calos_;

      // list of clusters
      std::vector<cluster> clusters_;

      // list of scenarios
      std::vector<scenario> scenarios_;

      //!Default constructor
      tracked_data()
      {
      }

      //!Default destructor
      virtual ~tracked_data(){}

      /*** dump ***/
      virtual void dump (std::ostream & a_out         = std::clog,
                         const std::string & a_title  = "",
                         const std::string & a_indent = "",
                         bool /*a_inherit*/          = false) const {
        std::string indent;
        if (! a_indent.empty ()) indent = a_indent;
        if (! a_title.empty ())
          {
            a_out << indent << a_title << std::endl;
          }

        a_out << indent << " number of cells : " << cells_.size() << std::endl;
        for(std::vector<cell>::const_iterator icell=cells_.begin(); icell!=cells_.end();++icell)
          icell->dump(a_out, "",indent + "     ");
        a_out << indent << " number of calos : " << calos_.size() << std::endl;
        for(std::vector<calorimeter_hit>::const_iterator icalo=calos_.begin(); icalo!=calos_.end();++icalo)
          icalo->dump(a_out, "",indent + "     ");
        a_out << indent << " number of clusters : " << clusters_.size() << std::endl;
        for(std::vector<cluster>::const_iterator icluster=clusters_.begin(); icluster!=clusters_.end();++icluster)
          icluster->dump(a_out, "",indent + "     ");
        a_out << indent << " number of scenarios : " << scenarios_.size() << std::endl;
        for(std::vector<scenario>::const_iterator iscenario=scenarios_.begin(); iscenario!=scenarios_.end();++iscenario)
          iscenario->dump(a_out, "",indent + "     ");
        a_out << indent << " ------------------- " << std::endl;

        return;
      }

      //! set cells
      void set_cells(const std::vector<cell> & cells)
      {
        cells_ = cells;
      }

      //! set calos
      void set_calos(const std::vector<calorimeter_hit> & calos)
      {
        calos_ = calos;
      }

      //! set clusters
      void set_clusters(const std::vector<cluster> & clusters)
      {
        clusters_ = clusters;
      }

      //! set scenarios
      void set_scenarios(const std::vector<scenario> & scenarios)
      {
        scenarios_ = scenarios;
      }

      /// Add a Geiger cell
      topology::cell & add_gg_cell();

      //! get cells
      std::vector<cell>& get_cells()
      {
        return cells_;
      }

      const std::vector<cell>& get_cells()const
      {
        return cells_;
      }

      /// Add a calorimeter cell
      topology::calorimeter_hit & add_calo_cell();

      //! get calos
      std::vector<calorimeter_hit>& get_calos()
      {
        return calos_;
      }

      const std::vector<calorimeter_hit>& get_calos()const
      {
        return calos_;
      }

      //! get clusters
      std::vector<cluster>& get_clusters()
      {
        return clusters_;
      }

      const std::vector<cluster>& get_clusters()const
      {
        return clusters_;
      }

      //! get scenarios
      std::vector<scenario>& get_scenarios()
      {
        return scenarios_;
      }

      const std::vector<scenario>& get_scenarios()const
      {
        return scenarios_;
      }

      void reset(){
        cells_.clear();
        clusters_.clear();
        scenarios_.clear();
      }

    };
  }
}

#endif
