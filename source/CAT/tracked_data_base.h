// -*- mode: c++ -*-

#ifndef CAT_TOPOLOGY_TRACKED_DATA_H
#define CAT_TOPOLOGY_TRACKED_DATA_H

// Standard library:
#include <iostream>
#include <cmath>
#include <vector>

// This project:
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

namespace CAT{

  namespace topology{

    /// \brief Tracked data
    ///
    /// A tracked data is composed of a list of cells
    /// a list of clusters
    /// and a list of scenarios.
    class tracked_data : public tracking_object
    {
    public:

      // list of cells
      std::vector<cell> cells_;

      // list of calos
      std::vector<calorimeter_hit> calos_;

      // list of clusters
      std::vector<cluster> clusters_;

      // list of scenarios
      std::vector<scenario> scenarios_;

      // list of true sequences
      std::vector<sequence> true_sequences_;

      // list of nemo3 sequences
      std::vector<sequence> nemo_sequences_;

      // is event selcted?
      bool selected_;

      // is event skipped?
      bool skipped_;

      //! Default constructor
      tracked_data()
      {
        selected_ = true;
        skipped_ = false;
      }

      //! Destructor
      virtual ~tracked_data(){}

      //! Constructor
      tracked_data(const std::vector<cell> &cells,
                   const std::vector<calorimeter_hit> &calos,
                   const std::vector<cluster> &clusters,
                   const std::vector<scenario> &scenarios,
                   const std::vector<sequence> &true_sequences,
                   const std::vector<sequence> &nemo_sequences,
                   double probmin=1.e-200)
      {
        set_probmin(probmin);
        cells_ = cells;
        calos_ = calos;
        clusters_ = clusters;
        scenarios_ = scenarios;
        true_sequences_ = true_sequences;
        nemo_sequences_ = nemo_sequences;
        selected_ = true;
        skipped_ = false;
      }

      /// Smart dump
      void dump(std::ostream & a_out         = std::clog,
                const std::string & a_title  = "",
                const std::string & a_indent = "",
                bool /*a_inherit*/          = false) const
      {
        std::string indent;
        if (! a_indent.empty ()) indent = a_indent;
        if (! a_title.empty ()) {
          a_out << indent << a_title << std::endl;
        }

        a_out << indent << " ------------------- " << std::endl;
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

      void set_cells(const std::vector<cell> & cells)
      {
        cells_ = cells;
      }

      void set_calos(const std::vector<calorimeter_hit> & calos)
      {
        calos_ = calos;
      }

      void set_clusters(const std::vector<cluster> & clusters)
      {
        clusters_ = clusters;
      }

      void set_scenarios(const std::vector<scenario> & scenarios)
      {
        scenarios_ = scenarios;
      }

      void set_true_sequences(const std::vector<sequence> & sequences)
      {
        true_sequences_ = sequences;
      }

      void set_nemo_sequences(const std::vector<sequence> & sequences)
      {
        nemo_sequences_ = sequences;
      }

      void set_selected(bool selected)
      {
        selected_ = selected;
      }

      void set_skipped(bool skipped)
      {
        skipped_ = skipped;
      }

      std::vector<cell>& get_cells()
      {
        return cells_;
      }

      const std::vector<cell>& get_cells() const
      {
        return cells_;
      }

      std::vector<calorimeter_hit>& get_calos()
      {
        return calos_;
      }

      const std::vector<calorimeter_hit>& get_calos() const
      {
        return calos_;
      }

      std::vector<cluster>& get_clusters()
      {
        return clusters_;
      }

      const std::vector<cluster>& get_clusters() const
      {
        return clusters_;
      }

      std::vector<scenario>& get_scenarios()
      {
        return scenarios_;
      }

      const std::vector<scenario>& get_scenarios() const
      {
        return scenarios_;
      }

      std::vector<sequence>& get_true_sequences()
      {
        return true_sequences_;
      }

      const std::vector<sequence>& get_true_sequences() const
      {
        return true_sequences_;
      }

      std::vector<sequence>& get_nemo_sequences()
      {
        return nemo_sequences_;
      }

      const std::vector<sequence>& get_nemo_sequences() const
      {
        return nemo_sequences_;
      }

      bool selected() const
      {
        return selected_;
      }

      bool skipped() const
      {
        return skipped_;
      }

      void reset()
      {
        cells_.clear();
        clusters_.clear();
        scenarios_.clear();
      }

    };

  } // namespace topology

} // namespace CAT

#endif // CAT_TOPOLOGY_TRACKED_DATA_H
