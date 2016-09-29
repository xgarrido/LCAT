// -*- mode: c++ -*-

#ifndef CAT_TOPOLOGY_CLUSTER_H
#define CAT_TOPOLOGY_CLUSTER_H

// Standard library:
#include <iostream>
#include <cmath>
#include <vector>

// This project:
#include <CAT/experimental_point.h>
#include <CAT/experimental_vector.h>
#include <CAT/cell_base.h>
#include <CAT/line.h>
#include <CAT/cell_couplet.h>
#include <CAT/cell_triplet.h>
#include <CAT/node.h>
#include <CAT/broken_line.h>


namespace CAT {

  namespace topology {

    /// \brief A cluster is composed of a list of nodes
    class cluster : public tracking_object
    {
    public:

      // list of nodes
      std::vector<node> nodes_;

      // status of cluster
      bool free_;

      //!Default constructor
      cluster();

      //!Default destructor
      virtual ~cluster();

      //! constructor from std::vector of nodes
      cluster(const std::vector<node> &nodes, double probmin=1.e-200);

      //! constructor from single node
      cluster(node &a_node, double probmin=1.e-200);

      /*** dump ***/
      virtual void dump (std::ostream & a_out         = std::clog,
                         const std::string & a_title  = "",
                         const std::string & a_indent = "",
                         bool a_inherit          = false) const ;
      //! set nodes
      void set_nodes(const std::vector<node> &nodes);

      //! set free level
      void set_free(bool free);

      //! get nodes
      const std::vector<node> & nodes()const;

      //! get free level
      bool Free()const;


    public:

      bool has_cell(const cell & c)const;


      cluster invert();

      topology::node node_of_cell(const topology::cell & c);


      void solve_ambiguities(std::vector< std::vector<topology::broken_line> > * sets_of_bl_alternatives);

      bool start_ambiguity(size_t i);

      bool end_ambiguity(size_t i);

      std::vector<topology::broken_line> solve_ambiguities_with_ends(size_t ifirst, size_t ilast);

      std::vector<topology::broken_line> solve_ambiguities_with_ends__1_node(size_t ifirst, size_t ilast, bool first_ambiguous_is_after_gap, bool first_ambiguous_is_second, bool last_ambiguous_is_begore_gap, bool last_ambiguous_is_last_but_one);

      std::vector<topology::broken_line> solve_ambiguities_with_ends__2_nodes(size_t ifirst, size_t ilast, bool first_ambiguous_is_after_gap, bool first_ambiguous_is_second, bool last_ambiguous_is_begore_gap, bool last_ambiguous_is_last_but_one);

      std::vector<topology::broken_line> solve_ambiguities_with_ends__3_nodes(size_t ifirst, size_t ilast, bool first_ambiguous_is_after_gap, bool first_ambiguous_is_second, bool last_ambiguous_is_begore_gap, bool last_ambiguous_is_last_but_one);

      std::vector<topology::broken_line> solve_ambiguities_with_ends__4_nodes(size_t ifirst, size_t ilast, bool first_ambiguous_is_after_gap, bool first_ambiguous_is_second, bool last_ambiguous_is_begore_gap, bool last_ambiguous_is_last_but_one);

      void solve_ambiguities_with_ends__more_than_4_nodes(topology::broken_line ACD[2][2][2], size_t ifirst, size_t ilast, bool first_ambiguous_is_after_gap, bool first_ambiguous_is_second);

      void solve_ambiguities_with_ends__more_than_4_nodes(topology::broken_line aACD[2][2][2][2], size_t ifirst, size_t ilast);

      void merge__more_than_4_nodes(topology::broken_line ACD[2][2][2], topology::broken_line aACD[2][2][2][2]);

      std::vector<topology::broken_line> finish__more_than_4_nodes(topology::broken_line ACD[2][2][2], size_t ipivot, size_t n_residuals);

    };

  } // namespace topology

} // namespace CAT

#endif // CAT_TOPOLOGY_CLUSTER_H
