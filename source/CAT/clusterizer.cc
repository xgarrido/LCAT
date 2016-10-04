// Ourselves
#include <CAT/clusterizer.h>
#include <CAT/tracked_data.h>

// Standard library
#include <set>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/exception.h>
#include <bayeux/datatools/utils.h>
#include <bayeux/datatools/clhep_units.h>
#include <bayeux/datatools/properties.h>

namespace CAT {

  void clusterizer::set_logging_priority(datatools::logger::priority priority_)
  {
    _logging_ = priority_;
    return;
  }

  datatools::logger::priority clusterizer::get_logging_priority() const
  {
    return _logging_;
  }

  void clusterizer::_set_initialized(bool i_)
  {
    _initialized_ = i_;
    return;
  }

  void clusterizer::_set_defaults()
  {
    _logging_ = datatools::logger::PRIO_WARNING;
    _tangent_phi_ = 20.0 * CLHEP::degree;
    _ratio_ = 10000.0;
    return;
  }

  bool clusterizer::is_initialized() const
  {
    return _initialized_;
  }

  // Default constructor :
  clusterizer::clusterizer()
  {
    _set_initialized(false);
    _set_defaults();
    return;
  }

  clusterizer::~clusterizer()
  {
    if (is_initialized()) {
      reset();
    }
    return;
  }

  void clusterizer::initialize(const datatools::properties & setup_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(is_initialized(), std::logic_error, "Already initialized !");

    //Logging priority
    datatools::logger::priority lp = datatools::logger::extract_logging_configuration(setup_);
    DT_THROW_IF(lp == datatools::logger::PRIO_UNDEFINED, std::logic_error,
                "Invalid logging priority level for clusterizer !");
    set_logging_priority(lp);

    std::string key;
    if (setup_.has_key(key = "tangent_phi")) {
      _tangent_phi_ = setup_.fetch_real(key);
      if (! setup_.has_explicit_unit(key)) {
        _tangent_phi_ *= CLHEP::degree;
      }
    }

    if (setup_.has_key(key = "ratio")) {
      _ratio_ = setup_.fetch_real(key);
    }

    _set_initialized(true);
    DT_LOG_TRACE(get_logging_priority(), "Entering.");
    return;
  }


  void clusterizer::reset()
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(! is_initialized(), std::logic_error, "Clusterizer is not initialized !");
    _set_defaults();
    _set_initialized(false);
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return;
  }

  void clusterizer::clusterize(tracked_data & tracked_data_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(! is_initialized(), std::logic_error, "Clusterizer is not initialized !");

    DT_LOG_DEBUG(get_logging_priority(), "Order cells");
    std::vector<cell> & the_cells = tracked_data_.grab_gg_hits();
    if (the_cells.empty()) return;
    std::sort(the_cells.begin(), the_cells.end(),
              [] (const cell & a_, const cell & b_) -> bool
              {
                if (a_.get_id() == b_.get_id()) return false;
                if (a_.get_id() > default_integer || b_.get_id() > default_integer) {
                  return false;
                }
                // side of foil
                if (a_.get_side() < 0 && b_.get_side() > 0) {
                  return false;
                }
                if (a_.get_side() > 0 && b_.get_side() < 0) {
                  return true;
                }
                // layer
                if (std::abs(a_.get_layer()) < std::abs(b_.get_layer())) {
                  return false;
                }
                // row
                if(a_.get_row() < b_.get_row()) {
                  return false;
                }
                return true;
              });

    DT_LOG_DEBUG(get_logging_priority(), "Fill clusters");
    std::vector<cluster> & the_clusters = tracked_data_.grab_clusters();
    the_clusters.clear();

    // List of cells already used
    std::set<int> flagged;
    for (const auto & icell : the_cells) {
      // Pick a cell that was never added
      if (flagged.count(icell.get_id())) continue;
      flagged.insert(icell.get_id());

      DT_LOG_DEBUG(get_logging_priority(), "Begin new cluster with cell " << icell.get_id());
      // Current cell will form a new cluster, i.e. a new list of nodes
      cluster cluster_connected_to_c;
      std::vector<node> nodes_connected_to_c;

      // Let's get the list of all the cells that can be reached from current
      // cell without jumps
      std::vector<cell> cells_connected_to_c;
      cells_connected_to_c.push_back(icell);

      // Loop on connected cells
      for (size_t i = 0; i < cells_connected_to_c.size(); i++) {
        // Take a connected cell by value since the array size will change and
        // thus the allocation space too : reference will change during the loop!
        const cell cconn = cells_connected_to_c[i];

        // The connected cell composes a new node
        node newnode(cconn);
        std::vector<cell_couplet> cc;

        // Get the list of cells near the connected cell
        DT_LOG_DEBUG(get_logging_priority(), "Filling list of cells near cell " << cconn.get_id());
        std::vector<cell> cells_near_iconn;
        cells_near_iconn.reserve(8);
        for (const auto & jcell : the_cells) {
          if (jcell.get_id() == cconn.get_id()) continue;
          const size_t nl = _near_level_(cconn, jcell);
          if (nl > 0) {
            cells_near_iconn.push_back(jcell);
          }
        }
        cc.reserve(cells_near_iconn.size());

        DT_LOG_DEBUG(get_logging_priority(), "Cluster " << the_clusters.size()
                     << " starts with " << icell.get_id() << " try to add cell " << cconn.get_id()
                     << " with n of neighbours = " << cells_near_iconn.size());
        cells_connected_to_c.reserve(cells_near_iconn.size());
        for (const auto & cnc : cells_near_iconn) {
          if (!_is_good_couplet_(cconn, cnc, cells_near_iconn)) continue;

          DT_LOG_DEBUG(get_logging_priority(), "Creating couplet " << cconn.get_id() << " -> " << cnc.get_id());
          cell_couplet ccnc(cconn, cnc);
          cc.push_back(ccnc);
          if (! flagged.count(cnc.get_id())) {
            flagged.insert(cnc.get_id());
            cells_connected_to_c.push_back(cnc);
          }
        }
        newnode.set_cc(cc);
        newnode.calculate_triplets(_ratio_, _tangent_phi_);
        nodes_connected_to_c.push_back(newnode);

        DT_LOG_DEBUG(get_logging_priority(), "Cluster started with " << icell.get_id()
                     << " has been given cell " << cconn.get_id() << " with " << cc.size() << " couplets ");
      }
      cluster_connected_to_c.set_nodes(nodes_connected_to_c);
      the_clusters.push_back(cluster_connected_to_c);
    }


    DT_LOG_DEBUG(get_logging_priority(), "There are " << the_clusters.size() << " clusters of cells");
    return;
  }

  bool clusterizer::_is_good_couplet_(const cell & c1_, const cell & c2_,
                                      const std::vector<cell> & c1_neighbors_) const
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    // the couplet mainc -> candidatec is good only if there is no other cell
    // that is near to both and can form a triplet between them

    for (const auto & c3 : c1_neighbors_) {
      if (c2_.get_id() == c3.get_id()) continue;

      const size_t level23 = _near_level_(c2_, c3);
      if (level23 == 0) continue;

      const size_t level12 = _near_level_(c1_, c2_);
      const size_t level13 = _near_level_(c1_, c3);
      if (level23 < level12 || level13 < level12)
        continue;  // cannot match c1->c3 or c3->c2 if c1->c2 is nearer

      DT_LOG_DEBUG(get_logging_priority(),
                   "... check if near node " << c3.get_id() << " has triplet " << c1_.get_id() << " <-> " << c2_.get_id());
      cell_triplet ccc(c1_, c3, c2_);
      ccc.calculate_joints(_ratio_, _tangent_phi_);
      if (ccc.joints().size() > 0) {
        DT_LOG_DEBUG(get_logging_priority(),
                     "... yes it does: so couplet " << c1_.get_id() << " and " << c2_.get_id() << " is not good");
        return false;
      }
    }
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return true;
  }

  size_t clusterizer::_near_level_(const cell & c1_, const cell & c2_) const
  {
    // returns 0 for far-away cell
    // 1 for diagonal cells
    // 2 for side-by-side cells

    // Use geiger locator for such research
    const int hit1_side  = c1_.get_side();  // -1, 1
    const int hit1_layer = std::abs(c1_.get_layer()); // 0, 1, ..., 8
    const int hit1_row   = c1_.get_row();  // -56, -55, ..., 55, 56

    const int hit2_side  = c2_.get_side();
    const int hit2_layer = std::abs(c2_.get_layer());
    const int hit2_row   = c2_.get_row();

    // Do not cross the foil
    if (hit1_side != hit2_side) return 0;

    // Check neighboring
    const unsigned int layer_distance = std::abs(hit1_layer - hit2_layer); // 1 --> side-by-side
    const unsigned int row_distance = std::abs(hit1_row - hit2_row);

    if (layer_distance == 0 && row_distance == 0){
      DT_LOG_DEBUG(get_logging_priority(), "CAT asking near level of cells with identical position ("
                   << hit1_side << ", " << hit1_layer << ", " << hit1_row << ") ("
                   << hit2_side << ", " << hit2_layer << ", " << hit2_row << ")");
      return 3;
    }
    else if (layer_distance == 1 && row_distance == 0) return 2;
    else if (layer_distance == 0 && row_distance == 1) return 2;
    else if (layer_distance == 1 && row_distance == 1) return 1;
    return 0;
  }
}
