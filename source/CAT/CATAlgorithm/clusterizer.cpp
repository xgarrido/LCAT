// Ourselves
#include <CATAlgorithm/clusterizer.h>

// Standard library
#include <set>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/exception.h>
#include <bayeux/datatools/utils.h>

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

  //! get cells
  const std::vector<topology::cell>& clusterizer::get_cells()const
  {
    return _cells_;
  }

  //! set cells
  void clusterizer::set_cells(const std::vector<topology::cell> & cells_)
  {
    _cells_ = cells_;
  }

  //! get clusters
  const std::vector<topology::cluster>& clusterizer::get_clusters()const
  {
    return _clusters_;
  }

  //! set clusters
  void clusterizer::set_clusters(const std::vector<topology::cluster> & clusters_)
  {
    _clusters_ = clusters_;
  }

  //! get calorimeter_hits
  const std::vector<topology::calorimeter_hit>& clusterizer::get_calorimeter_hits()const
  {
    return _calorimeter_hits_;
  }

  //! set calorimeter_hits
  void clusterizer::set_calorimeter_hits(const std::vector<topology::calorimeter_hit> & calorimeter_hits_)
  {
    _calorimeter_hits_ = calorimeter_hits_;
  }

  void clusterizer::_set_initialized(bool i_)
  {
    _initialized_ = i_;
    return;
  }

  void clusterizer::_set_defaults()
  {
    _logging_ = datatools::logger::PRIO_WARNING;

    datatools::invalidate(_tangent_phi_);
    datatools::invalidate(_tangent_theta_);
    datatools::invalidate(_quadrant_angle_);
    datatools::invalidate(_ratio_);

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

  void clusterizer::initialize()
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(is_initialized(), std::logic_error, "Already initialized !");
    _set_initialized(true);
    DT_LOG_TRACE(get_logging_priority(), "Entering.");
    return;
  }


  void clusterizer::reset()
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    _set_defaults();
    _set_initialized(false);
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return;
  }

  void clusterizer::clusterize(topology::tracked_data & tracked_data_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");

    DT_LOG_DEBUG(get_logging_priority(), "Fill clusters");
    _clusters_.clear();

    DT_LOG_DEBUG(get_logging_priority(), "Order cells");
    if (_cells_.empty()) return;
    std::sort(_cells_.begin(), _cells_.end());

    // List of cells already used
    std::set<int> flagged;

    for (const auto & icell : _cells_) {
      // Pick a cell that was never added
      if (flagged.count(icell.id())) continue;
      flagged.insert(icell.id());

      DT_LOG_DEBUG(get_logging_priority(), "Begin new cluster with cell " << icell.id());
      // Current cell will form a new cluster, i.e. a new list of nodes
      topology::cluster cluster_connected_to_c;
      std::vector<topology::node> nodes_connected_to_c;

      // Let's get the list of all the cells that can be reached from current
      // cell without jumps
      std::vector<topology::cell> cells_connected_to_c;
      cells_connected_to_c.push_back(icell);

      // Loop on connected cells
      for (size_t i = 0; i < cells_connected_to_c.size(); i++) {
        // Take a connected cell by value since the array size will change and
        // thus the allocation space too : reference will change during the loop!
        const topology::cell cconn = cells_connected_to_c[i];

        // The connected cell composes a new node
        topology::node newnode(cconn);
        std::vector<topology::cell_couplet> cc;

        // Get the list of cells near the connected cell
        std::vector<topology::cell> cells_near_iconn;
        get_near_cells(cconn, cells_near_iconn);
        DT_LOG_DEBUG(get_logging_priority(), "Cluster " << _clusters_.size()
                     << " starts with " << icell.id() << " try to add cell " << cconn.id()
                     << " with n of neighbours = " << cells_near_iconn.size());
        for (const auto & cnc : cells_near_iconn) {
          if (!_is_good_couplet_(cconn, cnc, cells_near_iconn)) continue;

          DT_LOG_DEBUG(get_logging_priority(), "Creating couplet " << cconn.id() << " -> " << cnc.id());
          topology::cell_couplet ccnc(cconn, cnc);
          cc.push_back(ccnc);

          if (! flagged.count(cnc.id())) {
            flagged.insert(cnc.id());
            cells_connected_to_c.push_back(cnc);
          }
        }
        newnode.set_cc(cc);
        newnode.calculate_triplets(_ratio_, _quadrant_angle_, _tangent_phi_, _tangent_theta_);
        nodes_connected_to_c.push_back(newnode);

        DT_LOG_DEBUG(get_logging_priority(), "Cluster started with " << icell.id()
                     << " has been given cell " << cconn.id() << " with " << cc.size() << " couplets ");
      }
      cluster_connected_to_c.set_nodes(nodes_connected_to_c);
      _clusters_.push_back(cluster_connected_to_c);
    }


    DT_LOG_DEBUG(get_logging_priority(), "There are " << _clusters_.size() << " clusters of cells");

    tracked_data_.set_cells(_cells_);
    tracked_data_.set_calos(_calorimeter_hits_);
    tracked_data_.set_clusters(_clusters_);

    return;
  }

  bool clusterizer::_is_good_couplet_(const topology::cell & main_cell_,
                                      const topology::cell & candidate_cell_,
                                      const std::vector<topology::cell> & cells_near_main_) const
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    // the couplet mainc -> candidatec is good only if
    // there is no other cell that is near to both and can form a triplet between them

    const topology::cell & a = main_cell_;
    const topology::cell & c = candidate_cell_;

    for (std::vector<topology::cell>::const_iterator icell = cells_near_main_.begin();
         icell != cells_near_main_.end(); ++icell) {
      const topology::cell & b = *icell;

      if (b.id() == c.id()) continue;
      if (near_level(b, c) == 0) continue;

      if (near_level(b, c) < near_level(a, c) ||
          near_level(b, a) < near_level(a, c))
        continue;  // cannot match a->b or b->c if a->c is nearer

      DT_LOG_DEBUG(get_logging_priority(),
                   "... check if near node " << b.id() << " has triplet " << a.id() << " <-> " << c.id());

      topology::cell_triplet ccc(a,b,c);
      ccc.calculate_joints(_ratio_, _quadrant_angle_, _tangent_phi_, _tangent_theta_);
      if (ccc.joints().size() > 0) {
        DT_LOG_DEBUG(get_logging_priority(),
                     "... yes it does: so couplet " << a.id() << " and " << c.id() << " is not good");
        return false;
      }
    }
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return true;
  }


  int clusterizer::cell_side(const topology::cell & c_) const
  {
    if (c_.ep().z().value() > 0.)
      return 1;

    return -1;
  }


  size_t clusterizer::near_level(const topology::cell & c1_, const topology::cell & c2_) const
  {
    // returns 0 for far-away cell
    // 1 for diagonal cells
    // 2 for side-by-side cells

    // Use geiger locator for such research
    const int hit1_side  = c1_.block();  // -1, 1
    const int hit1_layer = std::abs(c1_.layer()); // 0, 1, ..., 8
    const int hit1_row   = c1_.iid();  // -56, -55, ..., 55, 56

    const int hit2_side  = c2_.block();
    const int hit2_layer = std::abs(c2_.layer());
    const int hit2_row   = c2_.iid();

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


  void clusterizer::get_near_cells(const topology::cell & c_, std::vector<topology::cell> & cells_) const
  {
    DT_LOG_DEBUG(get_logging_priority(), "Filling list of cells near cell " << c_.id());

    for (const auto & icell : _cells_) {
      if (icell.id() == c_.id()) continue;

      if (cell_side(icell) != cell_side(c_)) continue;

      const size_t nl = near_level(c_, icell);
      if (nl > 0) {
        cells_.push_back(icell);
      }
    }
    return;
  }

  void clusterizer::set_tangent_phi(double phi_){
    _tangent_phi_ = phi_;
    return;
  }

  void clusterizer::set_tangent_theta(double theta_){
    _tangent_theta_ = theta_;
    return;
  }

  void clusterizer::set_quadrant_angle(double angle_){
    _quadrant_angle_ = angle_;
    return;
  }

  void clusterizer::set_ratio(double ratio_){
    _ratio_ = ratio_;
    return;
  }

}
