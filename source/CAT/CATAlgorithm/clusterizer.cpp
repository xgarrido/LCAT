#include <CATAlgorithm/clusterizer.h>
#include <sys/time.h>
#include <limits>
#include <cmath>
#include <map>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/exception.h>

namespace CAT {

  void clusterizer::set_logging_priority(datatools::logger::priority a_priority)
  {
    _logging_ = a_priority;
    return;
  }

  datatools::logger::priority clusterizer::get_logging_priority() const
  {
    return _logging_;
  }

  //! get cells
  const std::vector<topology::cell>& clusterizer::get_cells()const
  {
    return cells_;
  }

  //! set cells
  void clusterizer::set_cells(const std::vector<topology::cell> & cells)
  {
    cells_ = cells;
  }

  //! get clusters
  const std::vector<topology::cluster>& clusterizer::get_clusters()const
  {
    return clusters_;
  }

  //! set clusters
  void clusterizer::set_clusters(const std::vector<topology::cluster> & clusters)
  {
    clusters_ = clusters;
  }

  //! get calorimeter_hits
  const std::vector<topology::calorimeter_hit>& clusterizer::get_calorimeter_hits()const
  {
    return calorimeter_hits_;
  }

  //! set calorimeter_hits
  void clusterizer::set_calorimeter_hits(const std::vector<topology::calorimeter_hit> & calorimeter_hits)
  {
    calorimeter_hits_ = calorimeter_hits;
  }

  void clusterizer::_set_initialized(bool i_)
  {
    _initialized_ = i_;
    return;
  }

  void clusterizer::_set_defaults()
  {
    _logging_ = datatools::logger::PRIO_WARNING;
    num_blocks = -1;
    planes_per_block.clear ();
    gaps_Z.clear ();
    num_cells_per_plane = -1;
    calo_X = calo_Y = calo_Z = std::numeric_limits<double>::quiet_NaN ();
    pmax = std::numeric_limits<double>::quiet_NaN ();

    SmallRadius = std::numeric_limits<double>::quiet_NaN ();
    TangentPhi = std::numeric_limits<double>::quiet_NaN ();
    TangentTheta = std::numeric_limits<double>::quiet_NaN ();
    SmallNumber = std::numeric_limits<double>::quiet_NaN ();
    QuadrantAngle = std::numeric_limits<double>::quiet_NaN ();
    Ratio = std::numeric_limits<double>::quiet_NaN ();
    CompatibilityDistance = std::numeric_limits<double>::quiet_NaN ();
    MaxChi2 = std::numeric_limits<double>::quiet_NaN ();
    probmin = std::numeric_limits<double>::quiet_NaN ();
    nofflayers = 0;
    SuperNemo = true;
    SuperNemoChannel = false;
    NemoraOutput = false;
    N3_MC = false;
    _moduleNR.clear ();
    _MaxBlockSize = -1;

    nevent = 0;
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

  void clusterizer::prepare_event(topology::tracked_data & tracked_data_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");

    clusters_.clear();

    order_cells();

    tracked_data_.set_cells(cells_);
    tracked_data_.set_calos(calorimeter_hits_);


    return;
  }


  void clusterizer::clusterize(topology::tracked_data & tracked_data_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");

    DT_LOG_DEBUG(get_logging_priority(), "Fill clusters");

    if (cells_.empty()) return;

    std::map<int,bool> flags;
    for(size_t i=0; i<cells_.size(); i++)
      {
        flags[cells_[i].id()] = false;
      }

    for (std::vector<topology::cell>::const_iterator icell = cells_.begin();
         icell != cells_.end(); ++icell) {
      // pick a cell c that was never added
      const topology::cell & a_cell = *icell;
      if (flags[a_cell.id()]) continue;
      flags[a_cell.id()] = true;

      // cell c will form a new cluster, i.e. a new list of nodes
      topology::cluster cluster_connected_to_c;
      std::vector<topology::node> nodes_connected_to_c;
      DT_LOG_DEBUG(get_logging_priority(), "Begin new cluster with cell " << a_cell.id());

      // let's get the list of all the cells that can be reached from c without
      // jumps
      std::vector<topology::cell> cells_connected_to_c;
      cells_connected_to_c.push_back(a_cell);

      for( size_t i=0; i<cells_connected_to_c.size(); i++){ // loop on connected cells
        DT_LOG_WARNING(datatools::logger::PRIO_WARNING, "i=" << i);
        // take a connected cell (the first one is just c)
        topology::cell cconn = cells_connected_to_c[i];

        // the connected cell composes a new node
        topology::node newnode(cconn);
        std::vector<topology::cell_couplet> cc;

        // get the list of cells near the connected cell
        std::vector<topology::cell> cells_near_iconn = get_near_cells(cconn);

        DT_LOG_DEBUG(get_logging_priority(), "Cluster " << clusters_.size()
                     << " starts with " << a_cell.id() << " try to add cell " << cconn.id()
                     << " with n of neighbours = " << cells_near_iconn.size());
        for(std::vector<topology::cell>::const_iterator icnc = cells_near_iconn.begin();
            icnc != cells_near_iconn.end(); ++icnc) {
          const topology::cell & cnc = *icnc;

          if (!_is_good_couplet_(cconn, cnc, cells_near_iconn)) continue;

          topology::cell_couplet ccnc(cconn,cnc);
          cc.push_back(ccnc);

          DT_LOG_DEBUG(get_logging_priority(), "Creating couplet " << cconn.id() << " -> " << cnc.id());

          if(! flags[cnc.id()]) {
            flags[cnc.id()] = true;
            cells_connected_to_c.push_back(cnc);
          }
        }
        newnode.set_cc(cc);
        newnode.calculate_triplets(Ratio, QuadrantAngle, TangentPhi, TangentTheta);
        nodes_connected_to_c.push_back(newnode);

        DT_LOG_DEBUG(get_logging_priority(), "Cluster started with " << a_cell.id()
                     << " has been given cell " << cconn.id() << " with " << cc.size() << " couplets ");

      }

      cluster_connected_to_c.set_nodes(nodes_connected_to_c);
      clusters_.push_back(cluster_connected_to_c);
    }


    DT_LOG_DEBUG(get_logging_priority(), "There are " << clusters_.size() << " clusters of cells");

    tracked_data_.set_cells(cells_);
    tracked_data_.set_clusters(clusters_);

    return;
  }

  bool clusterizer::_is_good_couplet_(const topology::cell & main_cell_,
                                      const topology::cell & candidate_cell_,
                                      const std::vector<topology::cell> & cells_near_main_)
  {
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
      // ccc.set_probmin(probmin);
      ccc.calculate_joints(Ratio, QuadrantAngle, TangentPhi, TangentTheta);
      if (ccc.joints().size() > 0) {
        DT_LOG_DEBUG(get_logging_priority(),
                     "... yes it does: so couplet " << a.id() << " and " << c.id() << " is not good");
        return false;
      }
    }
    return true;
  }


  //*************************************************************
  int clusterizer::cell_side( const topology::cell & c){
    //*************************************************************

    if( c.ep().z().value() > 0. )
      return 1;

    return -1;

  }


  size_t clusterizer::near_level( const topology::cell & c1, const topology::cell & c2 ){

    // returns 0 for far-away cell
    // 1 for diagonal cells
    // 2 for side-by-side cells

    // side-by-side connection: distance = 1
    // diagonal connection: distance = sqrt(2) = 1.41
    // skip 1 connection, side: distance = 2
    // skip 1 connection, tilt: distance = sqrt(5) = 2.24
    // skip 1 connection, diag: distance = 2 sqrt(2) = 2.83

    topology::experimental_double distance = topology::experimental_vector(c1.ep(),c2.ep()).hor().length();

    // Use geiger locator for such research Warning: use integer
    // because uint32_t has strange behavior with absolute value
    // cmath::abs
    const int hit1_side  = c1.block();  // -1, 1
    const int hit1_layer = abs(c1.layer()); // 0, 1, ..., 8
    const int hit1_row   = c1.iid();  // -56, -55, ..., 55, 56

    const int hit2_side  = c2.block();
    const int hit2_layer = abs(c2.layer());
    const int hit2_row   = c2.iid();

    // Do not cross the foil
    if (hit1_side != hit2_side) return 0;

    // Check neighboring
    const unsigned int layer_distance = abs (hit1_layer - hit2_layer); // 1 --> side-by-side
    const unsigned int row_distance = abs (hit1_row - hit2_row);

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


  std::vector<topology::cell> clusterizer::get_near_cells(const topology::cell & c){


    DT_LOG_DEBUG(get_logging_priority(),
                 "Filling list of cells near cell " << c.id() << " fast " << c.fast() << " side " << cell_side(c));

    std::vector<topology::cell> cells;

    for(std::vector<topology::cell>::iterator kcell=cells_.begin(); kcell != cells_.end(); ++kcell){
      if( kcell->id() == c.id() ) continue;

      if( kcell->fast() != c.fast() ) continue;

      if( cell_side(*kcell) != cell_side(c) ) continue;

      size_t nl = near_level(c,*kcell);

      if( nl > 0 )
        {

          topology::cell ck = *kcell;
          cells.push_back(ck);
        }
    }


    return cells;

  }


  //*************************************************************
  void clusterizer::order_cells(){
    //*************************************************************


    //  std::sort( cells_.begin(), cells_.end(), topology::cell::compare );
    std::sort( cells_.begin(), cells_.end());


    return;

  }

  void clusterizer::set_num_blocks(int nb){
    if (nb > 0)
      {
        num_blocks = nb;
        planes_per_block.assign (num_blocks, 1);
      }
    else
      {
        std::cerr << "WARNING: CAT::clusterizer::set_num_blocks: "
                  << "Invalid number of GG layer blocks !" << std::endl;
        planes_per_block.clear ();
        num_blocks = -1; // invalid value
      }
    return;
  }

  void clusterizer::set_planes_per_block(int block, int nplanes){
    if (block< 0 || block>= (int)planes_per_block.size())
      {
        throw std::range_error ("CAT::clusterizer::set_planes_per_block: Invalid GG layer block index !");
      }
    if (nplanes > 0)
      {
        planes_per_block.at (block) = nplanes;
      }
    else
      {
        throw std::range_error ("CAT::clusterizer::set_planes_per_block: Invalid number of GG layers in block !");
      }
    return;
  }

  void clusterizer::set_num_cells_per_plane(int ncpp){
    if (ncpp <= 0)
      {
        num_cells_per_plane = -1; // invalid value
      }
    else
      {
        num_cells_per_plane = ncpp;
      }
    return;
  }

  void clusterizer::set_pmax(double v){
    if ( v <= 0.0)
      {
        pmax = std::numeric_limits<double>::quiet_NaN ();
      }
    else
      {
        pmax = v;
      }
    return;
  }

  void clusterizer::set_SmallRadius(double v){
    SmallRadius = v;
    return;
  }

  void clusterizer::set_TangentPhi(double v){
    TangentPhi = v;
    return;
  }

  void clusterizer::set_TangentTheta(double v){
    TangentTheta = v;
    return;
  }

  void clusterizer::set_SmallNumber(double v){
    SmallNumber = v;
    return;
  }

  void clusterizer::set_QuadrantAngle(double v){
    QuadrantAngle = v;
    return;
  }

  void clusterizer::set_Ratio(double v){
    Ratio = v;
    return;
  }

  void clusterizer::set_CompatibilityDistance(double v){
    CompatibilityDistance = v;
    return;
  }

  void clusterizer::set_MaxChi2(double v){
    MaxChi2 = v;
    return;
  }

  void clusterizer::set_probmin(double v){
    probmin = v;
    return;
  }

  void clusterizer::set_nofflayers(size_t v){
    nofflayers = v;
    return;
  }

  void clusterizer::set_SuperNemo(bool v){
    SuperNemo = v;
    return;
  }

  void clusterizer::set_SuperNemoChannel(bool v){
    if (v)
      {
        set_SuperNemo (true);
        SuperNemoChannel = true;
      }
    else
      {
        SuperNemoChannel = false;
      }
    return;
  }

  void clusterizer::set_bfield(double v){
    bfield = v;
    return;
  }

}
