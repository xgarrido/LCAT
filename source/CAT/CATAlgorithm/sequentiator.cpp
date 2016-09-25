#include "CATAlgorithm/sequentiator.h"
#include <vector>
#include <mybhep/system_of_units.h>
#include <sys/time.h>
#include <math.h>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/exception.h>

namespace CAT {

  void sequentiator::set_logging_priority(datatools::logger::priority priority_)
  {
    _logging_ = priority_;
    return;
  }

  datatools::logger::priority sequentiator::get_logging_priority() const
  {
    return _logging_;
  }

  void sequentiator::_set_initialized(bool i_)
  {
    _initialized_ = i_;
    return;
  }

  void sequentiator::_set_defaults ()
  {
    _logging_ = datatools::logger::PRIO_WARNING;
    bfield = std::numeric_limits<double>::quiet_NaN ();
    num_blocks = -1;
    planes_per_block.clear ();
    gaps_Z.clear ();
    num_cells_per_plane = -1;
    SOURCE_thick = std::numeric_limits<double>::quiet_NaN ();
    //    lastlayer = 0;
    vel  = std::numeric_limits<double>::quiet_NaN ();
    rad  = std::numeric_limits<double>::quiet_NaN ();
    len  = std::numeric_limits<double>::quiet_NaN ();
    CellDistance  = std::numeric_limits<double>::quiet_NaN ();
    xsize = ysize = zsize = std::numeric_limits<double>::quiet_NaN ();
    //    calo_X = calo_Y = calo_Z = std::numeric_limits<double>::quiet_NaN ();
    InnerRadius = OuterRadius= FoilRadius = std::numeric_limits<double>::quiet_NaN ();
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
    NOffLayers = 0;
    SuperNemo = true;
    SuperNemoChannel = false;
    NemoraOutput = false;
    N3_MC = false;
    MaxTime = std::numeric_limits<double>::quiet_NaN ();
    //    doDriftWires = true;
    //    DriftWires.clear ();

    return;
  }

  bool sequentiator::is_initialized() const
  {
    return _initialized_;
  }

  // Default constructor :
  sequentiator::sequentiator()
  {
    _set_initialized(false);
    _set_defaults();
    return;
  }

  sequentiator::~sequentiator()
  {
    if (is_initialized()) {
      reset();
    }
    return;
  }

  void sequentiator::initialize()
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(is_initialized(), std::logic_error, "Already initialized !");
    _set_initialized(true);
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return;
  }

  void sequentiator::reset()
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    // DT_THROW_IF(! is_initialized(), std::logic_error, "Sequentiator is not initialized !");
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return;
  }


  void sequentiator::sequentiate(topology::tracked_data & tracked_data_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(! is_initialized(), std::logic_error, "Sequentiator is not initialized !");

    clock.start(" sequentiator: sequentiate ","cumulative");
    clock.start(" sequentiator: sequentiation ","restart");

    // set_clusters(tracked_data_.get_clusters());
    std::vector<topology::cluster> & the_clusters = tracked_data_.get_clusters();

    NFAMILY = -1;
    NCOPY = 0;

    if (the_clusters.empty()) return;

    sequences_.clear();
    scenarios_.clear();

    tracked_data_.scenarios_.clear();

    for (std::vector<topology::cluster>::iterator
          icluster = the_clusters.begin();
        icluster != the_clusters.end(); ++icluster)
      {
        topology::cluster & a_cluster = *icluster;

        local_cluster_ = &(*icluster);

        sequentiate_cluster(a_cluster);
      }

    if (late())
      {
        tracked_data_.set_skipped(true);
        return;
      }

    clean_up_sequences();
    direct_out_of_foil();

    interpret_physics(tracked_data_.get_calos());
    make_families();

    match_gaps(tracked_data_.get_calos());
    clean_up_sequences();
    direct_out_of_foil();

    refine_sequences_near_walls(tracked_data_.get_calos());
    interpret_physics(tracked_data_.get_calos());

    if (late())
      {
        tracked_data_.set_skipped(true);
        return;
      }


    make_scenarios(tracked_data_);


    if (late())
      {
        tracked_data_.set_skipped(true);
        return;
      }

    clock.stop(" sequentiator: sequentiate ");
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return;
  }


  //*************************************************************
  void sequentiator::sequentiate_cluster(topology::cluster & cluster_) {
    //*************************************************************

    for (std::vector<topology::node>::iterator
          inode = cluster_.nodes_.begin();
        inode != cluster_.nodes_.end(); ++inode)
      {
        topology::node & a_node = *inode;

        DT_LOG_DEBUG(get_logging_priority(), "First node: " << inode->c().id());

        if (!good_first_node(a_node)) continue;

        make_new_sequence(a_node);

        if (late()) return;

        make_copy_sequence(a_node);
      }

    return;
  }


  //*************************************************************
  void sequentiator::make_new_sequence(topology::node & first_node){
    //*************************************************************

    if (late()) return;

    clock.start(" sequentiator: make new sequence ","cumulative");

    //  A node is added to the newsequence. It has the given cell but no other
    //  requirement. The free level is set to true.
    topology::sequence newsequence(first_node);
    newsequence.set_probmin(probmin);

    bool updated = true;

    while (updated)
      {
        updated = evolve(newsequence);
      }

    if (late()) return;

    NFAMILY++;
    NCOPY = 0;
    if (newsequence.nodes().size() != 2)
      {
        make_name(newsequence);

        sequences_.push_back(newsequence);
        DT_LOG_DEBUG(get_logging_priority(),"Finished track [" << sequences_.size()-1 << "] ");
        clean_up_sequences();
      }
    else
      {
        add_pair(newsequence);
      }

    clock.stop(" sequentiator: make new sequence ");

    return;
  }


  void sequentiator::make_name(topology::sequence & sequence_)
  {
    std::string number = mybhep::to_string(NFAMILY)+"_"+mybhep::to_string(NCOPY);
    std::string name = "track_"+number;
    sequence_.set_name(name);
    return;
  }

  void sequentiator::make_copy_sequence(topology::node & first_node)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");

    clock.start(" sequentiator: make copy sequence ","cumulative");

    size_t isequence;
    while (there_is_free_sequence_beginning_with(first_node.c(), &isequence))
      {
        if (late()) return;

        clock.start(" sequentiator: make copy sequence: part A ","cumulative");
        clock.start(" sequentiator: make copy sequence: part A: alpha ","cumulative");

        DT_LOG_DEBUG(get_logging_priority(), "Begin, with cell" << first_node.c().id() << ", parallel track "
                     << sequences_.size() << " to track" << isequence);

        clock.stop(" sequentiator: make copy sequence: part A: alpha ");
        clock.start(" sequentiator: copy to lfn ","cumulative");
        size_t ilink, ilfn;
        topology::sequence newcopy = sequences_[isequence].copy_to_last_free_node(&ilfn, &ilink);
        clock.stop(" sequentiator: copy to lfn ");

        clock.start(" sequentiator: make copy sequence: part A: beta ","cumulative");
        DT_LOG_DEBUG(get_logging_priority(), "Copied from sequence " << isequence);

        clock.stop(" sequentiator: make copy sequence: part A: beta ");
        clock.start(" sequentiator: make copy sequence: evolve ","cumulative");

        bool updated = true;
        while (updated)
          updated = evolve(newcopy);
        clock.stop(" sequentiator: make copy sequence: evolve ");

        if (late()) return;

        clock.stop(" sequentiator: make copy sequence: part A ");
        clock.start(" sequentiator: manage copy sequence ","cumulative");

        DT_LOG_DEBUG(get_logging_priority(), "Checking study case for sequence #" << isequence << " "
                     << "and node #" << first_node.c().id());;

        // not adding: case 1: new sequence did not evolve
        if (newcopy.nodes().size() == ilfn + 1)
          {
            DT_LOG_DEBUG(get_logging_priority(), "Not adding new sequence, since it couldn't evolve past lfn");
            clean_up_sequences();
          }
        else
          {
            // if copy has evolved past level ilfn, the link from node ilfn it used
            // is set to used in the original
            if( newcopy.nodes().size() > ilfn + 1 ){
              if( !sequences_[isequence].nodes().empty() ){
                clock.start(" sequentiator: get link index ","cumulative");
                size_t it1 = newcopy.get_link_index_of_cell(ilfn, newcopy.nodes()[ilfn + 1].c());
                clock.stop(" sequentiator: get link index ");
                DT_LOG_DEBUG(get_logging_priority(), "Setting as used original node " << ilfn << " cc " << it1);
                if( ilfn == 0 )
                  sequences_[isequence].nodes_[ilfn].cc_[it1].set_all_used();
                else
                  sequences_[isequence].nodes_[ilfn].ccc_[it1].set_all_used();
              }
              /*
                if( sequences_[isequence].nodes().size() > 1 && ilfn > 0){
                clock.start(" sequentiator: get link index ","cumulative");
                size_t it2 = newcopy.get_link_index_of_cell(1, newcopy.nodes()[2].c());
                clock.stop(" sequentiator: get link index ");
                m.message(" setting as used original node 1  ccc ", it2, mybhep::VVERBOSE);
                sequences_[isequence].nodes_[1].ccc_[it2].set_all_used();
                }
              */
              clock.start(" sequentiator: set free level ","cumulative");
              sequences_[isequence].set_free_level();
              clock.stop(" sequentiator: set free level ");
            }

          // not adding: case 2: new sequence contained
          if (newcopy.contained(sequences_[isequence]) && !newcopy.Free())  // new copy is contained in original
            {
              DT_LOG_DEBUG(get_logging_priority(), "Not adding new sequence, contained in " << isequence << " from which it was copied");
              clean_up_sequences();
            }
          else
            { // adding: case 3
              if (sequences_[isequence].contained(newcopy))
                { // original is contained in new copy
                  for (size_t k=0; k<ilfn; k++)
                    {
                      newcopy.nodes_[k].set_free( sequences_[isequence].nodes()[k].free());

                      for (std::vector<topology::cell_couplet>::iterator icc = sequences_[isequence].nodes_[k].cc_.begin();
                          icc != sequences_[isequence].nodes_[k].cc_.end(); ++icc)
                        {
                          newcopy.nodes_[k].cc_[icc - sequences_[isequence].nodes_[k].cc_.begin()].set_free( icc->free());
                          newcopy.nodes_[k].cc_[icc - sequences_[isequence].nodes_[k].cc_.begin()].set_begun( icc->begun());

                          for(std::vector<topology::line>::iterator itang = sequences_[isequence].nodes_[k].cc_[icc - sequences_[isequence].nodes_[k].cc_.begin()].tangents_.begin(); itang != sequences_[isequence].nodes_[k].cc_[icc - sequences_[isequence].nodes_[k].cc_.begin()].tangents_.end(); ++itang)
                            newcopy.nodes_[k].cc_[icc - sequences_[isequence].nodes_[k].cc_.begin()].tangents_[itang - sequences_[isequence].nodes_[k].cc_[icc - sequences_[isequence].nodes_[k].cc_.begin()].tangents_.begin()].set_used(itang->used() );

                        }

                      for (std::vector<topology::cell_triplet>::iterator iccc = sequences_[isequence].nodes_[k].ccc_.begin();
                           iccc != sequences_[isequence].nodes_[k].ccc_.end(); ++iccc)
                        {
                          newcopy.nodes_[k].ccc_[iccc - sequences_[isequence].nodes_[k].ccc_.begin()].set_free( iccc->free());
                          newcopy.nodes_[k].ccc_[iccc - sequences_[isequence].nodes_[k].ccc_.begin()].set_begun( iccc->begun());

                          for (std::vector<topology::joint>::iterator ijoint = sequences_[isequence].nodes_[k].ccc_[iccc - sequences_[isequence].nodes_[k].ccc_.begin()].joints_.begin(); ijoint != sequences_[isequence].nodes_[k].ccc_[iccc - sequences_[isequence].nodes_[k].ccc_.begin()].joints_.end(); ++ijoint)
                            newcopy.nodes_[k].ccc_[iccc - sequences_[isequence].nodes_[k].ccc_.begin()].joints_[ijoint - sequences_[isequence].nodes_[k].ccc_[iccc - sequences_[isequence].nodes_[k].ccc_.begin()].joints_.begin()].set_used(ijoint->used() );

                        }
                    }

                  if (ilfn < 2)
                    {
                      for (std::vector<topology::cell_couplet>::iterator icc = sequences_[isequence].nodes_[ilfn].cc_.begin();
                           icc != sequences_[isequence].nodes_[ilfn].cc_.end(); ++icc)
                        if ((size_t)(icc - sequences_[isequence].nodes_[ilfn].cc_.begin()) != ilink)
                          {
                            newcopy.nodes_[ilfn].cc_[icc - sequences_[isequence].nodes_[ilfn].cc_.begin()].set_free( icc->free());
                            newcopy.nodes_[ilfn].cc_[icc - sequences_[isequence].nodes_[ilfn].cc_.begin()].set_begun( icc->begun());

                            for (std::vector<topology::line>::iterator itang = sequences_[isequence].nodes_[ilfn].cc_[icc - sequences_[isequence].nodes_[ilfn].cc_.begin()].tangents_.begin(); itang !=sequences_[isequence].nodes_[ilfn].cc_[icc - sequences_[isequence].nodes_[ilfn].cc_.begin()].tangents_.end(); ++itang)
                              newcopy.nodes_[ilfn].cc_[icc - sequences_[isequence].nodes_[ilfn].cc_.begin()].tangents_[itang - sequences_[isequence].nodes_[ilfn].cc_[icc - sequences_[isequence].nodes_[ilfn].cc_.begin()].tangents_.begin()].set_used(itang->used() );
                          }
                    }
                  else
                    {
                      for (std::vector<topology::cell_triplet>::iterator iccc = sequences_[isequence].nodes_[ilfn].ccc_.begin();
                          iccc != sequences_[isequence].nodes_[ilfn].ccc_.end(); ++iccc)
                        if ((size_t)(iccc - sequences_[isequence].nodes_[ilfn].ccc_.begin()) != ilink )
                          {
                            newcopy.nodes_[ilfn].ccc_[iccc - sequences_[isequence].nodes_[ilfn].ccc_.begin()].set_free( iccc->free());
                            newcopy.nodes_[ilfn].ccc_[iccc - sequences_[isequence].nodes_[ilfn].ccc_.begin()].set_begun( iccc->begun());

                            for (std::vector<topology::joint>::iterator ijoint = sequences_[isequence].nodes_[ilfn].ccc_[iccc - sequences_[isequence].nodes_[ilfn].ccc_.begin()].joints_.begin(); ijoint !=sequences_[isequence].nodes_[ilfn].ccc_[iccc - sequences_[isequence].nodes_[ilfn].ccc_.begin()].joints_.end(); ++ijoint)
                              newcopy.nodes_[ilfn].ccc_[iccc - sequences_[isequence].nodes_[ilfn].ccc_.begin()].joints_[ijoint - sequences_[isequence].nodes_[ilfn].ccc_[iccc - sequences_[isequence].nodes_[ilfn].ccc_.begin()].joints_.begin()].set_used(ijoint->used() );
                          }
                    }

                  clock.start(" sequentiator: set free level ","cumulative");
                  newcopy.set_free_level();
                  clock.stop(" sequentiator: set free level ");

                  sequences_.erase(sequences_.begin()+isequence);
                  DT_LOG_DEBUG(get_logging_priority(), "Erased original sequence " << isequence << " contained in sequence "
                               << sequences_.size()+1 << " which was copied from it");
                  clean_up_sequences();
                }

              NCOPY++;
              if (newcopy.nodes().size() != 2)
                {
                  make_name(newcopy);
                  sequences_.push_back( newcopy );
                  DT_LOG_DEBUG(get_logging_priority(), "Finished track [" << sequences_.size()-1 << "]");
                  clean_up_sequences();
                }
              else
                {
                  add_pair(newcopy);
                }
            }// end of case 3
          }

        clock.stop(" sequentiator: manage copy sequence ");
      }

    NCOPY = 0;

    clock.stop(" sequentiator: make copy sequence ");

    return;
  }

  //*************************************************************
  bool sequentiator::late(void){
    //*************************************************************

    if( clock.read(" sequentiator: sequentiation ") >= MaxTime ){

      DT_LOG_DEBUG(get_logging_priority(), "Execution time " << clock.read(" sequentiator: sequentiation ") << " ms greater than MaxTime" << MaxTime << " quitting! ");

      //      clock.stop_all();

      /*
        if( !evt.find_property("MaxTime") )
        evt.add_property("MaxTime","1");
        else
        }
      */

      return true;
    }

    return false;

  }

  //*************************************************************
  bool sequentiator::evolve(topology::sequence & sequence){
    //*************************************************************

    DT_LOG_TRACE(get_logging_priority(), "Entering...");

    if (late()) return false;

    clock.start(" sequentiator: evolve ","cumulative");

    clock.start(" sequentiator: evolve: part A ","cumulative");

    const size_t sequence_size = sequence.nodes().size();

    DT_LOG_DEBUG(get_logging_priority(), "Sequence size = " << sequence_size);

    // protection
    if( sequence_size < 1 )
      {
        DT_LOG_DEBUG(get_logging_priority(), "Problem: sequence has length " << sequence_size << "... stop evolving");
        clock.stop(" sequentiator: evolve: part A ");
        clock.stop(" sequentiator: evolve ");
        return false;
      }

    if( sequence_size == 3 ){
      clock.start(" sequentiator: get link index ","cumulative");
      size_t it1 = sequence.get_link_index_of_cell(0, sequence.nodes()[1].c());
      if( it1 >= sequence.nodes_[0].cc_.size() ){
        DT_LOG_DEBUG(get_logging_priority(), "Problem: it1 " << it1 << " nodes size " << sequence.nodes_.size() << " cc size " << sequence.nodes_[0].cc_.size());
        clock.stop(" sequentiator: evolve: part A ");
        clock.stop(" sequentiator: evolve ");
        return false;
      }
      sequence.nodes_[0].cc_[it1].set_all_used();

      clock.stop(" sequentiator: get link index ");
    }

    clock.stop(" sequentiator: evolve: part A ");
    clock.start(" sequentiator: evolve: part B ","cumulative");

    // check if there is a possible link
    size_t ilink;
    topology::experimental_point newp;
    clock.start(" sequentiator: pick new cell ","cumulative");
    bool there_is_link = sequence.pick_new_cell(&ilink, &newp, *local_cluster_);
    clock.stop(" sequentiator: pick new cell ");

    DT_LOG_DEBUG(get_logging_priority(), "Is there a link = " << there_is_link);

    if( !there_is_link ){
      DT_LOG_DEBUG(get_logging_priority(), "No links could be added... stop evolving");
      clock.start(" sequentiator: evolve: part B: set free level ","cumulative");
      clock.start(" sequentiator: set free level ","cumulative");
      sequence.set_free_level();
      clock.stop(" sequentiator: set free level ");
      clock.stop(" sequentiator: evolve: part B: set free level ");
      clock.stop(" sequentiator: evolve: part B ");
      clock.stop(" sequentiator: evolve ");

      if (sequence.nodes().size() == 1){
        topology::experimental_point ep(sequence.nodes_[0].c().ep());
        ep.set_ex(sequence.nodes_[0].c().r().value());
        ep.set_ez(sequence.nodes_[0].c().r().value());
        sequence.nodes_[0].set_ep(ep);
      }


      return false;
    }

    topology::cell newcell = sequence.last_node().links()[ilink];
    clock.start(" sequentiator: evolve: part B: noc ","cumulative");
    topology::node newnode = local_cluster_->node_of_cell(newcell);
    //  topology::node newnode = local_cluster_->nodes()[local_cluster_->node_index_of_cell(newcell)];
    clock.stop(" sequentiator: evolve: part B: noc ");
    newnode.set_free(false); // standard initialization

    clock.stop(" sequentiator: evolve: part B ");
    clock.start(" sequentiator: evolve: part C ","cumulative");

    if( sequence_size == 1 )
      {
        // since it's the 2nd cell, only the four
        // standard alternatives are available (which will be treated out of "Evolve"), so
        // this link has no freedom left

        sequence.nodes_.push_back( newnode );
        clock.start(" sequentiator: set free level ","cumulative");
        sequence.set_free_level();
        clock.stop(" sequentiator: set free level ");

        clock.stop(" sequentiator: evolve: part C ");
        clock.stop(" sequentiator: evolve ");
        return true;
      }


    newnode.set_ep(newp);

    sequence.nodes_.push_back( newnode );
    DT_LOG_DEBUG(get_logging_priority(), "Points have been added");

    clock.start(" sequentiator: set free level ","cumulative");
    sequence.set_free_level();
    clock.stop(" sequentiator: set free level ");

    clock.stop(" sequentiator: evolve: part C ");
    clock.stop(" sequentiator: evolve ");
    return true;
  }



  //*************************************************************
  bool sequentiator::good_first_node(topology::node & node_) {
    //*************************************************************

    clock.start(" sequentiator: good first node ", "cumulative");

    const std::string type = node_.topological_type();

    // check that node is not in the middle of a cell_triplet
    if( type != "VERTEX" &&
        type != "MULTI_VERTEX" &&
        type != "ISOLATED" ){
      // clock.stop(" sequentiator: good first node ");
      DT_LOG_DEBUG(get_logging_priority(), "Not a good first node: type " << type);
      return false;
    }


    std::vector<size_t> done_connections;
    size_t connection_node;
    std::vector<size_t>::iterator fid;
    // check that node has never been added to a sequence
    for(std::vector<topology::sequence>::const_iterator iseq=sequences_.begin(); iseq!=sequences_.end(); ++iseq)
      {
        if( iseq->has_cell(node_.c()) ){
          if( type == "VERTEX" ){

            DT_LOG_DEBUG(get_logging_priority(), "Already used as vertex in sequence " << iseq - sequences_.begin());

            // clock.stop(" sequentiator: good first node ");
            return false;
          }
          else{ // multi-vertex
            if( iseq->nodes_.size() > 1 ){
              if( iseq->nodes_[0].c().id() == node_.c().id() )
                connection_node = 1;
              else if( iseq->last_node().c().id() == node_.c().id() ){
                connection_node = iseq->nodes_.size() - 2;
              }
              else{
                DT_LOG_DEBUG(get_logging_priority(), "Problem: multi-vertex " << node_.c().id() << " belongs to sequence "
                             << iseq->name() << " but not as first or last cell");
                continue;
              }
              // add to done_connections cell ids of those cells
              // that have already been connected to NODE in other sequences
              fid = std::find(done_connections.begin(),
                              done_connections.end(),
                              iseq->nodes_[connection_node].c().id());

              if( fid == done_connections.end())
                done_connections.push_back(iseq->nodes_[connection_node].c().id());
            }
          }
        }
      }

    size_t cc_index;
    if( type == "MULTI_VERTEX" ){
      for(size_t i=0; i<done_connections.size(); i++){
        cc_index = 0;
        if( !node_.has_couplet(done_connections[i],  &cc_index) ) {
          DT_LOG_DEBUG(get_logging_priority(), "Problem: multi-vertex " << node_.c().id()
                       << " should link to cell " << done_connections[i] << " but has not such couplet");
        } else{
          DT_LOG_DEBUG(get_logging_priority(), "Multi-vertex " << node_.c().id()
                       << " has already been added to a sequence connecting to cell " << done_connections[i]
                       << " so couplet " << cc_index << " will be erased");
          //node_.cc_.erase(node_.cc_.begin() + cc_index);
          node_.remove_couplet(cc_index);
        }
      }
    }

    clock.stop(" sequentiator: good first node ");
    return true;


  }


  //*************************************************************
  void sequentiator::make_families() {
    //*************************************************************

    clock.start(" sequentiator: make families ", "cumulative");

    families_.clear();

    bool found, added;
    size_t ifam;
    std::vector<size_t> Fam;
    for (std::vector<topology::sequence>::const_iterator iseq = sequences_.begin();
         iseq != sequences_.end(); ++iseq){

      const size_t ipart = iseq - sequences_.begin();
      const size_t this_family = mybhep::int_from_string(iseq->family());

      found = false;
      for(size_t i=0; i<families_.size(); i++)
        if( families_[i][0] == this_family )
          {
            found = true;
            ifam = i;
            break;
          }

      if( found )
        {
          added = false;
          for(size_t j=0; j<families_[ifam].size(); j++)
            {
              if(j == 0) continue;
              if( families_[ifam][j] == ipart )
                {
                  added = true;
                  break;
                }
            }
          if( !added )
            families_[ifam].push_back( ipart );
        }
      else
        {
          Fam.clear();
          Fam.push_back(this_family);

          Fam.push_back(ipart);
          families_.push_back( Fam );
        }
    }

    clock.stop(" sequentiator: make families ");

    return;

  }

  //*************************************************************
  bool sequentiator::make_scenarios(topology::tracked_data &td){
    //*************************************************************

    clock.start(" sequentiator: make scenarios ", "cumulative");

    for(std::vector<topology::sequence>::iterator iseq=sequences_.begin(); iseq!=sequences_.end(); ++iseq){
      if( late() ){
        td.set_skipped(true);
        return false;
      }


      DT_LOG_DEBUG(get_logging_priority(), "Begin scenario with sequence " << iseq->name());

      topology::scenario sc;
      sc.set_probmin(probmin);
      sc.sequences_.push_back(*iseq);
      sc.calculate_n_free_families(td.get_cells(), td.get_calos());
      sc.calculate_n_overlaps(td.get_cells(), td.get_calos());
      sc.calculate_chi2();

      size_t jmin = 0, nfree = 0, noverlaps = 0;
      double Chi2 = 0.0;
      int ndof = 0;
      while( can_add_family(sc, &jmin, &nfree, &Chi2, &noverlaps, &ndof, td) )
        {
          DT_LOG_DEBUG(get_logging_priority(), "Best sequence to add is " << jmin);
          DT_LOG_DEBUG(get_logging_priority(), "nfree " << nfree << " noverls " << noverlaps << " Chi2 " << Chi2);
          sc.sequences_.push_back(sequences_[jmin]);
          sc.set_n_free_families(nfree);
          sc.set_helix_chi2(Chi2);
          sc.set_ndof(ndof);
          sc.set_n_overlaps(noverlaps);
        }

      scenarios_.push_back(sc);


    }

    if( late() ){
      td.set_skipped(true);
      return false;
    }


    direct_scenarios_out_of_foil();

    if( scenarios_.size() > 0 ){

      size_t index_tmp = pick_best_scenario();

      DT_LOG_DEBUG(get_logging_priority(), "Made scenario");

      td.scenarios_.push_back(scenarios_[index_tmp]);

      clock.stop(" sequentiator: make scenarios ");
      return true;

    }

    DT_LOG_DEBUG(get_logging_priority(), "Not made scenario");
    clock.stop(" sequentiator: make scenarios ");

    return false;


  }


  //*************************************************************
  size_t sequentiator::pick_best_scenario(){
    //*************************************************************

    DT_LOG_TRACE(get_logging_priority(), "Entering...");

    size_t index = 0;

    for(std::vector<topology::scenario>::iterator sc=scenarios_.begin(); sc!=scenarios_.end(); ++sc){
      DT_LOG_DEBUG(get_logging_priority(), " ...scenario " << sc - scenarios_.begin()
                   << " nff " << sc->n_free_families() << " noverls " << sc->n_overlaps()
                   << " common vertexes " << sc->n_of_common_vertexes(2.*CellDistance)
                   << " n ends on wire " << sc->n_of_ends_on_wire() << " chi2 " << sc->helix_chi2()
                   << " prob " << sc->helix_Prob());

      if( sc->better_scenario_than( scenarios_[index] , 2.*CellDistance ) )
        {
          index = sc - scenarios_.begin();
        }
    }

    DT_LOG_DEBUG(get_logging_priority(), "Best scenario is " << index);

    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return index;
  }

  //*************************************************************
  bool sequentiator::can_add_family(topology::scenario &sc, size_t* jmin, size_t* nfree, double* Chi2, size_t* noverlaps, int* ndof, topology::tracked_data &td) {
    //*************************************************************

    if( late() )
      return false;

    clock.start(" sequentiator: can add family ", "cumulative");

    bool ok = false;

    if( sc.n_free_families() == 0 ){
      clock.stop(" sequentiator: can add family ");
      return false;
    }

#if 0
    clock.start(" sequentiator: copy logic scenario ", "cumulative");
    topology::logic_scenario tmpmin = topology::logic_scenario(sc);
    topology::logic_scenario tmpsave = tmpmin;
    clock.stop(" sequentiator: copy logic scenario ");
    topology::logic_scenario tmp;
#else
    clock.start(" sequentiator: copy scenario ", "cumulative");
    topology::scenario tmpmin = sc;
    clock.stop(" sequentiator: copy scenario ");
    topology::scenario tmp = sc;
#endif

    std::map<std::string,int> scnames;
    for(std::vector<topology::sequence>::iterator iseq=sc.sequences_.begin(); iseq!=sc.sequences_.end(); ++iseq)
      scnames[iseq->name()]=iseq-sc.sequences_.begin();

    for(std::vector<topology::sequence>::iterator jseq=sequences_.begin(); jseq!=sequences_.end(); ++jseq)
      {

        if( scnames.count(jseq->name()) ) continue;

#if 0
        clock.start(" sequentiator: copy logic scenario ", "cumulative");
        tmp = tmpsave;
        clock.stop(" sequentiator: copy logic scenario ");
        clock.start(" sequentiator: copy logic sequence ", "cumulative");
        tmp.sequences_.push_back(topology::logic_sequence(*jseq));
        clock.stop(" sequentiator: copy logic sequence ");
#else
        clock.start(" sequentiator: copy scenario ", "cumulative");
        tmp = sc;
        clock.stop(" sequentiator: copy scenario ");
        clock.start(" sequentiator: copy sequence ", "cumulative");
        tmp.sequences_.push_back(*jseq);
        clock.stop(" sequentiator: copy sequence ");
#endif

        clock.start(" sequentiator: calculate scenario ", "cumulative");
        tmp.calculate_n_free_families(td.get_cells(), td.get_calos());
        tmp.calculate_n_overlaps(td.get_cells(), td.get_calos());
        tmp.calculate_chi2();
        clock.stop(" sequentiator: calculate scenario ");

        DT_LOG_DEBUG(get_logging_priority(), " ...try to add sequence " << jseq->name());
        DT_LOG_DEBUG(get_logging_priority(), " ...nfree " << tmp.n_free_families() << " noverls " << tmp.n_overlaps()
                     << " chi2 " << tmp.helix_chi2() << " prob " << tmp.helix_Prob());

        clock.start(" sequentiator: better scenario ", "cumulative");
        if( tmp.better_scenario_than(tmpmin , 2.*CellDistance ) )
          {
            *jmin = jseq - sequences_.begin();
            *nfree = tmp.n_free_families();
            *noverlaps = tmp.n_overlaps();
            *Chi2 = tmp.helix_chi2();
            *ndof = tmp.ndof();
            tmpmin = tmp;
            ok = true;
          }
        clock.stop(" sequentiator: better scenario ");

      }

    clock.stop(" sequentiator: can add family ");
    return ok;

  }

  //*************************************************************
  topology::plane sequentiator::get_foil_plane(){
    //*************************************************************

    topology::experimental_point center(0., 0., 0., 0., 0., 0.);

    topology::experimental_vector norm(0.,0.,1.,0.,0.,0.);

    topology::experimental_vector sizes(xsize, ysize, 0.,
                                        0., 0., 0.);

    topology::plane pl(center, sizes, norm);
    pl.set_probmin(probmin);

    std::string the_type="Nemo3";
    if( SuperNemo )
      the_type="SuperNEMO";

    pl.set_type(the_type);

    return pl;


  }

  //*************************************************************
  void sequentiator::refine_sequences_near_walls(std::vector<topology::calorimeter_hit> & calos){
    //*************************************************************

    for( std::vector<topology::sequence>::iterator iseq = sequences_.begin(); iseq!=sequences_.end(); ++iseq){

      if( iseq->nodes_.size() < 3 ) continue;
      if( gap_number( iseq->second_last_node().c() ) == 0 &&
          iseq->phi_kink(iseq->nodes_.size()-2)*180./M_PI > 45 &&
          belongs_to_other_family(iseq->last_node().c(), &(*iseq)) ){
        DT_LOG_DEBUG(get_logging_priority(), "Removing last node " << iseq->last_node().c()
                     << " near foil of sequence " << iseq->name() << " (it belongs to other family and makes large kink)");
        iseq->remove_last_node();
      }

      if( iseq->nodes_.size() < 3 ) continue;
      if( gap_number( iseq->nodes_[1].c() ) == 0 &&
          iseq->phi_kink(1)*180./M_PI > 45 &&
          belongs_to_other_family(iseq->nodes_[0].c(), &(*iseq)) ){
        DT_LOG_DEBUG(get_logging_priority(), "Removing 1st node " << iseq->last_node().c() << " near foil of sequence "
                     << iseq->name() << " (it belongs to other family and makes large kink)");
        iseq->remove_first_node();
      }


      if( iseq->nodes_.size() < 3 ) continue;
      for(std::vector<topology::calorimeter_hit>::iterator ic=calos.begin(); ic != calos.end(); ++ic){
        if( near(iseq->nodes_[1].c(), *ic) &&
            iseq->phi_kink(1)*180./M_PI > 45 &&
            belongs_to_other_family(iseq->nodes_[0].c(), &(*iseq)) ){
          DT_LOG_DEBUG(get_logging_priority(), "Removing 1st node " << iseq->last_node().c() << " near calo of sequence "
                       << iseq->name() << " (it belongs to other family and makes large kink)");
          iseq->remove_first_node();
          break;
        }
      }

      if( iseq->nodes_.size() < 3 ) continue;
      for(std::vector<topology::calorimeter_hit>::iterator ic=calos.begin(); ic != calos.end(); ++ic){
        if( near(iseq->second_last_node().c(), *ic) &&
            iseq->phi_kink(iseq->nodes_.size()-2)*180./M_PI > 45 &&
            belongs_to_other_family(iseq->last_node().c(), &(*iseq)) ){
          DT_LOG_DEBUG(get_logging_priority(), "Removing last node " << iseq->last_node().c() << " near calo of sequence " <<
                       iseq->name() <<  " (it belongs to other family and makes large kink)");
          iseq->remove_last_node();
          break;
        }
      }
    }

    return;
  }


  //*************************************************************
  bool sequentiator::belongs_to_other_family(topology::cell c, topology::sequence *iseq){
    //*************************************************************

    for( std::vector<topology::sequence>::iterator jseq = sequences_.begin(); jseq!=sequences_.end(); ++jseq){
      if( iseq->same_families(*jseq) ) continue;
      if( !jseq->has_cell(c) ) continue;
      return true;
    }

    return false;

  }

  //*************************************************************
  void sequentiator::interpret_physics(std::vector<topology::calorimeter_hit> & calos){
    //*************************************************************

    clock.start(" sequentiator: interpret physics ", "cumulative");

    DT_LOG_DEBUG(get_logging_priority(), "Interpreting physics of " << sequences_.size() << " sequences with "
                 << calos.size() << " calorimeter hits");

    double helix_min_from_end = mybhep::default_min;
    size_t ihelix_min_from_end = mybhep::default_integer;
    double tangent_min_from_end = mybhep::default_min;
    size_t itangent_min_from_end = mybhep::default_integer;

    double helix_min_from_begin = mybhep::default_min;
    size_t ihelix_min_from_begin = mybhep::default_integer;
    double tangent_min_from_begin = mybhep::default_min;
    size_t itangent_min_from_begin = mybhep::default_integer;

    topology::experimental_point helix_extrapolation_from_end, helix_extrapolation_local_from_end;
    bool helix_found_from_end = false;
    topology::experimental_point helix_extrapolation_from_begin, helix_extrapolation_local_from_begin;
    bool helix_found_from_begin = false;

    topology::experimental_point tangent_extrapolation_from_end, tangent_extrapolation_local_from_end;
    bool tangent_found_from_end = false;
    topology::experimental_point tangent_extrapolation_from_begin, tangent_extrapolation_local_from_begin;
    bool tangent_found_from_begin = false;

    double dist_from_end, dist_from_begin;
    std::vector<topology::sequence>::iterator iseq = sequences_.begin();
    while( iseq != sequences_.end() )
      {
        DT_LOG_DEBUG(get_logging_priority(), " ... interpreting physics of sequence " << iseq->name());

        if( iseq->nodes().size() <= 2 ){
          ++iseq;
          continue;
        }

      // by default, conserve_clustering_from_removal_of_whole_cluster = false and conserve_clustering_from_reordering = true
        if( !iseq->calculate_helix(Ratio) && !iseq->has_kink() ){
          size_t index = iseq - sequences_.begin();
          DT_LOG_DEBUG(get_logging_priority(), "Erased sequence " << index << " not a good helix");
          sequences_.erase(iseq);
          iseq = sequences_.begin() + index;
          if( index + 1 >= sequences_.size() )
            break;
          continue;
        }
        iseq->calculate_charge();

        // match to calorimeter
        if (!calos.empty ())
          {

            DT_LOG_DEBUG(get_logging_priority(), "Extrapolate decay vertex with " << calos.size() << " calo hits");

            helix_min_from_end = mybhep::default_min;
            ihelix_min_from_end = mybhep::default_integer;
            tangent_min_from_end = mybhep::default_min;
            itangent_min_from_end = mybhep::default_integer;

            helix_min_from_begin = mybhep::default_min;
            ihelix_min_from_begin = mybhep::default_integer;
            tangent_min_from_begin = mybhep::default_min;
            itangent_min_from_begin = mybhep::default_integer;

            helix_found_from_end = false;
            helix_found_from_begin = false;

            tangent_found_from_end = false;
            tangent_found_from_begin = false;

            for(std::vector<topology::calorimeter_hit>::iterator ic=calos.begin(); ic != calos.end(); ++ic){

              DT_LOG_DEBUG(get_logging_priority(), "Trying to extrapolate to calo hit " << ic - calos.begin()
                           << " id " << ic->id() << " on view " << ic->pl_.view() << " energy " << ic->e().value());

              if( !near(iseq->last_node().c(), *ic) ){
                DT_LOG_DEBUG(get_logging_priority(), "End is not near");
              } else {
                if( !iseq->intersect_plane_from_end(ic->pl(), &helix_extrapolation_local_from_end) ){
                  DT_LOG_DEBUG(get_logging_priority(), "No helix intersection from end");
                } else {
                  dist_from_end = helix_extrapolation_local_from_end.distance(ic->pl_.face()).value();
                  if( dist_from_end < helix_min_from_end ){
                    helix_min_from_end = dist_from_end;
                    ihelix_min_from_end = ic->id();
                    helix_extrapolation_from_end = helix_extrapolation_local_from_end;
                    helix_found_from_end = true;
                    DT_LOG_DEBUG(get_logging_priority(), "New helix intersection from end with minimum distance "
                                 << dist_from_end << " position: "
                                 << helix_extrapolation_from_end.x().value() << ","
                                 << helix_extrapolation_from_end.y().value() << ","
                                 << helix_extrapolation_from_end.z().value());
                  }
                }


                if( !iseq->intersect_plane_with_tangent_from_end(ic->pl(), &tangent_extrapolation_local_from_end) ){
                  DT_LOG_DEBUG(get_logging_priority(), "No tangent intersection from end");
                }
                else{

                  dist_from_end = tangent_extrapolation_local_from_end.distance(ic->pl_.face()).value();
                  if( dist_from_end < tangent_min_from_end ){
                    tangent_min_from_end = dist_from_end;
                    itangent_min_from_end = ic->id();
                    tangent_extrapolation_from_end = tangent_extrapolation_local_from_end;
                    tangent_found_from_end = true;
                    DT_LOG_DEBUG(get_logging_priority(), "New tangent intersection from end with minimum distance "
                                 << dist_from_end << " position: "
                                 << tangent_extrapolation_from_end.x().value() << ","
                                 << tangent_extrapolation_from_end.y().value() << ","
                                 << tangent_extrapolation_from_end.z().value());
                  }
                }
              }

              if( !near(iseq->nodes_[0].c(), *ic) ){
                DT_LOG_DEBUG(get_logging_priority(), "Beginning is not near");
              }else if( ihelix_min_from_end == ic->id() || itangent_min_from_end == ic->id() ){
                DT_LOG_DEBUG(get_logging_priority(), "Beginning is near, but end was already extrapolated to same calo");
              }else{
                if( !iseq->intersect_plane_from_begin(ic->pl(), &helix_extrapolation_local_from_begin) ){
                  DT_LOG_DEBUG(get_logging_priority(), "No helix intersection from beginning");
                }
                else{

                  dist_from_begin = helix_extrapolation_local_from_begin.distance(ic->pl_.face()).value();
                  if( dist_from_begin < helix_min_from_begin ){
                    helix_min_from_begin = dist_from_begin;
                    ihelix_min_from_begin = ic->id();
                    helix_extrapolation_from_begin = helix_extrapolation_local_from_begin;
                    helix_found_from_begin = true;
                    DT_LOG_DEBUG(get_logging_priority(), "New helix intersection from beginning with minimum distance "
                                 << dist_from_begin << " position: "
                                 << helix_extrapolation_from_begin.x().value() << ","
                                 << helix_extrapolation_from_begin.y().value() << ","
                                 << helix_extrapolation_from_begin.z().value());
                  }
                }


                if( !iseq->intersect_plane_with_tangent_from_begin(ic->pl(), &tangent_extrapolation_local_from_begin) ){
                  DT_LOG_DEBUG(get_logging_priority(), "No tangent intersection from beginning");
                }
                else{
                  dist_from_begin = tangent_extrapolation_local_from_begin.distance(ic->pl_.face()).value();
                  if( dist_from_begin < tangent_min_from_begin ){
                    tangent_min_from_begin = dist_from_begin;
                    itangent_min_from_begin = ic->id();
                    tangent_extrapolation_from_begin = tangent_extrapolation_local_from_begin;
                    tangent_found_from_begin = true;
                    DT_LOG_DEBUG(get_logging_priority(), "New tangent intersection from beginning with minimum distance "
                                 << dist_from_begin << " position: "
                                 << tangent_extrapolation_from_begin.x().value()
                                 << tangent_extrapolation_from_begin.y().value()
                                 << tangent_extrapolation_from_begin.z().value());
                  }
                }
              }


            } // finish loop on calos

            if( helix_found_from_begin ){
              if( ihelix_min_from_begin >= calos.size() ){
                DT_LOG_DEBUG(get_logging_priority(), "Problem: calo hit of id "
                             << ihelix_min_from_begin << " but n of calo hits is " << calos.size());
              }
              else{
                DT_LOG_DEBUG(get_logging_priority(), "Track extrapolated by helix to calo " << ihelix_min_from_begin);
                iseq->set_helix_vertex(helix_extrapolation_from_begin, "calo", ihelix_min_from_begin);
              }
            }

            if( tangent_found_from_begin ){
              if( itangent_min_from_begin >= calos.size() ){
                DT_LOG_DEBUG(get_logging_priority(), "Problem: tangent calo hit of id " << itangent_min_from_begin
                             << " but n of calo hits is " << calos.size());
              }
              else{
                DT_LOG_DEBUG(get_logging_priority(), "Track extrapolated by tangent to calo " << itangent_min_from_begin);
                iseq->set_tangent_vertex(tangent_extrapolation_from_begin, "calo", itangent_min_from_begin);
              }
            }

            if( helix_found_from_end ){
              if( ihelix_min_from_end >= calos.size() ){
                DT_LOG_DEBUG(get_logging_priority(), "Problem: calo hit of id " << ihelix_min_from_end
                             << " but n of calo hits is " << calos.size());
              }
              else{
                DT_LOG_DEBUG(get_logging_priority(), "Track extrapolated by helix to calo " << ihelix_min_from_end);
                iseq->set_decay_helix_vertex(helix_extrapolation_from_end, "calo", ihelix_min_from_end);
              }
            }

            if( tangent_found_from_end ){
              if( itangent_min_from_end >= calos.size() ){
                DT_LOG_DEBUG(get_logging_priority(), "Problem: tangent calo hit of id " << itangent_min_from_end
                             << " but n of calo hits is " << calos.size());
              }
              else{
                DT_LOG_DEBUG(get_logging_priority(), "Track extrapolated by tangent to calo " << itangent_min_from_end);
                iseq->set_decay_tangent_vertex(tangent_extrapolation_from_end, "calo", itangent_min_from_end);
              }
            }

          }

        // match to foil
        if( !iseq->nodes_.empty() ){
          DT_LOG_DEBUG(get_logging_priority(), "Extrapolate vertex on foil:");

          if( gap_number(iseq->last_node().c() ) != 0 ){
            DT_LOG_DEBUG(get_logging_priority(), "End not near");
          }else{
            if( SuperNemo ){

              if( !iseq->intersect_plane_from_end(get_foil_plane(), &helix_extrapolation_from_end) ){
                DT_LOG_DEBUG(get_logging_priority(), "No helix intersection from end");
              }
              else{
                iseq->set_decay_helix_vertex(helix_extrapolation_from_end, "foil");

              }

              if( !iseq->intersect_plane_with_tangent_from_end(get_foil_plane(), &tangent_extrapolation_from_end) ){
                DT_LOG_DEBUG(get_logging_priority(), "No tangent intersection from end");
              }
              else
                iseq->set_decay_tangent_vertex(tangent_extrapolation_from_end, "foil");

            }
          }

          if( gap_number(iseq->nodes_[0].c() ) != 0 ){
            DT_LOG_DEBUG(get_logging_priority(), "Beginning not near");
          }else{
            if( SuperNemo ){

              if( !iseq->intersect_plane_from_begin(get_foil_plane(), &helix_extrapolation_from_begin) ){
                DT_LOG_DEBUG(get_logging_priority(), "No helix intersection from beginning");
              }
              else{
                iseq->set_helix_vertex(helix_extrapolation_from_begin, "foil");

              }

              if( !iseq->intersect_plane_with_tangent_from_begin(get_foil_plane(), &tangent_extrapolation_from_begin) ){
                DT_LOG_DEBUG(get_logging_priority(), "No tangent intersection from beginning");
              }
              else
                iseq->set_tangent_vertex(tangent_extrapolation_from_begin, "foil");
            }
          }
        }

	if( SuperNemo )
	  iseq->calculate_momentum(bfield, SuperNemo, get_foil_plane().center().z().value());
        iseq->calculate_length();

        if (get_logging_priority() >= datatools::logger::PRIO_TRACE) {
          DT_LOG_TRACE(get_logging_priority(), "Sequence " << iseq - sequences_.begin() << " has: ");
          DT_LOG_TRACE(get_logging_priority(), "Center"); iseq->center().dump();
          DT_LOG_TRACE(get_logging_priority(), "Radius"); iseq->radius().dump();
          DT_LOG_TRACE(get_logging_priority(), "Pitch"); iseq->pitch().dump();
          DT_LOG_TRACE(get_logging_priority(), "Momentum"); iseq->momentum().length().dump();
          DT_LOG_TRACE(get_logging_priority(), "Charge"); iseq->charge().dump();
          if( iseq->has_helix_vertex() ){
            DT_LOG_TRACE(get_logging_priority(), "Helix vertex " << iseq->helix_vertex_type()); iseq->helix_vertex().dump();
            if( iseq->helix_vertex_type() == "calo" ) {
              DT_LOG_TRACE(get_logging_priority(), "icalo " << iseq->helix_vertex_id());
            }
          }
          if( iseq->has_decay_helix_vertex() ){
            DT_LOG_TRACE(get_logging_priority(), "Decay helix vertex " << iseq->decay_helix_vertex_type()); iseq->decay_helix_vertex().dump();
            if( iseq->decay_helix_vertex_type() == "calo" ) {
              DT_LOG_TRACE(get_logging_priority(), "icalo " << iseq->calo_helix_id());
            }
          }
          if( iseq->has_tangent_vertex() ){
            DT_LOG_TRACE(get_logging_priority(), "Tangent_vertex " << iseq->tangent_vertex_type());
            if( iseq->tangent_vertex_type() == "calo" ) {
              DT_LOG_TRACE(get_logging_priority(), "icalo " << iseq->tangent_vertex_id());
            }
          }
          if( iseq->has_decay_tangent_vertex() ){
            DT_LOG_TRACE(get_logging_priority(), "Decay tangent vertex " << iseq->decay_tangent_vertex_type());
            if( iseq->decay_tangent_vertex_type() == "calo" )
              DT_LOG_TRACE(get_logging_priority(), "icalo " << iseq->calo_tangent_id());
          }
          if( iseq->has_tangent_length() ){
            DT_LOG_TRACE(get_logging_priority(), "Tangent length "); iseq->tangent_length().dump();
          }
          if( iseq->has_helix_length() ){
            DT_LOG_TRACE(get_logging_priority(), "Helix length "); iseq->helix_length().dump();
          }
        }
        ++iseq;
        continue;
      }


    clock.stop(" sequentiator: interpret physics ");

    return;

  }

  //*************************************************************
  void sequentiator::add_pair(const topology::sequence & newsequence){
    //*************************************************************
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    clock.start(" sequentiator: add pair ", "cumulative");

    if (newsequence.nodes().size() != 2){
      DT_LOG_DEBUG(get_logging_priority(), "Problem: pair has size " << newsequence.nodes().size());
      clock.stop(" sequentiator: add pair ");
      return;
    }

    topology::cell_couplet cc;
    if( !newsequence.nodes_[0].has_couplet(newsequence.nodes()[1].c(), &cc) ){
      DT_LOG_DEBUG(get_logging_priority(), "Problem: node " << newsequence.nodes_[0].c().id() << " has no pair "
                   << newsequence.nodes()[0].c().id() << "-" << newsequence.nodes()[1].c().id());
      clock.stop(" sequentiator: add pair ");
      return;
    }

    topology::node na = newsequence.nodes()[0];
    topology::node nb = newsequence.nodes()[1];

    std::vector<topology::node> nodes;
    nodes.push_back(na);
    nodes.push_back(nb);

    bool erased = true;

    DT_LOG_DEBUG(get_logging_priority(), "n of tangents: " << cc.tangents_.size());
    topology::sequence pair(nodes);
    pair.set_probmin(probmin);
    for(std::vector<topology::line>::iterator itangent=cc.tangents_.begin(); itangent != cc.tangents_.end(); ++itangent){

      pair.nodes_[0].set_ep(itangent->epa());
      pair.nodes_[0].set_free(false);

      pair.nodes_[1].set_ep(itangent->epb());
      pair.nodes_[1].set_free(false);

      if( itangent - cc.tangents_.begin() > 0 ){
        NCOPY ++;
        for(std::vector<topology::cell_couplet>::iterator icc=pair.nodes_[0].cc_.begin(); icc!=pair.nodes_[0].cc_.end(); ++icc)
          icc->set_all_used();
        for(std::vector<topology::cell_triplet>::iterator iccc=pair.nodes_[1].ccc_.begin(); iccc!=pair.nodes_[1].ccc_.end(); ++iccc)
          iccc->set_all_used();
      }

      clock.start(" sequentiator: set free level ","cumulative");
      pair.set_free_level();
      clock.stop(" sequentiator: set free level ");

      make_name(pair);
      sequences_.push_back(pair);

      DT_LOG_DEBUG(get_logging_priority(), "made track [" << sequences_.size()-1 << "] with cells "
                   << na.c().id() << " [" << pair.nodes_[0].ep().x().value() << ", " << pair.nodes_[0].ep().z().value() << "]"
                   <<  "and " << nb.c().id() << " [" << pair.nodes_[1].ep().x().value() << ", " << pair.nodes_[1].ep().z().value() << "]");
      if( erased )
        erased = clean_up_sequences();


    }

    clock.stop(" sequentiator: add pair ");
    return;

  }


  //*************************************************************
  bool sequentiator::clean_up_sequences(){
    //*************************************************************

    clock.start(" sequentiator: clean up sequences ", "cumulative");

    if( sequences_.size() < 2 ){
      clock.stop(" sequentiator: clean up sequences ");
      return false;
    }

    bool changed = false;

    std::vector<topology::sequence>::iterator iseq = sequences_.begin();

    while( iseq != sequences_.end() ){

      // compare iseq with last sequence
      if( (size_t)(iseq - sequences_.begin() + 1) >= sequences_.size() ){
        // iseq is last sequence
        // no need to check it against last sequence
        break;
      }

      DT_LOG_DEBUG(get_logging_priority(), "Should we erase last sequence [" << sequences_.size() - 1 << "] " << sequences_.back().name()
                   << " in favour of [" << iseq - sequences_.begin() << "] " << iseq->name() << " ? ");
      if( !sequences_.back().contained( *iseq )){
        DT_LOG_DEBUG(get_logging_priority(), "No, it's not contained");
        //}else if( sequences_.back().Free() && !sequences_.back().contained_same_extreme_quadrants( *iseq ) ){
        //      m.message("no, it's a free sequence with different extreme quadrants " ,mybhep::VVERBOSE); fflush(stdout);
      }else if( sequences_.back().Free() ){
        DT_LOG_DEBUG(get_logging_priority(), "No, it's a free sequence");
      }else
        {
          DT_LOG_DEBUG(get_logging_priority(), "Erased last sequence " << sequences_.size() - 1 << " contained in sequence" << iseq - sequences_.begin());
          sequences_.pop_back();
          changed =  true;
          // now check again the same iseq against the new last sequence
          continue;
        }

      DT_LOG_DEBUG(get_logging_priority(), "Should we erase sequence [" << iseq - sequences_.begin() << "] "
                   << iseq->name() << " in favour of the last [" << sequences_.size() - 1 << "] " << sequences_.back().name());
      if( !iseq->contained( sequences_.back() ) ){
        DT_LOG_DEBUG(get_logging_priority(), "No, it's not contained");
        //}else if( iseq->Free() && !iseq->contained_same_extreme_quadrants( sequences_.back() ) ){
        //m.message("no, it's a free sequence with different extreme quadrants " ,mybhep::VVERBOSE); fflush(stdout);
      }else if( iseq->Free() ){
        DT_LOG_DEBUG(get_logging_priority(), "No, it's a free sequence");
      }else
        {
          DT_LOG_DEBUG(get_logging_priority(), "Erased sequence " << iseq - sequences_.begin() << "contained in last sequence" << sequences_.size()-1);
          sequences_.erase(iseq);
          changed = true;
          iseq = sequences_.begin() + (iseq - sequences_.begin());
          // now check the new iseq against the same last sequence
          continue;
        }


      DT_LOG_DEBUG(get_logging_priority(), "Should we erase sequence [" << iseq - sequences_.begin() << "] " << iseq->name() << " because of connection out of range ?");
      for(std::vector<topology::node>::iterator in=iseq->nodes_.begin(); in!=iseq->nodes_.end();in++){
        if( changed ) continue;
        if( in-iseq->nodes_.begin() +1 >= (int) iseq->nodes_.size()) break;
        topology::node nA = *in;
        topology::node nB = iseq->nodes_[in-iseq->nodes_.begin()+1];
        if( !sequence_is_within_range(nA,nB,*iseq) ){
          DT_LOG_DEBUG(get_logging_priority(), "Erased sequence " << iseq - sequences_.begin() << " not in range");
          sequences_.erase(iseq);
          changed = true;
          iseq = sequences_.begin() + (iseq - sequences_.begin());
          // now check the new iseq
          break;
        }
      }



      // arrive here if no sequences have been deleted
      // now check (kseq, iseq, lastseq) for bridges

      std::vector<topology::sequence>::iterator kseq = sequences_.begin();

      while( kseq != sequences_.end() ){

        if( iseq==kseq ){
          ++ kseq;
          continue;
        }

        if( (size_t)(kseq - sequences_.begin() + 1) >= sequences_.size() ){
          // kseq is last sequence
          break;
        }

        if( (size_t)(iseq - sequences_.begin() + 1) >= sequences_.size() ){
          // iseq is last sequence
          break;
        }

        DT_LOG_DEBUG(get_logging_priority(), "Should we erase sequence ["<< kseq - sequences_.begin() << "] " << kseq->name()
                     << " as bridge between ["<< iseq - sequences_.begin() << "] " << iseq->name()
                     << " and last ["<< sequences_.size() - 1<< "] "<< sequences_.back().name() << " ? ");
        if( kseq->is_bridge(*iseq, sequences_.back() ) &&
            (*kseq).nodes().size() < (*iseq).nodes().size() &&
            (*kseq).nodes().size() < sequences_.back().nodes().size()  &&
            !kseq->Free())
          {
            DT_LOG_DEBUG(get_logging_priority(), "Erased sequence  ["<< kseq - sequences_.begin()<< "] " << kseq->name()
                         << " as bridge between ["<< iseq - sequences_.begin() << "] " << iseq->name()
                         << " and last ["<< sequences_.size() - 1<< "] "<< sequences_.back().name());
            sequences_.erase(kseq);
            changed = true;
            kseq = sequences_.begin() + (kseq - sequences_.begin());
            if( kseq - sequences_.begin() < iseq - sequences_.begin() )
              iseq = sequences_.begin() + (iseq - sequences_.begin() - 1);
            else
              iseq = sequences_.begin() + (iseq - sequences_.begin());
            // now check the new kseq
            continue;
          }


        DT_LOG_DEBUG(get_logging_priority(), "Should we erase last ["<< sequences_.size() - 1<< "] "<< sequences_.back().name()
                     << " as bridge between ["<< iseq - sequences_.begin() << "] " << iseq->name()
                     << " and  ["<< kseq - sequences_.begin()<< "] " << kseq->name() << " ? ");
        if( sequences_.back().is_bridge(*iseq, *kseq ) &&
            sequences_.back().nodes().size() < (*iseq).nodes().size() &&
            sequences_.back().nodes().size() < (*kseq).nodes().size()  &&
            !sequences_.back().Free())
          {
            DT_LOG_DEBUG(get_logging_priority(), "Erased  last ["<< sequences_.size() - 1<< "] "<< sequences_.back().name()
                         << " as bridge between ["<< iseq - sequences_.begin() << "] " << iseq->name()
                         << " and  ["<< kseq - sequences_.begin()<< "] " << kseq->name());
            sequences_.pop_back();
            changed =  true;
            iseq = sequences_.begin() + (iseq - sequences_.begin());
            // now check the same kseq with the new last sequence
            continue;
          }

        kseq ++;
        continue;
      }

      iseq++;
      continue;
    }


    clock.stop(" sequentiator: clean up sequences ");
    return changed;

  }

  //*************************************************************
  double sequentiator::distance_from_foil(const topology::experimental_point &ep){
    //*************************************************************

    if( SuperNemo )
      return std::abs(ep.z().value());

    return std::abs(ep.radius().value() - FoilRadius);

  }

  //*************************************************************
  bool sequentiator::direct_out_of_foil(void){
    //*************************************************************

    clock.start(" sequentiator: direct out of foil ", "cumulative");

    for(std::vector<topology::sequence>::iterator iseq = sequences_.begin(); iseq != sequences_.end(); ++iseq){

      if( distance_from_foil(iseq->nodes().front().ep()) >
          distance_from_foil(iseq->nodes().back().ep()) ){
        DT_LOG_DEBUG(get_logging_priority(), "Sequence " << iseq - sequences_.begin() << " will be directed out of foil");
        topology::sequence is = iseq->invert();
        std::swap(*iseq, is);
      }
    }

    /*
      for(size_t iseq=0; iseq < sequences_.size(); iseq++){

      if( distance_from_foil(sequences_[iseq].nodes().front().ep()) <
      distance_from_foil(sequences_[iseq].nodes().back().ep()) ){
      m.message(" sequence ", iseq, " will be directed out of foil ", mybhep::VVERBOSE);
      print_a_sequence(&sequences_[iseq]);
      sequences_[iseq] = sequences_[iseq].invert();
      print_a_sequence(&sequences_[iseq]);

      }
      }

    */
    clock.stop(" sequentiator: direct out of foil ");


    return true;

  }


  //*************************************************************
  bool sequentiator::direct_scenarios_out_of_foil(void){
    //*************************************************************

    clock.start(" sequentiator: direct scenarios out of foil ", "cumulative");

    for(std::vector<topology::scenario>::iterator isc = scenarios_.begin(); isc != scenarios_.end(); ++isc){

      for(std::vector<topology::sequence>::iterator iseq = isc->sequences_.begin(); iseq != isc->sequences_.end(); ++iseq){

        if( distance_from_foil(iseq->nodes().front().ep()) >
            distance_from_foil(iseq->nodes().back().ep()) ){
          DT_LOG_DEBUG(get_logging_priority(), "Sequence " << iseq - sequences_.begin() << " in scenario " << isc - scenarios_.begin() << " will be directed out of foil");
          topology::sequence is = iseq->invert();
          std::swap(*iseq, is);
        }
      }

    }
    clock.stop(" sequentiator: direct scenarios out of foil ");


    return true;

  }


  //*************************************************************
  bool sequentiator::there_is_free_sequence_beginning_with(const topology::cell &c, size_t *index) {
    //*************************************************************

    clock.start(" sequentiator: there is free sequence beginning with ", "cumulative");

    for(std::vector<topology::sequence>::iterator iseq=sequences_.begin(); iseq != sequences_.end(); ++iseq)
      if( iseq->nodes()[0].c().id() == c.id() )
        {
          if( iseq->Free() )
            {
              *index = iseq - sequences_.begin();
              clock.stop(" sequentiator: there is free sequence beginning with ");
              return true;
            }
        }

    clock.stop(" sequentiator: there is free sequence beginning with ");
    return false;

  }

  //*************************************************************
  bool sequentiator::near(const topology::cell &c, topology::calorimeter_hit &ch){
    //*************************************************************

    topology::plane pl = ch.pl();
    double calorimeter_layer = ch.layer();

    if( pl.view() == "x" ){

      DT_LOG_DEBUG(get_logging_priority(), "Matching cell " << c.id() << " with cell number " << c.cell_number() << " to calo on view " << pl.view()
                   << " max cell number " << cell_max_number << " plane norm x " << pl.norm().x().value());

      if( pl.norm().x().value() > 0. )
        return (std::abs(cell_max_number + c.cell_number()) < 1 + NOffLayers);

      return (std::abs(cell_max_number - c.cell_number()) < 1 + NOffLayers);
    }
    else if( pl.view() == "y" ){

      topology::experimental_vector distance(c.ep(), pl.face());
      distance = distance.hor();
      double size_z = CellDistance + pl.sizes().z().value();
      double size_x = CellDistance + pl.sizes().x().value();

      DT_LOG_DEBUG(get_logging_priority(), "Checking if cell " << c.id() << " is near plane: "
                   << pl.center().x().value() << ", " << pl.center().y().value() << ", " << pl.center().z().value()
                   << " on view " << pl.view() << " distance z " << distance.z().value()
                   << " size z " << size_z << " distance x " << distance.x().value() << " size x " << size_x);

      if( std::abs(distance.z().value()) > size_z/2. ) return false;
      if( std::abs(distance.x().value()) > size_x/2. ) return false;

      return true;

    }
    else if( pl.view() == "z" ){

      int g = gap_number(c);

      DT_LOG_DEBUG(get_logging_priority(), "Checking if cell " << c.id() << " on gap " << g << " is near plane: "
                   << pl.center().x().value() << ", " << pl.center().y().value() << ", " << pl.center().z().value());

      if( g <= 0 || std::abs(calorimeter_layer - c.layer()) > NOffLayers ) return false; // cell is not on a gap or is facing the foil

      if( g == 1 || std::abs(calorimeter_layer - c.layer()) <= NOffLayers ) return true;

      DT_LOG_DEBUG(get_logging_priority(), "Problem: can't match to calo on view " << pl.view());

      return false;
    }
    else if( pl.view() == "inner" || pl.view() == "outer"){

      int ln = c.layer();
      int g = gap_number(c);
      DT_LOG_DEBUG(get_logging_priority(), "Checking if cell " << c.id() << " layer and gap: " << ln << " " << g << " is near plane: "
                   << pl.center().x().value() << ", " << pl.center().y().value() << ", " << pl.center().z().value() << " on view "<< pl.view());
      if( ln < 0 && (g == 3 ||
                     std::abs(ln - calorimeter_layer) <= NOffLayers )) return true;
      return false;
    }
    else if( pl.view() == "top" ||  pl.view() == "bottom" ){
      int ln = c.layer();
      int g = gap_number(c);
      DT_LOG_DEBUG(get_logging_priority(), "Checking if cell " << c.id() << " on gap " << g << " is near plane: "
                   << pl.center().x().value() << ", " << pl.center().y().value() << ", " << pl.center().z().value() << " on view " << pl.view());
      if( ln > 0 && calorimeter_layer == 3.5 && (g == 1 || std::abs(ln - calorimeter_layer) <= NOffLayers + 0.5 ) ) return true;
      if( ln > 0 && calorimeter_layer == 5.5 && (g == 2 || std::abs(ln - calorimeter_layer) <= NOffLayers + 0.5 ) ) return true;
      if( ln < 0 && calorimeter_layer == -3.5 && (g == 1 || std::abs(ln - calorimeter_layer) <= NOffLayers + 0.5 ) ) return true;
      if( ln < 0 && calorimeter_layer == -5.5 && (g == 2 || std::abs(ln - calorimeter_layer) <= NOffLayers + 0.5 ) ) return true;
      return false;
    }

    DT_LOG_DEBUG(get_logging_priority(), "Problem: can't match to calo on view " << pl.view());

    return false;

  }


  //*************************************************************
  int sequentiator::gap_number(const topology::cell &c){
    //*************************************************************

    // returns the index of the gap on which the hit is facing: 1, 2, 3
    // ... returns -1 if not on a gap
    // ... returns 0 if the hit is on layer 0, facing the foil

    size_t ln = abs(c.layer());

    DT_LOG_DEBUG(get_logging_priority(), "Cell " << c.id() << " layer " << c.layer());

    if( ln == 0 ) return 0;

    size_t counter = 0;
    for(size_t i=0; i<planes_per_block.size(); i++){
      counter += (int)planes_per_block[i];  // 4, 6, 9

      if( ln == counter - 1 )  // layer = 3, 5, 8
        return (i + 1);  // gap = 1, 2, 3
      if( ln == counter )  // layer = 4, 6, 9
        return (i + 1);  // gap = 1, 2, 3

    }

    return -1;

  }

  //*************************************************************
  bool sequentiator::sequence_is_within_range(topology::node nodeA, topology::node nodeB, topology::sequence seq){
    //*************************************************************

    if( gaps_Z.size() == 0 ) return true;

    topology::experimental_point epA = nodeA.ep();
    topology::experimental_point epB = nodeB.ep();
    topology::cell cA = nodeA.c();
    topology::cell cB = nodeB.c();

    int gnA = gap_number(cA);
    int gnB = gap_number(cB);
    int blockA = cA.block();
    int blockB = cB.block();

    double rmin, rmax;

    topology::experimental_point ep_maxr, ep_minr;

    // get the index of the gap through which the matching should occurr
    if( blockA == blockB &&
        ((gnA == -1 && gnB == -1) ||
         gnA != gnB )){
      // matching is within a single block and not through a gap
      // so helix should be contained in the layers of the two cells
      rmin = std::min(cA.ep().radius().value(),cB.ep().radius().value()) - CellDistance;
      rmax = std::max(cA.ep().radius().value(),cB.ep().radius().value()) + CellDistance;
    }else{
      if( blockA != blockB ){

        /*
        if( gnA != -1 && gnB != -1 && gnA != gnB ){
          if( level >= mybhep::VVERBOSE ){
            std::clog << " connection between cells " << nodeA.c().id() << " and " << nodeB.c().id() << " blocks " << blockA << " and " << blockB << " gaps " << gnA << " and " << gnB << " is forbidden (cells face different gaps) " << std::endl;
          }
          return false;
        }
        */

        // matching is between different blocks
        rmin = std::min(cA.ep().radius().value(),cB.ep().radius().value()) - CellDistance;
        rmax = std::max(cA.ep().radius().value(),cB.ep().radius().value()) + CellDistance;
      }
      else if( blockA == blockB && gnA == -1 && gnB != -1 ){ // B is on gap, A is inside
        size_t gn = gnB;
        if( cA.ep().radius().value() > cB.ep().radius().value() ){ // gap - B - A
          rmax = cA.ep().radius().value() + CellDistance;
          rmin = cB.ep().radius().value() - CellDistance - gaps_Z[gn];
        }else{ // A - B - gap
          rmin = cA.ep().radius().value() - CellDistance;
          rmax = cB.ep().radius().value() + CellDistance + gaps_Z[gn];
        }
      }
      else if( blockA == blockB && gnA != -1 && gnB == -1 ){ // A is on gap, B is inside
        size_t gn = gnA;
        if( cA.ep().radius().value() < cB.ep().radius().value() ){ // gap - A - B
          rmax = cB.ep().radius().value() + CellDistance;
          rmin = cA.ep().radius().value() - CellDistance - gaps_Z[gn];
        }else{ // B - A - gap
          rmin = cB.ep().radius().value() - CellDistance;
          rmax = cA.ep().radius().value() + CellDistance + gaps_Z[gn];
        }
      }
      else if( blockA == blockB && gnA != -1 && gnB != -1 && gnA == gnB ){ // A and B on same gap

        // get the index of the gap through which the matching should occurr

        size_t gn = gnA;
        size_t gaplayer=0;
        for(size_t i=0; i<gn; i++)
          gaplayer += (size_t)(planes_per_block[i]+0.5);  // 0, 4, 6, 9
        if( abs(cA.layer()) < (int)gaplayer ){ // foil - A,B - gap
          if( blockA > 0 ){ // origin - A, B - gap
            rmin = cA.ep().radius().value() - CellDistance;
            rmax = cA.ep().radius().value() + CellDistance + gaps_Z[gn];
          }else{ // origin - gap - A,B
            rmin = cA.ep().radius().value() - CellDistance - gaps_Z[gn];
            rmax = cA.ep().radius().value() + CellDistance;
          }
        }else{ // foil - gap - A, B
          if( blockA > 0 ){ // origin - gap - A,B
            rmax = cA.ep().radius().value() + CellDistance;
            rmin = cA.ep().radius().value() - CellDistance - gaps_Z[gn];
          }else{ // origin - A,B - gap
            rmin = cA.ep().radius().value() - CellDistance;
            rmax = cA.ep().radius().value() + CellDistance + gaps_Z[gn];
          }
        }

      }
      else{
        DT_LOG_DEBUG(get_logging_priority(), "Problem: blockA " << blockA << " blockB " << blockB << " gnA " << gnA << " gnB " << gnB);
        return true;
      }
    }

    seq.point_of_max_min_radius(epA, epB, &ep_maxr, &ep_minr);

    if( ep_maxr.radius().value() > rmax || ep_minr.radius().value() < rmin ){
      DT_LOG_DEBUG(get_logging_priority(), "Sequence penetrates outside of admissible range between cells " << nodeA.c().id()
                   << " layer " << nodeA.c().layer() << " r " << nodeA.c().ep().radius().value() << " and " << nodeB.c().id()
                   << " layer " << nodeB.c().layer() << " r " << nodeB.c().ep().radius().value()
                   << ": point of max radius has radius " << ep_maxr.radius().value() << " rmax " << rmax
                   << " point of min radius has radius " << ep_minr.radius().value() << " rmin " << rmin);
      return false;
    }

    DT_LOG_DEBUG(get_logging_priority(), "Sequence blockA " << blockA << " blockB " << blockB << " gnA " << gnA << " gnB " << gnB
                 << " connecting cells " << nodeA.c().id() << " and " << nodeB.c().id()
                 << ": point of max radius has radius " << ep_maxr.radius().value() << " rmax " << rmax
                 << " point of min radius has radius " << ep_minr.radius().value() << " rmin " << rmin);
    return true;

  }

  //*************************************************************
  bool sequentiator::good_first_to_be_matched(topology::sequence& seq){
    //*************************************************************

    if( !seq.fast() ) return false;

    if( seq.nodes().size() < 2 )
      {
        DT_LOG_DEBUG(get_logging_priority(), "SQ user: not a good first, because hits size is" << seq.nodes().size());
        return false;
      }

    return true;

  }

  //*************************************************************
  bool sequentiator::match_gaps(std::vector<topology::calorimeter_hit> & calos){
    //*************************************************************

    //if( gaps_Z.size() <= 1 ) return true;
    if( sequences_.size() < 2 ) return true;

    clock.start(" sequentiator: match gaps ", "cumulative");

    std::vector<bool> matched;
    matched.assign (families_.size(), false);

    std::vector<topology::sequence> newseqs;

    size_t ifam;
    size_t jmin;
    bool invertA, invertB, first;
    for(std::vector<topology::sequence>::iterator iseq=sequences_.begin(); iseq!=sequences_.end(); ++iseq){

      if( late() )
        return false;

      ifam = mybhep::int_from_string(iseq->family());

      if( matched[ifam] ) continue;

      if( !good_first_to_be_matched(*iseq) ){
        continue;
      }

      DT_LOG_DEBUG(get_logging_priority(), "Begin matching with sequence " << iseq->name());

      first = true;
      topology::sequence newseq = *iseq;
      int with_kink = 0;
      int cells_to_delete = 0;
      while( can_match(newseq, &jmin, invertA, invertB, with_kink, cells_to_delete, calos) )
        {
          DT_LOG_DEBUG(get_logging_priority(), "Best matching is " << sequences_[jmin].name()
                       << " invertA " << invertA << " invertB " << invertB << " with kink " << with_kink << " delete cells " << cells_to_delete);

          bool ok;
          newseq = newseq.match(sequences_[jmin], invertA, invertB, &ok, with_kink,cells_to_delete, Ratio);

          if( !ok && !with_kink ){
            DT_LOG_DEBUG(get_logging_priority(), "... no good helix match");
            continue;
          }

          if( first ){
            //      matched[ifam] = true;
            first = false;
          }
          size_t ifa = mybhep::int_from_string(sequences_[jmin].family());
          matched[ifa] = true;
          DT_LOG_DEBUG(get_logging_priority(), "Setting family " << ifa << " as used for matching");

        }

      if( !first )
        newseqs.push_back(newseq);

    }
    if( late() )
      return false;

    DT_LOG_DEBUG(get_logging_priority(), "Made matching through gaps");

    for(std::vector<topology::sequence>::iterator iseq=sequences_.begin(); iseq!=sequences_.end(); ++iseq){
      if( !matched[mybhep::int_from_string(iseq->family())] )
        newseqs.push_back(*iseq);
    }

    set_sequences(newseqs);
    make_families();

    //free(matched);

    return true;

  }

  //*************************************************************
  bool sequentiator::can_match(topology::sequence &s, size_t* jmin, bool& bestinvertA, bool& bestinvertB, int &with_kink, int &cells_to_delete_best, std::vector<topology::calorimeter_hit> & calos) {
    //*************************************************************

    if( late() )
      return false;

    clock.start(" sequentiator: can match ", "cumulative");

    bool ok = false;
    double limit_diagonal = sqrt(2.)*cos(M_PI/8.)*CellDistance;

    double probmax = -1.;
    double chi2min = mybhep::default_min;
    int ndofbest = 1;
    int cells_to_delete;

    DT_LOG_DEBUG(get_logging_priority(), "Try to match sequence" << s.name() << " of chi2 = " << chi2min << " ndof " << ndofbest << " prob " << probmax);
    bool invertA, invertB, acrossGAP;
    double p;
    double c;
    int n;
    topology::node nodeA, nodeB;

    for(std::vector<topology::sequence>::iterator jseq=sequences_.begin(); jseq!=sequences_.end(); ++jseq)
      {

        cells_to_delete = 0;
        with_kink=0;

        DT_LOG_DEBUG(get_logging_priority(), "Try to match sequence " << s.name() << " to " << jseq->name());

        // check that jseq's family is not the same as s
        if( s.same_families(*jseq) ){
          DT_LOG_DEBUG(get_logging_priority(), "... forbidden, same family " << jseq->family());
          continue;
        }

        // match along helix
        bool ok_match = s.good_match(*jseq, invertA, invertB, NOffLayers);
        bool ok_kink_match=false;
        bool ok_kink_match_chi2 = false;
        topology::sequence news;
        if( ok_match )
          news = s.match(*jseq, invertA, invertB, &ok_match,with_kink,0, Ratio);

        if( !invertA )
          nodeA = s.last_node();
        else
          nodeA = s.nodes_[0];

        if( !invertB )
          nodeB = jseq->nodes_[0];
        else
          nodeB = jseq->last_node();

        ok_match = ok_match && sequence_is_within_range(nodeA, nodeB, news);

        if( ok_match ){
          p = news.helix_Prob();
          c = news.helix_chi2();
          n = news.ndof();

          DT_LOG_DEBUG(get_logging_priority(), "... matched to " << jseq->name() << ", chi2 =" << c << " ndof " << n << " prob "<< p);
        }

        if( ok_match && (p > news.probmin())) {
          DT_LOG_DEBUG(get_logging_priority(), "Good helix match");
        } else {
          DT_LOG_DEBUG(get_logging_priority(), " ... no good helix match, try to match with kink");

          ok_kink_match= s.good_match_with_kink(*jseq, invertA, invertB, acrossGAP, limit_diagonal, NOffLayers, cells_to_delete);
          if( !ok_kink_match ){
            DT_LOG_DEBUG(get_logging_priority(), "... obviously no good match with kink");
          }
          else{
            // do not match with kink if the extreme is near a calo or near the foil

            if( !invertA ){
              if( cells_to_delete == 0 || cells_to_delete == 1 )
                nodeA = s.last_node();
              else
                nodeA = s.second_last_node();
            }
            else{
              if( cells_to_delete == 0 || cells_to_delete == 1 )
                nodeA = s.nodes_[0];
              else
                nodeA = s.nodes_[1];
            }

            if( !invertB ){
              if( cells_to_delete != 1 )
                nodeB = jseq->nodes_[0];
              else
                nodeB = jseq->nodes_[1];
            }else{
              if( cells_to_delete != 1 )
                nodeB = jseq->last_node();
              else
                nodeB = jseq->second_last_node();
            }

            DT_LOG_DEBUG(get_logging_priority(), "Possible match with kink, across GAP" << acrossGAP << ", cells to delete " << cells_to_delete << ", try to extrapolate");
            for(std::vector<topology::calorimeter_hit>::iterator ic=calos.begin(); ic != calos.end(); ++ic){
              if( near(nodeA.c(), *ic) ||  near(nodeB.c(), *ic) ){
                ok_kink_match = false;
                DT_LOG_DEBUG(get_logging_priority(), "  will not match with kink because end cell is near calo " << ic - calos.begin());
                break;
              }
            }
            if( gap_number(nodeA.c() ) == 0 || gap_number(nodeB.c() ) == 0 ){
              ok_kink_match = false;
              DT_LOG_DEBUG(get_logging_priority(), "  will not match with kink because end cell is near foil");
            }

            topology::experimental_point kink_point;
            ok_kink_match = ok_kink_match && s.intersect_sequence(*jseq, invertA, invertB, acrossGAP, &kink_point, limit_diagonal, &with_kink, cells_to_delete, Ratio);

            if( ok_kink_match ){

              news = s.match(*jseq, invertA, invertB, &ok_kink_match_chi2, with_kink,cells_to_delete, Ratio);
              ok_kink_match = sequence_is_within_range(nodeA, nodeB, news);

              if( ok_kink_match )
                DT_LOG_DEBUG(get_logging_priority(), "Good kink match (" << kink_point.x().value() << ", " << kink_point.y().value() << ", " << kink_point.z().value() << ") to sequence " << jseq->name());

            }else{
              DT_LOG_DEBUG(get_logging_priority(), "  no good kink match to sequence " << jseq->name());
            }
          }
        }

        if( !ok_match && !ok_kink_match ) continue;

        p = news.helix_Prob();
        c = news.helix_chi2();
        n = news.ndof();

        DT_LOG_DEBUG(get_logging_priority(), " ... matched to " << jseq->name()
                     << ", chi2 =" << c << " ndof " << n << " prob " << p
                     << " with kink " << with_kink << " cells_to_delete " << cells_to_delete
                     << " probmax " << probmax << " chi2min " << chi2min);

        if( (p > probmax || (p == probmax && c < chi2min ))  &&
            ((ok_match && p > news.probmin()) || ok_kink_match ) )
          {
            *jmin = jseq - sequences_.begin();
            probmax = p;
            chi2min = c;
            ndofbest = n;
            bestinvertA = invertA;
            bestinvertB = invertB;
            cells_to_delete_best = cells_to_delete;
            ok = true;
          }
      }

    if( ok ){
      DT_LOG_DEBUG(get_logging_priority(), "Sequence " << s.name() << " can be matched to " << sequences_[*jmin].name()
                   << ", chi2 =" << chi2min << " ndof " << ndofbest << " prob " << probmax << " cells_to_delete " << cells_to_delete);
    }

    clock.stop(" sequentiator: can match ");
    return ok;

  }

  size_t sequentiator::near_level( const topology::cell & c1, const topology::cell & c2 ){

    // returns 0 for far-away cell
    // 1 for diagonal cells
    // 2 for side-by-side cells

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
        DT_LOG_DEBUG(get_logging_priority(), "Problem: cat asking near level of cells with identical posiion ("
                     << hit1_side << ", " << hit1_layer << ", " << hit1_row << ") ("
                     << hit2_side << ", " << hit2_layer << ", " << hit2_row << ")");
        return 3;
      }
      else if (layer_distance == 1 && row_distance == 0) return 2;
      else if (layer_distance == 0 && row_distance == 1) return 2;
      else if (layer_distance == 1 && row_distance == 1) return 1;
      return 0;
  }

  void sequentiator::reassign_cells_based_on_helix( topology::sequence * seq ){

    if( seq->nodes_.size() < 3 ) return;

    topology::experimental_point helix_pos;
    topology::experimental_double distance;
    size_t index;
    std::vector<topology::node>::iterator inode = seq->nodes_.begin();
    inode ++;
    while( (size_t)(inode - seq->nodes_.begin()) < seq->nodes_.size()-1 && seq->nodes_.size() >= 3){

      helix_pos = seq->get_helix().position(inode->ep());
      distance = helix_pos.distance(inode->ep());
      index = inode - seq->nodes_.begin();

      DT_LOG_DEBUG(get_logging_priority(), "Sequence of " << seq->nodes_.size() << " nodes, index " << index
                   << " id " << inode->c().id() << " distance from helix " << distance.value() << " +- " << distance.error());

      if( distance.value() > CellDistance/2. &&
	  near_level(seq->nodes_[index-1].c(), seq->nodes_[index+1].c()) ){
        DT_LOG_DEBUG(get_logging_priority(), "Sequence of " << seq->nodes_.size() << " nodes, index " << index
                     << " id " << inode->c().id() << " distance from helix " << distance.value() << " +- " << distance.error() << " remove node");
	seq->nodes_.erase(inode);
	inode = seq->nodes_.begin() + (inode - seq->nodes_.begin());
	index = inode - seq->nodes_.begin();
	inode = seq->nodes_.begin() + index;
	if( index == 0 || index + 2 >= seq->nodes_.size() )
	  break;
	continue;
      }

      ++inode;
      continue;
    }

  }


} // end of namespace CAT
