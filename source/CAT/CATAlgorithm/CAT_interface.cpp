/* -*- mode: c++ -*- */
// CAT_interface.cpp

// Standard library:
#include <stdexcept>
#include <limits>
#include <map>
#include <vector>
#include <sstream>

// Third party:
// - Boost:
#include <boost/algorithm/string.hpp>
// - Bayeux/datatools:
#include <bayeux/datatools/exception.h>

// This project:
#include <CATAlgorithm/CAT_interface.h>

namespace CAT {

  const std::string & setup_data::get_error_message () const
  {
    return _error_message;
  }

  void setup_data::_set_error_message (const std::string & message_)
  {
    std::ostringstream oss;
    oss << "CAT::setup_data: ";
    oss << message_;
    _error_message = oss.str ();
    return;
  }

  setup_data::setup_data ()
  {
    _set_defaults ();
    return;
  }

  void setup_data::reset ()
  {
    _set_defaults ();
    return;
  }

  void setup_data::_set_defaults ()
  {
    _error_message.clear ();
    level                 = datatools::logger::PRIO_WARNING;
    SuperNemo             = true;
    MaxTime               = 5000.0 * CLHEP::ms;
    SmallRadius           =    2.0 * CLHEP::mm;
    TangentPhi            =   20.0 * CLHEP::degree;
    TangentTheta          =  160.0 * CLHEP::degree;
    SmallNumber           =    0.1 * CLHEP::mm;
    QuadrantAngle         =   90.0 * CLHEP::degree;
    Ratio                 =   10000.0;
    CompatibilityDistance = 4.0 * CLHEP::mm;
    MaxChi2               = 3.;
    probmin                = 0.;
    nofflayers            = 1;
    len                   = 2503. * CLHEP::mm;
    rad                   = 30.   * CLHEP::mm;
    vel                   = 0.06  * CLHEP::mm;
    CellDistance          = 30.   * CLHEP::mm;
    FoilRadius            = 0.;

    // SuperNEMO geometry default parameters :
    num_blocks = 1;
    planes_per_block.clear ();
    planes_per_block.push_back (num_blocks);
    planes_per_block.at(0) = 9;
    num_cells_per_plane    = 113;
    cell_size              = 44.0  * CLHEP::mm;
    bfield                 = 0.0025; // Tesla
    xsize                  = 2500. * CLHEP::mm; // this is y in SnWare coordinates
    ysize                  = 1350. * CLHEP::mm; // this is z in SnWare coordinates
    zsize                  = 450. * CLHEP::mm; // this is x in SnWare coordinates

    return;
  }

  bool setup_data::check () const
  {
    setup_data * mutable_this = const_cast<setup_data *> (this);
    return mutable_this->_check_snemo ();
  }

  bool setup_data::_check_snemo ()
  {
    if (SmallRadius <= 0.0)
      {
        _set_error_message ("Invalid 'SmallRadius'");
        return false;
      }
    if (TangentPhi <= 0.0)
      {
        _set_error_message ("Invalid 'TangentPhi'");
        return false;
      }
    if (TangentTheta <= 0.0)
      {
        _set_error_message ("Invalid 'TangentTheta'");
        return false;
      }
    if (SmallNumber <= 0.0)
      {
        _set_error_message ("Invalid 'SmallNumber'");
        return false;
      }
    if (QuadrantAngle <= 0.0)
      {
        _set_error_message ("Invalid 'QuadrantAngle'");
        return false;
      }
    if (Ratio <= 0.0)
      {
        _set_error_message ("Invalid 'Ratio'");
        return false;
      }
    if (CompatibilityDistance <= 0.0)
      {
        _set_error_message ("Invalid 'CompatibilityDistance'");
        return false;
      }
    if (MaxChi2 <= 0.0)
      {
        _set_error_message ("Invalid 'MaxChi2'");
        return false;
      }
    if (probmin < 0.0)
      {
        _set_error_message ("Invalid 'probmin'");
        return false;
      }
    if (nofflayers < 0.0)
      {
        _set_error_message ("Invalid 'nofflayers'");
        return false;
      }
    if (num_blocks < 1)
      {
        _set_error_message ("Invalid 'num_blocks'");
        return false;
      }
    if ((int) planes_per_block.size () != num_blocks)
      {
        _set_error_message ("Invalid size of 'planes_per_block'");
        return false;
      }
    if (num_cells_per_plane < 1)
      {
        _set_error_message ("Invalid 'num_cells_per_plane'");
        return false;
      }
    if (cell_size <= 0.0)
      {
        _set_error_message ("Invalid 'cell_size'");
        return false;
      }
    if (len < 0.0)
      {
        _set_error_message ("Invalid 'len'");
        return false;
      }
    if (rad < 0.0)
      {
        _set_error_message ("Invalid 'rad'");
        return false;
      }
    if (vel < 0.0)
      {
        _set_error_message ("Invalid 'vel'");
        return false;
      }
    if (CellDistance < 0.0)
      {
        _set_error_message ("Invalid 'CellDistance'");
        return false;
      }
    if (FoilRadius < 0.0)
      {
        _set_error_message ("Invalid 'FoilRadius'");
        return false;
      }

    return true;
  }

  void clusterizer_configure (clusterizer & czer_,
                              const setup_data & setup_)
  {
    DT_THROW_IF(! setup_.check(), std::logic_error,
                "Invalid setup data :" << setup_.get_error_message());

    // General parameters :
    czer_.set_logging_priority(setup_.level);

    // Algorithm parameters :
    czer_.set_tangent_phi(setup_.TangentPhi / CLHEP::degree);
    czer_.set_tangent_theta(setup_.TangentTheta / CLHEP::degree);
    czer_.set_quadrant_angle(setup_.QuadrantAngle / CLHEP::degree);
    czer_.set_ratio(setup_.Ratio);
    return;
  }

  void sequentiator_configure (sequentiator & stor_,
                               const setup_data & setup_)
  {
    if (! setup_.check ())
      {
        std::ostringstream emess;
        emess << "ERROR: CAT::sequentiator_configure: Invalid setup data :"
             << setup_.get_error_message ();
        throw std::logic_error(emess.str());
      }

    // General parameters :
    stor_.set_PrintMode (false);
    stor_.set_MaxTime (setup_.MaxTime / CLHEP::ms);
    // std::string leveltmp = setup_.level;
    // boost::to_upper(leveltmp);

    // stor_.set_level (leveltmp); //mybhep::get_info_level (leveltmp));

    // Algorithm parameters :
    stor_.set_SmallRadius (setup_.SmallRadius / CLHEP::mm);
    stor_.set_TangentPhi (setup_.TangentPhi / CLHEP::degree);
    stor_.set_TangentTheta (setup_.TangentTheta / CLHEP::degree);
    stor_.set_SmallNumber (setup_.SmallNumber / CLHEP::mm);
    stor_.set_QuadrantAngle (setup_.QuadrantAngle / CLHEP::degree);
    stor_.set_Ratio (setup_.Ratio);
    stor_.set_CompatibilityDistance (setup_.CompatibilityDistance);
    stor_.set_MaxChi2 (setup_.MaxChi2);
    stor_.set_probmin (setup_.probmin);
    stor_.set_nofflayers (setup_.nofflayers);
    stor_.set_len (setup_.len);
    stor_.set_rad (setup_.rad);
    stor_.set_vel (setup_.vel);
    stor_.set_CellDistance (setup_.CellDistance);
    stor_.set_FoilRadius (setup_.FoilRadius);

    stor_.set_bfield (setup_.bfield);
    stor_.set_xsize (setup_.xsize);
    stor_.set_ysize (setup_.ysize);
    stor_.set_zsize (setup_.zsize);

    // Geometry description :
    if (setup_.SuperNemo)
      {

        /// Activate the special new mode :
        stor_.set_SuperNemoChannel(true);

        // Layout of the tracking chamber :
        stor_.set_num_blocks (setup_.num_blocks);
        for (int i = 0; i < setup_.num_blocks; i++)
          {
            stor_.set_planes_per_block (i, (int)(setup_.planes_per_block.at (i)+0.5));
          }
        stor_.set_num_cells_per_plane (setup_.num_cells_per_plane);
        stor_.set_GG_CELL_pitch (setup_.cell_size / CLHEP::mm);
      }
    else
      {
        throw std::logic_error ("CAT::sequentiator_configure: Only SuperNEMO setup is supported !");
      }

    return;
  }


  /***********************************************************/

  topology::cell & input_data::add_cell ()
  {
    if (cells.size () == 0)
      {
        // memory preallocation at the first cell
        cells.reserve (50);
      }
    {
      topology::cell tmp;
      cells.push_back (tmp);
    }
    return cells.back ();
  }

  topology::calorimeter_hit & input_data::add_calo_cell ()
  {

    if (calo_cells.size () == 0)
      {
        // memory preallocation at the first calo_cell
        calo_cells.reserve (50);
      }
    {
      topology::calorimeter_hit tmp;
      calo_cells.push_back (tmp);
    }
    return calo_cells.back ();
  }

  bool input_data::check () const{
    return gg_check() && calo_check();
  }

  bool input_data::gg_check () const
  {
    // A map would be better to check cell IDs :
    std::map<int,bool> mids;
    for (int i = 0; i < (int) cells.size (); i++)
      {
        const topology::cell & c = cells.at(i);
        int cell_id = c.id();
        if (cell_id < 0 || cell_id > 10000)
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "Out of range cell ID '" <<  cell_id << "' !"
                      << std::endl;
            return false;
          }
        if (mids.find (cell_id) != mids.end ())
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "mids Duplicate cell ID '" <<  cell_id << "' !"
                      << std::endl;
            return false;
          }
        mids[cell_id] = true;
      }

    // Duplicate test for now :
    std::vector<bool> ids;
    ids.assign (cells.size (), false);
    for (int i = 0; i < (int) cells.size (); i++)
      {
        const topology::cell & c = cells.at(i);
        int cell_id = c.id();
        if ((cell_id < 0) || (cell_id >= (int) cells.size ()))
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "Invalid cell ID '" <<  cell_id << "' !"
                      << std::endl;
            return false;
          }
        if (ids[cell_id])
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "ids Duplicate cell ID '" <<  cell_id << "' !"
                      << std::endl;
            return false;
          }
        ids[cell_id] = true;
      }
    for (int i = 0; i < (int) ids.size (); i++)
      {
        if (! ids[i])
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "Cell ID '" << i << "' is not used ! There are some missing cells !"
                      << std::endl;
            return false;
          }
      }
    return true;
  }

  bool input_data::calo_check () const
  {
    // A map would be better to check cell IDs :
    std::map<int,bool> mids;
    for (int i = 0; i < (int) calo_cells.size (); i++)
      {
        const topology::calorimeter_hit & c = calo_cells.at(i);
        int calo_cell_id = c.id();
        if (calo_cell_id < 0 || calo_cell_id > 10000)
          {
            std::cerr << "ERROR: CAT::input_data::calo_check: "
                      << "Out of range calo_cell ID '" <<  calo_cell_id << "' !"
                      << std::endl;
            return false;
          }
        if (mids.find (calo_cell_id) != mids.end ())
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "Duplicate calo_cell ID '" <<  calo_cell_id << "' !"
                      << std::endl;
            return false;
          }
        mids[calo_cell_id] = true;
      }

    // Duplicate test for now :
    std::vector<bool> ids;
    ids.assign (calo_cells.size (), false);
    for (int i = 0; i < (int) calo_cells.size (); i++)
      {
        const topology::calorimeter_hit & c = calo_cells.at(i);
        int calo_cell_id = c.id();
        if ((calo_cell_id < 0) || (calo_cell_id >= (int) calo_cells.size ()))
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "Invalid calo_cell ID '" <<  calo_cell_id << "' !"
                      << std::endl;
            return false;
          }
        if (ids[calo_cell_id])
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "Duplicate calo_cell ID '" <<  calo_cell_id << "' !"
                      << std::endl;
            return false;
          }
        ids[calo_cell_id] = true;
      }
    for (int i = 0; i < (int) ids.size (); i++)
      {
        if (! ids[i])
          {
            std::cerr << "ERROR: CAT::input_data::check: "
                      << "Calo_Cell ID '" << i << "' is not used ! There are some missing calo_cells !"
                      << std::endl;
            return false;
          }
      }
    return true;
  }

  input_data::input_data ()
  {
    return;
  }

  /***********************************************************/

  output_data::output_data ()
  {
    return;
  }

}

// end of CAT_interface.cpp
