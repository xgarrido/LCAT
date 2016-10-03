// Ourselves
#include <CAT/tracked_data.h>

// Standard library
#include <sstream>

// Third party
// - Bayeux/datatools:
#include <bayeux/datatools/i_tree_dump.h>

namespace CAT {

  tracked_data::tracked_data()
  {
    return;
  }

  tracked_data::~tracked_data()
  {
    reset();
    return;
  }

  void tracked_data::reset()
  {
    _gg_hits_.clear();
    _calo_hits_.clear();
    _clusters_.clear();
    _scenarios_.clear();
    return;
  }

  void tracked_data::tree_dump(std::ostream & out_,
                               const std::string & title_,
                               const std::string & indent_,
                               bool /*a_inherit*/) const
  {
    std::string indent;
    if (! indent_.empty()) indent = indent_;
    if (! title_.empty()) {
      out_ << indent << title_ << std::endl;
    }

    out_ << indent << datatools::i_tree_dumpable::tag << "Geiger hits : " << _gg_hits_.size() << std::endl;
    for (size_t i = 0; i < _gg_hits_.size(); i++) {
      out_ << indent << datatools::i_tree_dumpable::skip_tag;
      std::ostringstream indent_oss;
      indent_oss << indent;
      indent_oss << datatools::i_tree_dumpable::skip_tag;

      if (i == _gg_hits_.size() - 1) {
        out_ << datatools::i_tree_dumpable::last_tag;
        indent_oss << datatools::i_tree_dumpable::last_skip_tag;
      } else {
        out_ << datatools::i_tree_dumpable::tag;
        indent_oss << datatools::i_tree_dumpable::skip_tag;
      }
      out_ << "Hit #" << i << std::endl;
      _gg_hits_[i].tree_dump(out_, "", indent_oss.str());
    }
    // out_ << indent << " number of calos : " << calos_.size() << std::endl;
    // for(std::vector<calorimeter_hit>::const_iterator icalo=calos_.begin(); icalo!=calos_.end();++icalo)
    //   icalo->dump(out_, "",indent + "     ");
    // out_ << indent << " number of clusters : " << clusters_.size() << std::endl;
    // for(std::vector<cluster>::const_iterator icluster=clusters_.begin(); icluster!=clusters_.end();++icluster)
    //   icluster->dump(out_, "",indent + "     ");
    // out_ << indent << " number of scenarios : " << scenarios_.size() << std::endl;
    // for(std::vector<scenario>::const_iterator iscenario=scenarios_.begin(); iscenario!=scenarios_.end();++iscenario)
    //   iscenario->dump(out_, "",indent + "     ");
    // out_ << indent << " ------------------- " << std::endl;

    return;
  }

  cell & tracked_data::add_gg_hit()
  {
    if (_gg_hits_.empty()) {
      // memory preallocation at the first cell
      _gg_hits_.reserve(50);
    }
    cell tmp;
    _gg_hits_.push_back(tmp);
    return _gg_hits_.back();
  }

  std::vector<cell> & tracked_data::grab_gg_hits()
  {
    return _gg_hits_;
  }

  const std::vector<cell> & tracked_data::get_gg_hits() const
  {
    return _gg_hits_;
  }

  calorimeter_hit & tracked_data::add_calo_hit()
  {
    if (_calo_hits_.empty()) {
      // memory preallocation at the first cell
      _calo_hits_.reserve(10);
    }
    calorimeter_hit tmp;
    _calo_hits_.push_back(tmp);
    return _calo_hits_.back();
  }

  std::vector<calorimeter_hit> & tracked_data::grab_calo_hits()
  {
    return _calo_hits_;
  }

  const std::vector<calorimeter_hit> & tracked_data::get_calo_hits() const
  {
    return _calo_hits_;
  }

  std::vector<cluster> & tracked_data::grab_clusters()
  {
    return _clusters_;
  }

  const std::vector<cluster> & tracked_data::get_clusters() const
  {
    return _clusters_;
  }

  std::vector<scenario> & tracked_data::grab_scenarios()
  {
    return _scenarios_;
  }

  const std::vector<scenario> & tracked_data::get_scenarios() const
  {
    return _scenarios_;
  }

}
