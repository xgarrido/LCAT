// Ourselves
#include <CAT/tracked_data.h>

namespace CAT {
  namespace topology {

    cell & tracked_data::add_gg_cell()
    {
      if (cells_.empty()) {
        // memory preallocation at the first cell
        cells_.reserve(50);
      }
      cell tmp;
      cells_.push_back(tmp);
      return cells_.back();
    }

    calorimeter_hit & tracked_data::add_calo_cell()
    {
      if (calos_.empty()) {
        // memory preallocation at the first cell
        calos_.reserve(10);
      }
      calorimeter_hit tmp;
      calos_.push_back(tmp);
      return calos_.back();
    }

  }
}
