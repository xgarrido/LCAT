#include <CAT/cell.h>

int main()
{
  // First geiger cell with default contructor
  {
    CAT::cell a_cell;
    a_cell.tree_dump(std::clog, "Default cell:");
  }
  // Another geiger cell with fake data
  {
    CAT::cell a_cell;
    a_cell.set_id(666);
    a_cell.set_side(1);
    a_cell.set_layer(2);
    a_cell.set_row(3);
    a_cell.set_prompt(true);
    a_cell.set_position(geomtools::vector_3d(1.0*CLHEP::cm, 2.0*CLHEP::cm, 3.0*CLHEP::cm));
    a_cell.set_radius(1.5*CLHEP::cm);
    a_cell.set_radius_error(0.3*CLHEP::cm);
    a_cell.tree_dump(std::clog, "A fake cell:");
  }
}
