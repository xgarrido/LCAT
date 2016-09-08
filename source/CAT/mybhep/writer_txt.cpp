/* -*- mode: c++ -*- */
/*
 *
 * Copyright 2006
 * J.J. Gomez-Cadenas
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <mybhep/writer_txt.h>

namespace mybhep{


  using namespace std;
  //! constructor
  writer_txt::writer_txt() :
    sequential_writer()
  {
  }
  //! destructor
  writer_txt::~writer_txt()
  {
    close();
  }

  void writer_txt::open_file(std::string fileName)
  {
    os_.open(fileName.c_str(), ios::out);
  }


  void writer_txt::close_file()
  {
    os_.close();
  }

    //! write the event as a record
  void writer_txt::write_record(std::string record)
  {
    os_ << record << endl;
  }


}
