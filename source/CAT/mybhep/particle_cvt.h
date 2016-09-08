/* -*- mode: c++ -*- */
//
/*
 *
 * Copyright (C) 2004 J.J. Gomez Cadenas
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
 */


#ifndef __PARTICLE_CVT_H__
#define __PARTICLE_CVT_H__

#include <string>
#include <map>
#include <vector>
#include <mybhep/particle.h>
#include <mybhep/ibconverter.h>


namespace mybhep{

  //! class particle_cvt
  /*!
    Converter to persistency
  */

  typedef particle particle;

  class particle_cvt : public ibconverter
  {

  private:
    //! pointer to the particle
    particle*   nparticle_;

    //! const pointer to the particle
    const particle*   particle_;

    //! integer for mother particle
    int imother_;

    //! integer for mirror particle
    std::vector<int> imirror_;

    //! integer for daguther particle
    std::vector<int> idaughters_;

    //! std::vector of track integers
    std::vector<int> tracks_;

    //! std::vector of detectors
    std::vector<std::string> idets_;

    //! map of hit integers
    std::multimap<std::string,int> hits_;

    //! label for mirrors
    std::string mirror_label_;

    //! label for tracks
    std::string track_label_;

    //! label for hits
    std::string hit_label_;

    //! label for hits size
    std::string hit_size_label_;

    size_t  number_of_hits( std::string& def );
  public:
    //! construct from a particle and an index
    particle_cvt( const particle& p, int inx);
    //! construct from a string
    particle_cvt( std::string );
    //! destructor
    virtual ~particle_cvt();

    //! store class in a string
    virtual std::string store();
    //! restore from a string
    virtual void restore( std::string );
    //! complete restore
    virtual void complete_restore();

    //! retrieve a constant reference to the particle
    const particle& retrieve() const;

    //! retrieve a constant reference to the particle
    particle& reference() ;

    //! instantiate a new particle
    particle* create();


  };
}
#endif









