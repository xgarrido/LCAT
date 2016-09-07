// Ourselves
#include <CAT/clusterizer.h>

// #include <mybhep/system_of_units.h>
// #include <sys/time.h>
// #include <limits>
// #include <cmath>
// #include <map>

// #if CAT_WITH_DEVEL_ROOT == 1
// #include <TApplication.h>
// #include <TROOT.h>
// #include <TChain.h>
// #include <TH2.h>
// #include <TH1.h>
// #include <TGraph.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TFile.h>
// #include <TMath.h>
// #include <TBox.h>
// #include <TMarker.h>
// #endif // CAT_WITH_DEVEL_ROOT == 1

namespace CAT {

  // //! get cells
  // const std::vector<topology::cell>& clusterizer::get_cells()const
  // {
  //   return cells_;
  // }

  // //! set cells
  // void clusterizer::set_cells(const std::vector<topology::cell> & cells)
  // {
  //   //cells_.clear();
  //   cells_ = cells;
  // }

  // //! get clusters
  // const std::vector<topology::cluster>& clusterizer::get_clusters()const
  // {
  //   return clusters_;
  // }

  // //! set clusters
  // void clusterizer::set_clusters(const std::vector<topology::cluster> & clusters)
  // {
  //   //clusters_.clear();
  //   clusters_ = clusters;
  // }

  // //! get calorimeter_hits
  // const std::vector<topology::calorimeter_hit>& clusterizer::get_calorimeter_hits()const
  // {
  //   return calorimeter_hits_;
  // }

  // //! set calorimeter_hits
  // void clusterizer::set_calorimeter_hits(const std::vector<topology::calorimeter_hit> & calorimeter_hits)
  // {
  //   calorimeter_hits_.clear();
  //   calorimeter_hits_ = calorimeter_hits;
  // }

  // void clusterizer::_set_defaults ()
  // {

  //   level = mybhep::NORMAL;
  //   m = mybhep::messenger(level);
  //   num_blocks = -1;
  //   planes_per_block.clear ();
  //   gaps_Z.clear ();
  //   GG_CELL_pitch = std::numeric_limits<double>::quiet_NaN ();
  //   GG_GRND_diam = std::numeric_limits<double>::quiet_NaN ();
  //   GG_CELL_diam = std::numeric_limits<double>::quiet_NaN ();
  //   CHAMBER_X = std::numeric_limits<double>::quiet_NaN ();
  //   GG_BLOCK_X = std::numeric_limits<double>::quiet_NaN ();
  //   num_cells_per_plane = -1;
  //   SOURCE_thick = std::numeric_limits<double>::quiet_NaN ();
  //   lastlayer = 0;
  //   vel  = std::numeric_limits<double>::quiet_NaN ();
  //   rad  = std::numeric_limits<double>::quiet_NaN ();
  //   len  = std::numeric_limits<double>::quiet_NaN ();
  //   CellDistance  = std::numeric_limits<double>::quiet_NaN ();
  //   xsize = ysize = zsize = std::numeric_limits<double>::quiet_NaN ();
  //   calo_X = calo_Y = calo_Z = std::numeric_limits<double>::quiet_NaN ();
  //   InnerRadius = OuterRadius= FoilRadius = std::numeric_limits<double>::quiet_NaN ();
  //   pmax = std::numeric_limits<double>::quiet_NaN ();

  //   SmallRadius = std::numeric_limits<double>::quiet_NaN ();
  //   TangentPhi = std::numeric_limits<double>::quiet_NaN ();
  //   TangentTheta = std::numeric_limits<double>::quiet_NaN ();
  //   SmallNumber = std::numeric_limits<double>::quiet_NaN ();
  //   QuadrantAngle = std::numeric_limits<double>::quiet_NaN ();
  //   Ratio = std::numeric_limits<double>::quiet_NaN ();
  //   CompatibilityDistance = std::numeric_limits<double>::quiet_NaN ();
  //   MaxChi2 = std::numeric_limits<double>::quiet_NaN ();
  //   probmin = std::numeric_limits<double>::quiet_NaN ();
  //   nofflayers = 0;
  //   first_event_number = 0;
  //   PrintMode = false;
  //   SuperNemo = true;
  //   SuperNemoChannel = false;
  //   NemoraOutput = false;
  //   N3_MC = false;
  //   MaxTime = std::numeric_limits<double>::quiet_NaN ();
  //   doDriftWires = true;
  //   DriftWires.clear ();
  //   eman = 0;
  //   _moduleNR.clear ();
  //   _MaxBlockSize = -1;
  //   event_number=0;
  //   hfile.clear ();

  //   nevent = 0;
  //   InitialEvents = 0;
  //   SkippedEvents = 0;
  //   run_list.clear ();
  //   run_time = std::numeric_limits<double>::quiet_NaN ();
  //   first_event=true;

  //   return;
  // }



  clusterizer::clusterizer()
  {
    //_set_defaults();
    return;
  }

  clusterizer::~clusterizer()
  {
    return;
  }

  void clusterizer::initialize()
  {

    // m.message("CAT::clusterizer::initialize: Entering...",mybhep::NORMAL);

    // m.message("CAT::clusterizer::initialize: Beginning algorithm clusterizer \n",mybhep::VERBOSE);

    // //----------- read dst param -------------//

    // readDstProper();

    // //------- end of read pram -----------//

    // _initialize();

    // m.message("CAT::clusterizer::initialize: Done.",mybhep::NORMAL);

    return;
  }


  // //*************************************************************
  // bool clusterizer::finalize() {
  //   //*************************************************************

  //   clock.start(" clusterizer: finalize ");

  //   m.message("CAT::clusterizer::finalize: Ending algorithm clusterizer...",mybhep::NORMAL);

  //   m.message("CAT::clusterizer::finalize: Initial events: ", InitialEvents, mybhep::NORMAL);
  //   m.message("CAT::clusterizer::finalize: Skipped events: ", SkippedEvents, "(", 100.*SkippedEvents/InitialEvents, "%)", mybhep::NORMAL);

  //   if( PrintMode )
  //     {
  //       finalizeHistos();
  //     }
  //   clock.stop(" clusterizer: finalize ");

  //   if( level >= mybhep::NORMAL ){
  //     clock.dump();
  //   }

  //   _set_defaults ();
  //   return true;
  // }




//   //*************************************************************
//   double clusterizer::long_resolution(double Z, double d[3])const{
//     //*************************************************************

//     double xp = abs(Z*2./len);
//     double kx = k0 + k1*xp + k2*xp*xp + k3*xp*xp*xp;
//     double thx = th0 + th1*xp + th2*xp*xp + th3*xp*xp*xp;
//     double tgth2=(d[0]*d[0] + d[2]*d[2])/(d[1]*d[1]);
//     double sigma=sqrt(kx*sigma0*sigma0*(1-xp*xp)+thx*thx/tgth2);

//     return sigma;
//   }


//   //*************************************************************
//   double clusterizer::long_resolution_1cthd(double Zdist)const{
//     //*************************************************************
//     double L = Zdist;
//     double sigma = l0 + l1*L;

//     return sigma;
//   }


//   //*************************************************************
//   double clusterizer::GetYError( double y, float tup, float tdown, double direction[3]){
//     //*************************************************************

//     double erry;

//     if( tup != 0. && tdown != 0.){
//       erry = long_resolution(y, direction);
//     }
//     else if( tup != 0. && tdown == 0.){
//       erry = long_resolution_1cthd(len/2.-y);
//     }
//     else if( tup == 0. && tdown != 0.){
//       erry = long_resolution_1cthd(len/2.+y);
//     }
//     else{
//       erry = len/2.;
//     }

//     return erry;

//   }


//   //*************************************************************
//   void clusterizer::FillYPositions( mybhep::event& evt ){
//     //*************************************************************

//     for(size_t i=0; i<evt.digi_particles().size(); i++){
//       if( evt.digi_particles()[i]->find_property("SUNAMI") )
//         FillYPositions( evt.digi_particles()[i] );
//     }

//     return;

//   }

//   //*************************************************************
//   void clusterizer::FillYPositions( mybhep::particle* p ){
//     //*************************************************************

//     for(size_t i=0; i<p->hits("trk").size(); i++){
//       FillYPosition(p->hits("trk")[i]);
//     }

//     return;

//   }

//   //*************************************************************
//   void clusterizer::FillYPosition( mybhep::hit* h ){
//     //*************************************************************

//     double erry = 0.;

//     if( SuperNemo ){
//       float tup, tdown, tdelay;

//       std::string ttime = h->fetch_property("TTIME");
//       std::string btime = h->fetch_property("BTIME");
//       std::string atime = h->fetch_property("ATIME");

//       tup = mybhep::float_from_string(ttime);
//       tdown = mybhep::float_from_string(btime);
//       tdelay = mybhep::float_from_string(atime);

//       if( std::isnan(tup) || std::isinf(tup) )
//         tup = -1.;
//       if( std::isnan(tdown) || std::isinf(tdown) )
//         tdown = -1.;
//       if( std::isnan(tdelay) || std::isinf(tdelay) )
//         tdelay = -1.;


//       std::vector<float> cellpos;
//       mybhep::vector_from_string(h->fetch_property("CELL_POS"), cellpos);

//       if( cellpos.size() != 3 ){
//         m.message("CAT::clusterizer::FillYPosition: problem: cell_pos size is ", cellpos.size(), mybhep::MUTE);
//         return;
//       }

//       ///////////////////////////////////////////////////////////
//       //                                                       //
//       //  len/2 - y = tup*v                                    //
//       //  y - (-len/2) = tdown*v                               //
//       //                                                       //
//       //  so:                                                  //
//       //                                                       //
//       //  (tup + tdown)*v = len                                //
//       //  (tup - tdown)*v = -2y                                //
//       //                                                       //
//       ///////////////////////////////////////////////////////////


//       if( tup != 0. && tdown != 0.)
//         cellpos[1] = (-1.)*((tup - tdown)/(tup + tdown))*(len/2.);
//       else if( tup != 0. && tdown == 0.)
//         cellpos[1] = (len/2.) - tup*vel;
//       else if( tup == 0. && tdown != 0.)
//         cellpos[1] = tdown*vel - (len/2.);
//       else
//         cellpos[1] = 0.;

//       h->change_property("CELL_POS", mybhep::vector_to_string(cellpos));

//       double dir[3];
//       dir[0] = 0.;
//       dir[1] = 0.;
//       dir[2] = 1.;

//       erry = GetYError((double)cellpos[1], tup, tdown, dir);

//       if( std::isnan( erry ) || std::isinf(erry) )
//         erry = -1.;

//       if( !h->find_property("ERRY") )
//         h->add_property("ERRY",mybhep::to_string(erry));
//     }
//     else{
//       erry =  mybhep::double_from_string(h->fetch_property("ERRY"));
//     }

//     return;

//   }

//   //*************************************************************
//   void clusterizer::FillTrueVertexes( mybhep::event& evt ){
//     //*************************************************************

//     for(size_t i=0; i<evt.true_particles().size(); i++){
//       FillTrueVertexes( evt.true_particles()[i] );
//     }

//     return;

//   }

//   //*************************************************************
//   void clusterizer::FillTrueVertexes( mybhep::particle* tp ){
//     //*************************************************************

//     if( tp->find_property("foil_vertex") )
//       tp->change_property("foil_vertex", "1");
//     else
//       tp->add_property("foil_vertex", "1");

//     return;

//   }


//   //*************************************************************
//   void clusterizer::GenerateWires( void ){
//     //*************************************************************
//     m.message("CAT::clusterizer::GenerateWires: Entering...",mybhep::NORMAL);

//     //clock.start(" clusterizer: generate wires ");

//     DriftWires.clear();

//     if( SuperNemo )
//       {
//         m.message("CAT::clusterizer::GenerateWires: SuperNemo geometry...",mybhep::NORMAL);

//         // TRACKING GG BLOCKS
//         double theta = M_PI/8.;
//         m.message("CAT::clusterizer::GenerateWires: rad=",
//                   rad,
//                   mybhep::NORMAL);
//         if (std::isnan (GG_CELL_pitch))
//           {
//             GG_CELL_pitch = rad*cos(theta);
//           }
//         m.message("CAT::clusterizer::GenerateWires: GG_CELL_pitch=",
//                   GG_CELL_pitch,
//                   mybhep::NORMAL);
//         m.message("CAT::clusterizer::GenerateWires: GG_GRND_diam=",
//                   GG_GRND_diam,
//                   mybhep::NORMAL);
//         m.message("CAT::clusterizer::GenerateWires: CHAMBER_X=",
//                   CHAMBER_X,
//                   mybhep::NORMAL);
//         if (num_cells_per_plane <= 0)
//           {
//             num_cells_per_plane=(int)((CHAMBER_X-GG_GRND_diam)/GG_CELL_pitch);
//           }
//         m.message("CAT::clusterizer::GenerateWires: num_cells_per_plane=",
//                   num_cells_per_plane,
//                   mybhep::NORMAL);
//         GG_BLOCK_X = num_cells_per_plane*GG_CELL_pitch+GG_GRND_diam;
//         m.message("CAT::clusterizer::GenerateWires: GG_BLOCK_X=",
//                   GG_BLOCK_X,
//                   mybhep::NORMAL);

//         std::vector<double> GG_BLOCK_thick;
//         std::vector<double> GG_BLOCK_posZ;
//         GG_BLOCK_thick.assign (num_blocks, std::numeric_limits<double>::quiet_NaN ());
//         GG_BLOCK_posZ.assign (num_blocks, std::numeric_limits<double>::quiet_NaN ());
//         double distance = SOURCE_thick/2.;
//         m.message("CAT::clusterizer::GenerateWires: SOURCE_thick=",
//                   SOURCE_thick,
//                   mybhep::NORMAL);
//         //Calculating thickness and positions of blocks
//         for(int i=0; i<num_blocks; i++){
//           GG_BLOCK_thick[i]=(planes_per_block[i]*GG_CELL_pitch+GG_GRND_diam);
//           GG_BLOCK_posZ[i] = distance + gaps_Z[i] + GG_BLOCK_thick[i]/2.;
//           distance = GG_BLOCK_posZ[i] + GG_BLOCK_thick[i]/2.;
//         }
//         m.message("CAT::clusterizer::GenerateWires: DEVEL: **** STEP 1",mybhep::NORMAL);

//         int sign[2];
//         sign[0] = 1;
//         sign[1] = -1;

//         for(size_t isign=0; isign<2; isign++)
//           for (int iblock=0; iblock<num_blocks;iblock++){

//             double block_pos[3];
//             block_pos[2] = GG_BLOCK_posZ[iblock];

//             // loop over planes in block

//             double plane_pos_z0 = block_pos[2]-GG_BLOCK_thick[iblock]/2;

//             for(int iplane=0; iplane<planes_per_block[iblock];iplane++){

//               double plane_pos_Z;
//               plane_pos_Z=GG_GRND_diam/2.+GG_CELL_pitch/2.+iplane*GG_CELL_pitch;
//               plane_pos_Z = plane_pos_Z + plane_pos_z0;
//               plane_pos_Z *= sign[isign];

//               for(int iwire=0; iwire<num_cells_per_plane; iwire++)
//                 {
//                   double increment = GG_CELL_pitch*((double)iwire);
//                   double xpos = -(CHAMBER_X-GG_GRND_diam)/2.+6.+increment+GG_CELL_pitch/2.;

//                   POINT point;
//                   point.x = xpos;
//                   point.z = plane_pos_Z;
//                   DriftWires.push_back( point );
//                 }


//             }//end of loop over planes

//           }//end of loop over blocks
//       }
//     else
//       {

//         size_t NOfWires[18];
//         double FirstWireX[18];
//         double FirstWireZ[18];
//         double FirstWirePhi[18];
//         double LayerRadius[18];

//         NOfWires[0] = 320;
//         NOfWires[1] = 320;
//         NOfWires[2] = 320;
//         NOfWires[3] = 320;
//         NOfWires[4] = 280;
//         NOfWires[5] = 280;
//         NOfWires[6] = 240;
//         NOfWires[7] = 240;
//         NOfWires[8] = 240;
//         NOfWires[9] = 360;
//         NOfWires[10] = 360;
//         NOfWires[11] = 360;
//         NOfWires[12] = 360;
//         NOfWires[13] = 400;
//         NOfWires[14] = 400;
//         NOfWires[15] = 460;
//         NOfWires[16] = 460;
//         NOfWires[17] = 460;

//         FirstWireX[0] = 1521.;
//         FirstWireX[1] = 1493.;
//         FirstWireX[2] = 1465.;
//         FirstWireX[3] = 1437.;
//         FirstWireX[4] = 1271.;
//         FirstWireX[5] = 1243.;
//         FirstWireX[6] = 1077.;
//         FirstWireX[7] = 1049.;
//         FirstWireX[8] = 1021.;
//         FirstWireX[9] = 1579.;
//         FirstWireX[10] = 1607.;
//         FirstWireX[11] = 1635.;
//         FirstWireX[12] = 1663.;
//         FirstWireX[13] = 1829.;
//         FirstWireX[14] = 1857.;
//         FirstWireX[15] = 2023.;
//         FirstWireX[16] = 2051.;
//         FirstWireX[17] = 2079.;

//         FirstWireZ[0] = 10.57;
//         FirstWireZ[1] = 10.375;
//         FirstWireZ[2] = 10.181;
//         FirstWireZ[3] = 9.9871;
//         FirstWireZ[4] = 10.167;
//         FirstWireZ[5] = 9.943;
//         FirstWireZ[6] = 10.180;
//         FirstWireZ[7] = 9.915;
//         FirstWireZ[8] = 9.651;
//         FirstWireZ[9] = 9.811;
//         FirstWireZ[10] = 9.985;
//         FirstWireZ[11] = 10.158;
//         FirstWireZ[12] = 10.333;
//         FirstWireZ[13] = 10.129;
//         FirstWireZ[14] = 10.284;
//         FirstWireZ[15] = 9.804;
//         FirstWireZ[16] = 9.939;
//         FirstWireZ[17] = 10.075;


//         for(size_t i=0; i<18; i++)
//           {
//             LayerRadius[i] = sqrt(FirstWireX[i]*FirstWireX[i] + FirstWireZ[i]*FirstWireZ[i]);
//             FirstWirePhi[i] = acos(FirstWireX[i]/LayerRadius[i]);

//             for(size_t j=0; j<NOfWires[i]; j++)
//               {
//                 double layerphi = 2.*M_PI/NOfWires[i];
//                 double ph = FirstWirePhi[i] + j*layerphi;

//                 POINT point;
//                 point.x = LayerRadius[i]*cos(ph);
//                 point.z = LayerRadius[i]*sin(ph);

//                 DriftWires.push_back( point );
//               }
//           }
//       }


//     //clock.stop(" clusterizer: generate wires ");
//     m.message("CAT::clusterizer::GenerateWires: Done.",mybhep::NORMAL);

//     return;

//   }


//   //*******************************************************************
//   size_t clusterizer::get_true_hit_index(mybhep::hit& hit, bool print){
//     //*******************************************************************

//     topology::node tn(hit, 0, SuperNemo, level, probmin);

//     for(std::vector<topology::cell>::iterator ic=cells_.begin(); ic!=cells_.end(); ++ic){
//       if( ic->same_cell(tn.c()) )
//         return ic->id();
//     }

//     if( print )
//       m.message("CAT::clusterizer::get_true_hit_index: warning: can't find corresponding reco hit for true hit (", tn.c().ep().x().value(), ", ", tn.c().ep().y().value(), ", ", tn.c().ep().z().value(), ") layer", tn.c().layer(), mybhep::VVERBOSE);

//     return 0;

//   }

//   //*******************************************************************
//   size_t clusterizer::get_nemo_hit_index(mybhep::hit& hit, bool print){
//     //*******************************************************************

//     topology::node tn(hit, 0, SuperNemo, level, probmin);

//     for(std::vector<topology::cell>::iterator ic=cells_.begin(); ic!=cells_.end(); ++ic){
//       if( ic->same_cell(tn.c()) )
//         return ic->id();
//     }

//     if( print )
//       m.message("CAT::clusterizer::get_nemo_hit_index: warning: can't find corresponding reco hit for nemo hit (", tn.c().ep().x().value(), ", ", tn.c().ep().y().value(), ", ", tn.c().ep().z().value(), ") layer", tn.c().layer(), mybhep::VVERBOSE);

//     return 0;

//   }

//   //*******************************************************************
//   size_t clusterizer::get_calo_hit_index(const topology::calorimeter_hit & c){
//     //*******************************************************************

//     for(std::vector<topology::calorimeter_hit>::iterator ic=calorimeter_hits_.begin(); ic!=calorimeter_hits_.end(); ++ic){
//       if( ic->same_calo(c) )
//         return ic->id();
//     }

//     m.message("CAT::clusterizer::get_calo_hit_index: warning: can't find corresponding calo hit for nemo calo hit (", c.pl().center().x().value(), ", ", c.pl().center().y().value(), ", ", c.pl().center().z().value(), ") layer", c.layer(), mybhep::VVERBOSE);

//     return 0;

//   }


//   //*******************************************************************
//   bool clusterizer::read_event(mybhep::event& event_ref, topology::tracked_data & tracked_data_){
//     //*******************************************************************

//     clock.start(" clusterizer: read event ","cumulative");

//     m.message("CAT::clusterizer::read_event: local_tracking: reading event", mybhep::VERBOSE);

//     cells_.clear();
//     parts.clear();
//     clusters_.clear();
//     calorimeter_hits_.clear();
//     true_sequences_.clear();
//     nemo_sequences_.clear();

//     bool bhep_input = true;

//     if( bhep_input ){

//       // read digi particles
//       std::vector<mybhep::particle*> digi_parts;

//       if( !SuperNemo && first_event)
//         {
//           first_event=false;

//           N3_MC = true;
//           if( event_ref.find_property("DATA") )
//             if( mybhep::int_from_string(event_ref.fetch_property("DATA")) == 1 )
//               N3_MC = false;  // data

//           if( N3_MC )
//             m.message("CAT::clusterizer::read_event: Nemo3 MC ", mybhep::NORMAL);
//           else
//             m.message("CAT::clusterizer::read_event: Nemo3 data ", mybhep::NORMAL);
//         }

//       FillYPositions(event_ref);

//       if( SuperNemo ){
//         FillTrueVertexes(event_ref);
//       }

//       fill_fast_information(event_ref);

//       event_ref.filter(mybhep::DIGI,"SUNAMI","1",digi_parts); // scans all event.digi_particles(),
//       // and copies to "digi_parts" those ones having property "SUNAMI" = 1

//       if ( digi_parts.empty() )
//         {
//           m.message("CAT::clusterizer::read_event: Problem: there are no gg hits in this event.\n",mybhep::MUTE);

//           return false;
//         }

//       if (digi_parts[0]->find_property("mID"))               // BAR design??
//         {
//           for (size_t pnr=0; pnr < digi_parts.size(); pnr++)
//             if (digi_parts[pnr]->fetch_property("mID") == _moduleNR )
//               parts.push_back(digi_parts[pnr]);
//         }
//       else                                               //NOT BAR DESIGN and OLD data
//         {
//           //    for (size_t pnr=0; pnr < digi_parts.size(); pnr++)
//           //      parts.push_back(digi_parts[pnr]);
//           parts.push_back(digi_parts[0]);
//         }

//       if( parts.size() != 1 )
//         return false;

//       const std::vector<mybhep::hit*>& hits = parts[0]->hits("trk");

//       for (size_t ihit=0; ihit<hits.size();ihit++){
//         topology::cell c(*hits[ihit],ihit, SuperNemo, level, probmin);
//         c.set_small_radius(SmallRadius);
//         cells_.push_back(c);
//       }

//       clock.start(" clusterizer: make calo hit ","cumulative");
//       const std::vector<mybhep::hit*>& chits = parts[0]->hits("cal");
//       for (size_t ihit=0; ihit<chits.size();ihit++){
//         topology::calorimeter_hit ch = make_calo_hit(*chits[ihit], ihit);
//         calorimeter_hits_.push_back(ch);
//       }
//       clock.stop(" clusterizer: make calo hit ");

//       if( level >= mybhep::VVERBOSE )
//         print_calos();



//       read_true_sequences(event_ref);

//       if( !SuperNemo )
//         read_nemo_sequences(event_ref);

//     }


//     order_cells();

//     setup_cells();

//     tracked_data_.set_cells(cells_);

//     tracked_data_.set_calos(calorimeter_hits_);

//     tracked_data_.set_true_sequences(true_sequences_);

//     tracked_data_.set_nemo_sequences(nemo_sequences_);

//     clock.stop(" clusterizer: read event ");


//     return true;


//   }


//   //*******************************************************************
//   bool clusterizer::prepare_event(topology::tracked_data & tracked_data_){
//     //*******************************************************************

//     clock.start(" clusterizer: prepare event ","cumulative");

//     event_number ++;
//     m.message("CAT::clusterizer::prepare_event: local_tracking: preparing event", event_number, mybhep::VERBOSE);

//     if( event_number < first_event_number ){
//       m.message("CAT::clusterizer::prepare_event: local_tracking: skip event", event_number, " first event is ", first_event_number,  mybhep::VERBOSE);
//       return false;
//     }

//     parts.clear();
//     clusters_.clear();

//     order_cells();
//     setup_cells();

//     tracked_data_.set_cells(cells_);
//     tracked_data_.set_calos(calorimeter_hits_);

//     clock.stop(" clusterizer: prepare event ");


//     return true;


//   }


//   //*******************************************************************
//   void clusterizer::read_true_sequences(mybhep::event& event_ref){
//     //*******************************************************************
//     const std::vector<mybhep::particle*>& truep = event_ref.true_particles();
//     m.message("CAT::clusterizer::read_true_sequences: number of true particles is", truep.size(), mybhep::VVERBOSE);
//     if (truep.size()!=0){

//       for( size_t ipart=0; ipart<truep.size(); ipart++)
//         {

//           mybhep::particle& tp = *truep[ipart];
//           const std::vector<mybhep::hit*>& thits = tp.hits("trk");
//           topology::sequence trueseq;
//           std::vector<topology::node> truenodes;
//           for(size_t i=0; i<thits.size(); i++)
//             {
//               size_t index = get_true_hit_index(*thits[i], tp.primary());
//               topology::node tn(*thits[i], index, SuperNemo, level, probmin);
//               truenodes.push_back(tn);
//             }
//           trueseq.set_nodes(truenodes);
//           if( tp.find_property("charge"))
//             trueseq.set_charge(topology::experimental_double(tp.charge(),0.));//particle charge
//           else if( tp.name() == "e+")
//             trueseq.set_charge(topology::experimental_double(1.,0.));
//           else if( tp.name() == "e-")
//             trueseq.set_charge(topology::experimental_double(-1.,0.));
//           else
//             trueseq.set_charge(topology::experimental_double(-1.,0.));

//           if( tp.find_property("length"))
//             trueseq.set_helix_length(topology::experimental_double(mybhep::double_from_string(tp.fetch_property("length")), 0.));

//           trueseq.set_helix_vertex( topology::experimental_point( tp.vertex().x(), tp.vertex().y(), tp.vertex().z(), 0., 0., 0.),"true");

//           size_t cindex = 0;
//           const std::vector<mybhep::hit*>& chits = tp.hits("cal");
//           if( !chits.empty() ){
//             topology::calorimeter_hit ch = make_calo_hit(*chits[0], 0);
//             cindex = get_calo_hit_index(ch);
//             trueseq.set_decay_helix_vertex( topology::experimental_point( tp.decay_vertex().x(), tp.decay_vertex().y(), tp.decay_vertex().z(), 0., 0., 0.),"calo", cindex);
//           }

//           trueseq.set_name( tp.name() );

//           trueseq.set_momentum( topology::experimental_vector(tp.p3().x(), tp.p3().y(), tp.p3().z(), 0., 0., 0.));

//           trueseq.set_primary( tp.primary());

//           true_sequences_.push_back(trueseq);



//         }

//       if( level >= mybhep::VVERBOSE ){
//         print_true_sequences();
//       }

//     }

//     return;

//   }

//   //*******************************************************************
//   void clusterizer::read_nemo_sequences(mybhep::event& event_ref){
//     //*******************************************************************

//     std::vector<mybhep::particle*> nemo_parts;

//     event_ref.filter(mybhep::DIGI,"NEMO","1",nemo_parts);

//     if (! nemo_parts.empty()){
//       m.message("CAT::clusterizer::read_nemo_sequences: number of nemo particles is", nemo_parts.size(), mybhep::VVERBOSE);

//       for( size_t ipart=0; ipart<nemo_parts.size(); ipart++)
//         {

//           mybhep::particle& tp = *nemo_parts[ipart];
//           const std::vector<mybhep::hit*>& thits = tp.hits("trk");
//           topology::sequence nemoseq;
//           std::vector<topology::node> nemonodes;
//           for(size_t i=0; i<thits.size(); i++)
//             {
//               size_t index = get_nemo_hit_index(*thits[i], tp.primary());
//               topology::node tn(*thits[i], index, SuperNemo, level, probmin);
//               nemonodes.push_back(tn);
//             }
//           nemoseq.set_nodes(nemonodes);
//           if( tp.find_property("charge"))
//             nemoseq.set_charge(topology::experimental_double(tp.charge(),0.));//particle charge

//           if( tp.find_property("length"))
//             nemoseq.set_helix_length(topology::experimental_double(mybhep::double_from_string(tp.fetch_property("length")), 0.));

//           nemoseq.set_helix_vertex( topology::experimental_point( tp.vertex().x(), tp.vertex().y(), tp.vertex().z(), 0., 0., 0.),"nemo");


//           size_t cindex = 0;
//           const std::vector<mybhep::hit*>& chits = tp.hits("cal");
//           if( !chits.empty() ){
//             topology::calorimeter_hit ch = make_calo_hit(*chits[0], 0);
//             cindex = get_calo_hit_index(ch);
//             nemoseq.set_decay_helix_vertex( topology::experimental_point( tp.decay_vertex().x(), tp.decay_vertex().y(), tp.decay_vertex().z(), 0., 0., 0.),"calo", cindex);
//           }

//           nemoseq.set_name( tp.name() );

//           nemoseq.set_momentum( topology::experimental_vector(tp.p3().x(), tp.p3().y(), tp.p3().z(), 0., 0., 0.));

//           nemoseq.set_primary( tp.primary());

//           nemo_sequences_.push_back(nemoseq);



//         }

//       if( level >= mybhep::VVERBOSE ){
//         print_nemo_sequences();
//       }

//     }

//     return;

//   }

//   //*******************************************************************
//   void clusterizer::print_cells(void)const{
//     //*******************************************************************

//     for(std::vector<topology::cell>::const_iterator icell=cells_.begin(); icell!=cells_.end();++icell){
//       icell->dump();
//     }

//     return;
//   }


//   //*******************************************************************
//   void clusterizer::print_calos(void)const{
//     //*******************************************************************

//     for(std::vector<topology::calorimeter_hit>::const_iterator icalo=calorimeter_hits_.begin(); icalo!=calorimeter_hits_.end();++icalo){
//       icalo->dump();
//     }

//     return;
//   }

//   //*******************************************************************
//   void clusterizer::print_true_sequences(void)const{
//     //*******************************************************************

//     for(std::vector<topology::sequence>::const_iterator iseq=true_sequences_.begin(); iseq != true_sequences_.end(); ++iseq){
//       iseq->dump();
//     }

//     return;
//   }


//   //*******************************************************************
//   void clusterizer::print_nemo_sequences(void)const{
//     //*******************************************************************

//     for(std::vector<topology::sequence>::const_iterator iseq=nemo_sequences_.begin(); iseq != nemo_sequences_.end(); ++iseq){
//       iseq->dump();
//     }

//     return;
//   }



//   //*******************************************************************
//   void clusterizer::print_clusters(void)const{
//     //*******************************************************************

//     for(std::vector<topology::cluster>::const_iterator icluster=clusters_.begin(); icluster != clusters_.end(); ++icluster){
//       icluster->dump();
//     }

//     return;
//   }



//   //*******************************************************************
//   void clusterizer::clusterize(topology::tracked_data & tracked_data_){
//     //*******************************************************************

//     if( cells_.empty() ) return;

//     if( !select_true_tracks(tracked_data_) ){
//       m.message("CAT::clusterizer::clusterize: event is not selected at true level ", mybhep::NORMAL);
//       tracked_data_.set_selected(false);
//       return;
//     }

//     float side[2]; // loop on two sides of the foil
//     side[0] =  1.;
//     side[1] = -1.;

//     bool fast[2]; // loop on fast and slow hits
//     fast[0] = true;
//     fast[1] = false;

//     std::map<int,unsigned int > flags;

  //  Plus nécessaire car le travail est fait par le preclustering
//     for(size_t ip=0; ip<2; ip++)  // loop on two sides of the foil
//       {
//         for(size_t iq=0; iq<2; iq++) // loop on fast and slow hits
//           {
//             for(size_t i=0; i<cells_.size(); i++)
//               {
//                 flags[cells_[i].id()] = 0;
//               }

//             for(std::vector<topology::cell>::const_iterator icell=cells_.begin(); icell!=cells_.end(); ++icell){
//               // pick a cell c that was never added
//               const topology::cell & c = *icell;
//               if( (cell_side(c) * side[ip]) < 0) continue;
//               if( c.fast() != fast[iq] ) continue;
//               if( flags[c.id()] == 1 ) continue;
//               flags[c.id()] = 1;

//               // cell c will form a new cluster, i.e. a new list of nodes
//               topology::cluster cluster_connected_to_c;
//               std::vector<topology::node> nodes_connected_to_c;
//               m.message("CAT::clusterizer::clusterize: begin new cluster with cell ", c.id(), mybhep::VERBOSE);

//               // let's get the list of all the cells that can be reached from c
//               // without jumps
//               std::vector<topology::cell> cells_connected_to_c;
//               cells_connected_to_c.push_back(c);

//               for( size_t i=0; i<cells_connected_to_c.size(); i++){ // loop on connected cells
//                 // take a connected cell (the first one is just c)
//                 topology::cell cconn = cells_connected_to_c[i];

//                 // the connected cell composes a new node
//                 topology::node newnode(cconn, level, probmin);
//                 std::vector<topology::cell_couplet> cc;

//                 // get the list of cells near the connected cell
//                 std::vector<topology::cell> cells_near_iconn = get_near_cells(cconn);

//                 m.message("CAT::clusterizer::clusterize: cluster ", clusters_.size(), " starts with ", c.id(), " try to add cell ", cconn.id(), " with n of neighbours = ", cells_near_iconn.size(), mybhep::VERBOSE);
//                 for(std::vector<topology::cell>::const_iterator icnc=cells_near_iconn.begin(); icnc!=cells_near_iconn.end(); ++icnc){

//                   topology::cell cnc = *icnc;

//                   if( !is_good_couplet(& cconn, cnc, cells_near_iconn) ) continue;

//                   topology::cell_couplet ccnc(cconn,cnc,level,probmin);
//                   cc.push_back(ccnc);

//                   m.message("CAT::clusterizer::clusterize: ... creating couplet ", cconn.id(), " -> ", cnc.id(), mybhep::VERBOSE);

//                   if( flags[cnc.id()] != 1 )
//                     {
//                       flags[cnc.id()] = 1 ;
//                       cells_connected_to_c.push_back(cnc);
//                     }
//                 }
//                 newnode.set_cc(cc);
//                 newnode.calculate_triplets(Ratio, QuadrantAngle, TangentPhi, TangentTheta);
//                 nodes_connected_to_c.push_back(newnode);

//                 m.message("CAT::clusterizer::clusterize: cluster started with ", c.id(), " has been given cell ", cconn.id(), " with ", cc.size(), " couplets ", mybhep::VERBOSE);

//               }

//               cluster_connected_to_c.set_nodes(nodes_connected_to_c);

//               clusters_.push_back(cluster_connected_to_c);
//             }

//           }
//       }


//     setup_clusters();

//     m.message("CAT::clusterizer::clusterize: there are ", clusters_.size(), " clusters of cells ", mybhep::VVERBOSE);

//     if( PrintMode )
//       make_plots(tracked_data_);

//     if( level >= mybhep::VVERBOSE ){
//       print_clusters();
//     }

//     tracked_data_.set_cells(cells_);
//     tracked_data_.set_clusters(clusters_);

//     clock.stop(" clusterizer: clusterize ");


//     return;

//   }

//   //*******************************************************************
//   void clusterizer::clusterize_after_sultan(topology::tracked_data & tracked_data_){
//     //*******************************************************************

//     if( event_number < first_event_number ) return;

//     clock.start(" clusterizer: clusterize_after_sultan ","cumulative");

//     m.message("CAT::clusterizer::clusterize_after_sultan: local_tracking: fill clusters ", mybhep::VERBOSE);

//     if( cells_.empty() ) return;

//     if( !select_true_tracks(tracked_data_) ){
//       m.message("CAT::clusterizer::clusterize_after_sultan: event is not selected at true level ", mybhep::NORMAL);
//       tracked_data_.set_selected(false);
//       return;
//     }


//     m.message("CAT::clusterizer::clusterize_after_sultan: sultan has made ", tracked_data_.clusters_.size(), " clusters of cells ", mybhep::VVERBOSE);

//     // loop on sultan clusters
//     for(std::vector<topology::cluster>::iterator iclu=tracked_data_.clusters_.begin(); iclu!=tracked_data_.clusters_.end(); ++iclu){


//       m.message("CAT::clusterizer::clusterize_after_sultan: cluster [", iclu - tracked_data_.clusters_.begin(), "] has ", iclu->nodes_.size(), " cells ", mybhep::VVERBOSE);

//       if( level >= mybhep::VVERBOSE ){
// 	std::clog << "[";
// 	for( std::vector<topology::node>::iterator inode = iclu->nodes_.begin(); inode != iclu->nodes_.end(); ++inode )
// 	  std::clog << inode->c().id() << " ";
// 	std::clog << "]" << std::endl;
//       }

//       // loop on nodes in sultan cluster
//       for( std::vector<topology::node>::iterator inode = iclu->nodes_.begin(); inode != iclu->nodes_.end(); ++inode ){

// 	topology::cell c = inode->c();

// 	std::vector<topology::cell_couplet> cc;
// 	std::vector<topology::cell> links;
// 	if( inode - iclu->nodes_.begin() > 0 ){
// 	  topology::cell pc = iclu->nodes_[inode - iclu->nodes_.begin() - 1].c();
// 	  topology::cell_couplet pcc(c, pc);
// 	  m.message("CAT::clusterizer::clusterize_after_sultan: adding couplet", c.id(), "-", pc.id(), mybhep::VVERBOSE);
// 	  cc.push_back(pcc);
// 	}
// 	if( (size_t)(inode - iclu->nodes_.begin() + 1) < iclu->nodes_.size() ){
// 	  topology::cell nc = iclu->nodes_[inode - iclu->nodes_.begin() + 1].c();
// 	  topology::cell_couplet ncc(c,nc);
// 	  m.message("CAT::clusterizer::clusterize_after_sultan: adding couplet", c.id(), "-", nc.id(), mybhep::VVERBOSE);
// 	  cc.push_back(ncc);
// 	  links.push_back(nc);
// 	}

// 	m.message("CAT::clusterizer::clusterize_after_sultan: node [", inode->c().id(), "] has ", cc.size(), " couplets ", mybhep::VVERBOSE);

// 	inode->set_cc(cc);
// 	inode->set_links(links);
// 	inode->calculate_triplets_after_sultan(Ratio);
//       }
//     }

//     clusters_ = tracked_data_.clusters_;

//     setup_clusters();

//     if( PrintMode )
//       make_plots(tracked_data_);

//     if( level >= mybhep::VVERBOSE ){
//       print_clusters();
//     }

//     clock.stop(" clusterizer: clusterize_after_sultan ");


//     return;

//   }

//   //*************************************************************
//   void clusterizer::make_plots(topology::tracked_data & /* tracked_data_ */){
//     //*************************************************************
//     /*
//       if( PrintMode ){
//       for(std::vector<topology::cluster>::iterator icluster = clusters_.begin(); icluster!=clusters_.end(); ++icluster){
//       for(std::vector<topology::node>::iterator inode = (*icluster).nodes_.begin(); inode != (*icluster).nodes_.end(); ++inode){
//       for(std::vector<topology::cell_triplet>::iterator iccc = (*inode).ccc_.begin(); iccc != (*inode).ccc_.end(); ++iccc){
//       topology::cell_triplet ccc = *iccc;
//       for(std::vector<double>::const_iterator ichi = ccc.chi2s().begin(); ichi != ccc.chi2s().end(); ++ichi){
//       hman.fill("chi2_triplet", *ichi);
//       }
//       for(std::vector<double>::const_iterator iprob = ccc.probs().begin(); iprob != ccc.probs().end(); ++iprob){
//       hman.fill("prob_triplet", *iprob);
//       }
//       }
//       }

//       }

//       std::vector<topology::sequence> true_sequences = tracked_data_.get_true_sequences();
//       for(std::vector<topology::sequence>::iterator iseq=true_sequences.begin(); iseq != true_sequences.end(); ++iseq){
//       topology::node n;
//       double phi;
//       if( iseq->largest_kink_node(n, phi)){
//       phi *= 180./M_PI;
//       std::clog << " largest kink on true sequence " << iseq - true_sequences.begin() << " is " << phi << " on cell " << n.c().id() << std::endl;
//       hman.fill("largest_true_kink", phi);
//       if( phi > 20.){
//       topology::experimental_vector dist(n.c().ep(), n.ep());
//       hman.fill("largest_true_kink_position", dist.z().value(), dist.x().value());
//       }
//       }

//       }




//       }

//     */
//     return;


//   }

//   //*************************************************************
//   bool clusterizer::is_good_couplet(topology::cell * mainc,
//                                     const topology::cell &candidatec,
//                                     const std::vector<topology::cell> & nearmain){
//     //*************************************************************

//     // the couplet mainc -> candidatec is good only if
//     // there is no other cell that is near to both and can form a triplet between them

//     clock.start(" clusterizer: is good couplet ","cumulative");

//     topology::cell a=*mainc;


//     for(std::vector<topology::cell>::const_iterator icell=nearmain.begin(); icell != nearmain.end(); ++icell){

//       topology::cell b=*icell;
//       if( b.id() == candidatec.id()) continue;

//       if(near_level(b, candidatec) == 0 ) continue;

//       if(near_level(b, candidatec) < near_level(a, candidatec) ||
//          near_level(b, a) < near_level(a, candidatec) )
//         continue;  // cannot match a->b or b->c if a->c is nearer

//       //    if( icell->intersect(candidatec) || icell->intersect(mainc) ) continue;
//       // don't reject candidate based on a cell that intersects it

//       m.message("CAT::clusterizer::is_good_couplet: ... ... check if near node ", b.id(), " has triplet ", a.id(), " <-> ", candidatec.id(), mybhep::VERBOSE);

//       topology::cell_triplet ccc(a,b,candidatec, level, probmin);
//       ccc.calculate_joints(Ratio, QuadrantAngle, TangentPhi, TangentTheta);
//       if(ccc.joints().size() > 0 ){
//         m.message("CAT::clusterizer::is_good_couplet: ... ... yes it does: so couplet ", a.id(), " and ", candidatec.id(), " is not good",  mybhep::VERBOSE);
//         clock.stop(" clusterizer: is good couplet ");
//         return false;
//       }

//     }


//     clock.stop(" clusterizer: is good couplet ");
//     return true;

//   }


//   //*************************************************************
//   void clusterizer::fill_fast_information( mybhep::event& evt ){
//     //*************************************************************

//     for(size_t i=0; i<evt.digi_particles().size(); i++){
//       if( evt.digi_particles()[i]->find_property("SUNAMI") )
//         fill_fast_information( evt.digi_particles()[i] );
//     }


//     return;

//   }

//   //*************************************************************
//   void clusterizer::fill_fast_information( mybhep::particle* p ){
//     //*************************************************************

//     for(size_t i=0; i<p->hits("trk").size(); i++){
//       fill_fast_information(p->hits("trk")[i]);
//     }

//     return;

//   }

//   //*************************************************************
//   void clusterizer::fill_fast_information( mybhep::hit* h ){
//     //*************************************************************

//     double radius = mybhep::double_from_string(h->fetch_property("DIST"));
//     double tdelay = mybhep::double_from_string(h->fetch_property("ATIME"));
//     if( std::isnan(tdelay) || std::isinf(tdelay) )
//       tdelay = -1.;

//     if( radius != 0. && tdelay == 0. ) // fast hit
//       h->add_property("FAST",mybhep::to_string(1));
//     else if( radius == 0. && tdelay != 0. ){ // delayed hit
//       h->add_property("SLOW",mybhep::to_string(1));
//     }
//     else
//       {
//         m.message("CAT::clusterizer::fill_fast_information: Problem: hit has radius", radius, "delay time", tdelay, mybhep::NORMAL);
//       }

//     return;

//   }

//   //*************************************************************
//   int clusterizer::cell_side( const topology::cell & c){
//     //*************************************************************

//     if( SuperNemo )
//       {
//         if( c.ep().z().value() > 0. )
//           return 1;

//         return -1;
//       }


//     if( c.ep().radius().value() > FoilRadius )
//       return 1;

//     return -1;

//   }


//   size_t clusterizer::near_level( const topology::cell & c1, const topology::cell & c2 ){

//     // returns 0 for far-away cell
//     // 1 for diagonal cells
//     // 2 for side-by-side cells

//     // side-by-side connection: distance = 1
//     // diagonal connection: distance = sqrt(2) = 1.41
//     // skip 1 connection, side: distance = 2
//     // skip 1 connection, tilt: distance = sqrt(5) = 2.24
//     // skip 1 connection, diag: distance = 2 sqrt(2) = 2.83

//     topology::experimental_double distance = topology::experimental_vector(c1.ep(),c2.ep()).hor().length();

//     if( SuperNemo ){  // use side, layer and row

//       // Use geiger locator for such research Warning: use integer
//       // because uint32_t has strange behavior with absolute value
//       // cmath::abs
//       const int hit1_side  = c1.block();  // -1, 1
//       const int hit1_layer = abs(c1.layer()); // 0, 1, ..., 8
//       const int hit1_row   = c1.iid();  // -56, -55, ..., 55, 56

//       const int hit2_side  = c2.block();
//       const int hit2_layer = abs(c2.layer());
//       const int hit2_row   = c2.iid();

//       // Do not cross the foil
//       if (hit1_side != hit2_side) return 0;

//       // Check neighboring
//       const unsigned int layer_distance = abs (hit1_layer - hit2_layer); // 1 --> side-by-side
//       const unsigned int row_distance = abs (hit1_row - hit2_row);

//       if (layer_distance == 0 && row_distance == 0){
//         if( level >= mybhep::NORMAL ){
//           std::clog << "CAT::clusterizer::near_level: problem: cat asking near level of cells with identical posiion (" << hit1_side << ", " << hit1_layer << ", " << hit1_row << ") (" << hit2_side << ", " << hit2_layer << ", " << hit2_row << ")" << std::endl;
//         }
//         return 3;
//       }
//       else if (layer_distance == 1 && row_distance == 0) return 2;
//       else if (layer_distance == 0 && row_distance == 1) return 2;
//       else if (layer_distance == 1 && row_distance == 1) return 1;
//       return 0;

//     }else{ // use physical distance

//       double limit_side;
//       double limit_diagonal;
//       if (SuperNemo && SuperNemoChannel)
// 	{
// 	  limit_side = GG_CELL_pitch;
// 	  limit_diagonal = sqrt(2.)*GG_CELL_pitch;
// 	}
//       else
// 	{
// 	  double factor = cos(M_PI/8.); // 0.923879532511287 // octogonal factor = 0.92
// 	  limit_side = factor*CellDistance;
// 	  limit_diagonal = sqrt(2.)*factor*CellDistance; // new factor = 1.31
// 	}
//       double precision = 0.15*limit_side;

//       if( level >= mybhep::VVERBOSE )
// 	std::clog << "CAT::clusterizer::near_level: (c " << c2.id() << " d " << distance.value() << " )"
// 		  << std::endl;

//       if( std::abs(distance.value() - limit_side) < precision )
// 	return 2;

//       if( std::abs(distance.value() - limit_diagonal) < precision )
// 	return 1;

//       return 0;
//     }


//   }


//   std::vector<topology::cell> clusterizer::get_near_cells(const topology::cell & c){

//     clock.start(" clusterizer: get near cells ","cumulative");

//     m.message("CAT::clusterizer::get_near_cells: filling list of cells near cell ", c.id(), " fast ", c.fast(), " side ", cell_side(c), mybhep::VVERBOSE);

//     std::vector<topology::cell> cells;

//     for(std::vector<topology::cell>::iterator kcell=cells_.begin(); kcell != cells_.end(); ++kcell){
//       if( kcell->id() == c.id() ) continue;

//       if( kcell->fast() != c.fast() ) continue;

//       if( cell_side(*kcell) != cell_side(c) ) continue;

//       size_t nl = near_level(c,*kcell);

//       if( nl > 0 )
//         {
//           if( level >= mybhep::VVERBOSE ){
//             std::clog << "*";
//           }

//           topology::cell ck = *kcell;
//           cells.push_back(ck);
//         }
//     }

//     if( level >= mybhep::VVERBOSE )
//       std::clog << " " << std::endl;

//     clock.stop(" clusterizer: get near cells ");

//     return cells;

//   }


//   //*************************************************************
//   void clusterizer::setup_cells(){
//     //*************************************************************

//     for(std::vector<topology::cell>::iterator icell=cells_.begin(); icell!=cells_.end(); ++icell){
//       icell->set_print_level(level);
//       icell->set_probmin(probmin);
//     }

//     return;

//   }



//   //*************************************************************
//   void clusterizer::setup_clusters(){
//     //*************************************************************

//     clock.start(" clusterizer: setup_clusters ","cumulative");

//     // loop on clusters
//     for(std::vector<topology::cluster>::iterator icl=clusters_.begin(); icl != clusters_.end(); ++icl){
//       icl->set_print_level(level);
//       icl->set_probmin(probmin);

//       // loop on nodes
//       for(std::vector<topology::node>::iterator inode=(*icl).nodes_.begin(); inode != (*icl).nodes_.end(); ++inode){
//         inode->set_print_level(level);
//         inode->set_probmin(probmin);

//         for(std::vector<topology::cell_couplet>::iterator icc=(*inode).cc_.begin(); icc != (*inode).cc_.end(); ++icc){
//           icc->set_print_level(level);
//           icc->set_probmin(probmin);
//         }

//         for(std::vector<topology::cell_triplet>::iterator iccc=(*inode).ccc_.begin(); iccc != (*inode).ccc_.end(); ++iccc){
//           iccc->set_print_level(level);
//           iccc->set_probmin(probmin);
//         }

//       }

//     }

//     clock.stop(" clusterizer: setup_clusters ");

//     return;
//   }


//   //*************************************************************
//   int clusterizer::get_effective_layer(const mybhep::hit &hit){
//     //*************************************************************

//     // calculates effective layer of a hit in calorimeter block on +-X wall
//     // returns -8, -7, ..., -1, 0, 0, 1, ..., 7, 8
//     // returns N if the center of calo block is within layers N-1 and N
//     // for instance, if block center is between layers 3 and 4, returns 4

//     std::vector<double> block_pos;

//     //this deals with true hits
//     if (hit.find_property("Ini_Ekin")){
//       if(hit.find_property("BLK_Pos")){
//         mybhep::vector_from_string(hit.fetch_property("BLK_Pos"), block_pos);
//       }
//     }
//     //this deals with digi hits
//     else if (hit.find_property("E")){
//       std::string bp = hit.fetch_property("BLK_POS");
//       mybhep::vector_from_string(bp, block_pos);
//     }

//     double pos = GG_CELL_pitch;
//     bool found = false;
//     int counter = 0;

//     for(size_t i=0; i<planes_per_block.size(); i++){
//       if( found )
//         break;

//       pos += gaps_Z[i];

//       for(size_t j=0; j<(size_t)planes_per_block[i]; j++){

//         if( pos > std::abs(block_pos[2]) ){
//           found = true;
//           break;
//         }
//         pos += GG_CELL_pitch ;
//         counter ++;

//       }
//     }

//     int layer = (int)counter;
//     if( block_pos[2] < 0. )
//       layer *= -1;

//     return layer;

//   }

//   //*******************************************************************
//   topology::calorimeter_hit clusterizer::make_calo_hit(const mybhep::hit &ahit, size_t _id){
//     //*******************************************************************

//     std::vector<double> block_pos;
//     double en, time;
//     //this deals with true hits
//     if (ahit.find_property("Ini_Ekin")){
//       if(ahit.find_property("BLK_Pos")){
//         mybhep::vector_from_string(ahit.fetch_property("BLK_Pos"), block_pos);
//       }
//       en = mybhep::double_from_string(ahit.fetch_property("E_dep"));
//       time = mybhep::double_from_string(ahit.fetch_property("TOF"));
//     }
//     //this deals with digi hits
//     else if (ahit.find_property("E")){
//       std::string bp = ahit.fetch_property("BLK_POS");
//       mybhep::vector_from_string(bp, block_pos);
//       en = mybhep::double_from_string(ahit.fetch_property("E"));
//       time = mybhep::double_from_string(ahit.fetch_property("TIME"));
//     }

//     topology::experimental_point center(block_pos[0], block_pos[1], block_pos[2], 0., 0., 0.);

//     std::string block_type = ahit.fetch_property("BLK");
//     std::string plane;

//     topology::experimental_vector norm(0.,0.,0.,0.,0.,0.);

//     double layer=0;

//     double local_calo_X = calo_X;

//     if( SuperNemo )
//       {
//         double x_offset = 110.;
//         plane =  block_type.substr(0,11);

//         if (plane=="CALO_WRAP+X"){
//           norm.set_x(topology::experimental_double(-1.,0.));
//           layer = (double)get_effective_layer(ahit);
//           local_calo_X -= x_offset;
//         }
//         else  if (plane=="CALO_WRAP-X"){
//           norm.set_x(topology::experimental_double(1.,0.));
//           layer = (double)get_effective_layer(ahit);
//           local_calo_X -= x_offset;
//         }
//         else  if (plane=="CALO_WRAP-Z"){
//           norm.set_z(topology::experimental_double(1.,0.));
//         }
//         else  if (plane=="CALO_WRAP+Z"){
//           norm.set_z(topology::experimental_double(-1.,0.));
//         }
//         else{
//           m.message("CAT::clusterizer::make_calo_hit: problem: calo block not recognized",plane,mybhep::MUTE);
//           exit(1);
//         }
//       }
//     else
//       {
//         int block,planeid,id,idd;
//         sscanf(block_type.c_str(),"%d_%d_%d_%d",&block,&planeid,&id,&idd);



//         if (planeid==0 || planeid == 1){

//           topology::experimental_point origin(0.,0.,0.,0.,0.,0.);
//           norm = (topology::experimental_vector(origin, center)).hor().unit();

//           if (planeid==0){
//             //  plane = "inner"; //inner
//             layer = -(lastlayer - 1);
//           }
//           else  if (planeid==1){
//             plane = "outer"; //outer
//             //  *caloparam = OuterRadius;
//             layer = (lastlayer - 1);
//             norm *= -1.;
//           }

//           topology::experimental_vector mysizes(local_calo_X, calo_Y, calo_Z,
//                                                 0., 0., 0.);
//           double offset = std::abs((norm * mysizes).value())/2.;
//           center = (topology::experimental_vector(center) - norm * offset).point_from_vector();
//         }
//         else if( planeid==2 || planeid ==3){  // bottom and top, for Nemo3
//           if(planeid==2){
//             norm.set_y(topology::experimental_double(1.,0.));
//             //    center = topology::experimental_point(0,-ysize/2.,0.,0.,0.,0.);
//           }
//           else{
//             norm.set_y(topology::experimental_double(-1.,0.));
//             //    center = topology::experimental_point(0,ysize/2.,0.,0.,0.,0.);
//           }
//           if( id==0 ) layer=-5.5;
//           else if(id==1) layer=-3.5;
//           else if(id==2) layer=3.5;
//           else if(id==3) layer=5.5;
//           else{
//             m.message("CAT::clusterizer::make_calo_hit: CAL wall not recognized",block_type,mybhep::NORMAL);
//             exit(1);
//           }
//         }
//         else{
//           m.message("CAT::clusterizer::make_calo_hit: CAL wall not recognized",block_type,mybhep::NORMAL);
//           exit(1);
//         }
//       }

//     topology::experimental_vector sizes(local_calo_X, calo_Y, calo_Z,
//                                         0., 0., 0.);

//     double time_resol = 0.250; // ns
//     double e_resol = 0.240; // MeV
//     topology::experimental_double e(en, e_resol);
//     topology::experimental_double t(time, time_resol);
//     topology::plane pl(center, sizes, norm, level, probmin);
//     topology::calorimeter_hit ch(pl, e, t, _id, layer, level, probmin);

//     return ch;
//   }

//   //*************************************************************
//   void clusterizer::order_cells(){
//     //*************************************************************

//     clock.start(" clusterizer: order cells ","cumulative");

//     if( cells_.size() ){
//       if( level >= mybhep::VVERBOSE ){
//         std::clog << "CAT::clusterizer::order_cells: printing cells " << cells_.size() << std::endl;
//         print_cells();
//         std::clog << "CAT::clusterizer::order_cells: sorting cells " << std::endl;
//       }

//       //  std::sort( cells_.begin(), cells_.end(), topology::cell::compare );
//       std::sort( cells_.begin(), cells_.end());
//     }

//     clock.stop(" clusterizer: order cells ");

//     return;

//   }

//   //*************************************************************
//   bool clusterizer::select_true_tracks(topology::tracked_data & tracked_data_){
//     //*************************************************************

//     return true;

//     std::vector<topology::sequence> true_sequences = tracked_data_.get_true_sequences();

//     m.message("CAT::clusterizer::select_true_tracks: selecting events based on true tracks ", mybhep::VVERBOSE);

//     size_t counter = 0;
//     topology::sequence sa;
//     topology::sequence sb;

//     size_t n_tracks_ = 1;
//     double emin = 0.2;

//     for(std::vector<topology::sequence>::iterator iseq=true_sequences.begin(); iseq != true_sequences.end(); ++iseq){
//       if( !iseq->primary() ) continue;
//       if( counter == 0 )
//         sa = *iseq;
//       if( counter == 1 )
//         sb = *iseq;
//       counter ++;
//     }


//     if( counter != n_tracks_ ){
//       m.message("CAT::clusterizer::select_true_tracks: reject: n primary tracks = ", counter, mybhep::NORMAL);
//       return false;
//     }
//     m.message("CAT::clusterizer::select_true_tracks: n primary tracks = ", counter, mybhep::VVERBOSE);

//     if( !sa.one_side() ){
//       m.message("CAT::clusterizer::select_true_tracks: reject: 1st track crosses the foil ", mybhep::NORMAL);
//       return false;
//     }

//     if( n_tracks_ == 2 ){
//       if( !sb.one_side() ){
//         m.message("CAT::clusterizer::select_true_tracks: reject: 2nd track crosses the foil ", mybhep::NORMAL);
//         return false;
//       }
//     }

//     if( !sa.has_decay_helix_vertex() ){
//       m.message("CAT::clusterizer::select_true_tracks: reject: 1st track has no calo ", mybhep::NORMAL);
//       return false;
//     }


//     if( n_tracks_ == 2 ){
//       if( !sb.has_decay_helix_vertex() ){
//         m.message("CAT::clusterizer::select_true_tracks: reject: 2nd track has no calo ", mybhep::NORMAL);
//         return false;
//       }
//     }

//     if( n_tracks_ == 2 ){
//       if( sa.calo_helix_id() == sb.calo_helix_id() ){
//         m.message("CAT::clusterizer::select_true_tracks: reject: same calo ", sa.calo_helix_id(), mybhep::NORMAL);
//         return false;
//       }
//     }

//     std::vector<topology::calorimeter_hit> calos = tracked_data_.get_calos();

//     topology::calorimeter_hit caloA = calos[sa.calo_helix_id()];


//     if( caloA.e().value() < emin ){
//       m.message("CAT::clusterizer::select_true_tracks: reject: 1st calo has energy ", caloA.e().value(), mybhep::NORMAL);
//       return false;

//     }

//     if( n_tracks_ == 2 ){
//       topology::calorimeter_hit caloB = calos[sb.calo_helix_id()];
//       if( caloB.e().value() < emin ){
//         m.message("CAT::clusterizer::select_true_tracks: reject: 2nd calo has energy ", caloA.e().value(), mybhep::NORMAL);
//         return false;
//       }
//     }

//     /*
//       if( sa.nodes().size() > 15 ){
//       m.message(" reject: 1st tracks has nhits = ", sa.nodes().size(), mybhep::NORMAL);
//       return false;
//       }

//       if( sb.nodes().size() > 15 ){
//       m.message(" reject: 2nd tracks has nhits = ", sb.nodes().size(), mybhep::NORMAL);
//       return false;
//       }
//     */

//     size_t calo_counter = 0;
//     size_t calo_counter_thr = 0;

//     for(std::vector<topology::calorimeter_hit>::iterator iseq=calos.begin(); iseq != calos.end(); ++iseq){
//       calo_counter ++;
//       if( iseq->e().value() <= emin ) continue;
//       calo_counter_thr ++;
//     }

//     if( calo_counter != n_tracks_ ){
//       m.message("CAT::clusterizer::select_true_tracks: reject: n calos is ", calo_counter, mybhep::NORMAL);
//       return false;
//     }

//     if( calo_counter_thr != n_tracks_ ){
//       m.message("CAT::clusterizer::select_true_tracks: reject: n calos with E > ", emin, " is ", calo_counter_thr, mybhep::NORMAL);
//       return false;
//     }

// #if 0 // cut on kinks

//     size_t all_kinks_counter = 0;
//     for(std::vector<topology::sequence>::iterator iseq=true_sequences.begin(); iseq != true_sequences.end(); ++iseq){
//       if( !iseq->primary() ) continue;
//       topology::node n;
//       double phi;
//       size_t kinks_counter = 0;
//       int layer = -100;
//       if( iseq->largest_kink_node(n, phi)){
//         phi *= 180./M_PI;
//         if( phi > 20.){
//           kinks_counter ++;
//           all_kinks_counter ++;
//           layer = n.c().layer();
//         }
//       }

//       if( kinks_counter > 1 ){
//         m.message("CAT::clusterizer::select_true_tracks: reject: sequence ", iseq - true_sequences_.begin(), " has ", kinks_counter, " kinks ", mybhep::NORMAL);
//         return false;
//       }

//       if( kinks_counter > 0 && (layer == 0 || layer == 8 || layer == -8) ){
//         m.message("CAT::clusterizer::select_true_tracks: reject: sequence ", iseq - true_sequences_.begin(), " has kink on plane ", layer, mybhep::NORMAL);
//         return false;
//       }

//     }

//     if( all_kinks_counter == 0 ){
//       m.message("CAT::clusterizer::select_true_tracks: reject: no kinks in the event ", mybhep::NORMAL);
//       return false;
//     }


// #endif

//     return true;

//   }

//   void clusterizer::setDoDriftWires(bool ddw){
//     doDriftWires=ddw;
//     return;
//   }

//   void clusterizer::compute_lastlayer(){
//     lastlayer = 0;
//     for(size_t i=0; i<planes_per_block.size(); i++){
//       lastlayer += (int)planes_per_block[i];
//     }
//     return;
//   }

//   void clusterizer::set_GG_GRND_diam (double ggd){
//     GG_GRND_diam = ggd;
//     return;
//   }

//   void clusterizer::set_GG_CELL_diam (double ggcd){
//     GG_CELL_diam = ggcd;
//     return;
//   }

//   void clusterizer::set_lastlayer(int ll_){
//     lastlayer = ll_;
//     return;
//   }

//   void clusterizer::set_num_blocks(int nb){
//     if (nb > 0)
//       {
//         num_blocks = nb;
//         planes_per_block.assign (num_blocks, 1);
//       }
//     else
//       {
//         std::cerr << "WARNING: CAT::clusterizer::set_num_blocks: "
//                   << "Invalid number of GG layer blocks !" << std::endl;
//         planes_per_block.clear ();
//         num_blocks = -1; // invalid value
//       }
//     return;
//   }

//   void clusterizer::set_planes_per_block(int block, int nplanes){
//     if (block< 0 || block>= (int)planes_per_block.size())
//       {
//         throw std::range_error ("CAT::clusterizer::set_planes_per_block: Invalid GG layer block index !");
//       }
//     if (nplanes > 0)
//       {
//         planes_per_block.at (block) = nplanes;
//       }
//     else
//       {
//         throw std::range_error ("CAT::clusterizer::set_planes_per_block: Invalid number of GG layers in block !");
//       }
//     return;
//   }

//   void clusterizer::set_num_cells_per_plane(int ncpp){
//     if (ncpp <= 0)
//       {
//         num_cells_per_plane = -1; // invalid value
//       }
//     else
//       {
//         num_cells_per_plane = ncpp;
//       }
//     return;
//   }

//   void clusterizer::set_SOURCE_thick(double st){
//     if (st <= 0.0)
//       {
//         SOURCE_thick = std::numeric_limits<double>::quiet_NaN ();
//       }
//     else
//       {
//         SOURCE_thick = st;
//       }
//     return;
//   }

//   // What is it ?
//   void clusterizer::set_module_nr(const std::string & mID){
//     _moduleNR=mID;
//     return;
//   }

//   // What is it ?
//   int clusterizer::get_module_nr(void){
//     return _MaxBlockSize;
//   }

//   void clusterizer::set_MaxBlockSize(int mbs){
//     _MaxBlockSize=mbs;
//     return;
//   }

//   void clusterizer::set_pmax(double v){
//     if ( v <= 0.0)
//       {
//         pmax = std::numeric_limits<double>::quiet_NaN ();
//       }
//     else
//       {
//         pmax = v;
//       }
//     return;
//   }

//   void clusterizer::set_MaxTime(double v){
//     MaxTime = v;
//     return;
//   }

//   void clusterizer::set_PrintMode(bool v){
//     PrintMode = v;
//     return;
//   }

//   void clusterizer::set_SmallRadius(double v){
//     SmallRadius = v;
//     return;
//   }

//   void clusterizer::set_TangentPhi(double v){
//     TangentPhi = v;
//     return;
//   }

//   void clusterizer::set_TangentTheta(double v){
//     TangentTheta = v;
//     return;
//   }

//   void clusterizer::set_SmallNumber(double v){
//     SmallNumber = v;
//     return;
//   }

//   void clusterizer::set_QuadrantAngle(double v){
//     QuadrantAngle = v;
//     return;
//   }

//   void clusterizer::set_Ratio(double v){
//     Ratio = v;
//     return;
//   }

//   void clusterizer::set_CompatibilityDistance(double v){
//     CompatibilityDistance = v;
//     return;
//   }

//   void clusterizer::set_MaxChi2(double v){
//     MaxChi2 = v;
//     return;
//   }

//   void clusterizer::set_hfile(std::string v){
//     hfile = v;
//     return;
//   }

//   void clusterizer::set_probmin(double v){
//     probmin = v;
//     return;
//   }

//   void clusterizer::set_nofflayers(size_t v){
//     nofflayers = v;
//     return;
//   }

//   void clusterizer::set_first_event(int v){
//     first_event_number = v;
//     return;
//   }

//   void clusterizer::set_level(std::string v){
//     level = mybhep::get_info_level(v);
//     m = mybhep::messenger(level);
//     return;
//   }

//   void clusterizer::set_len(double v){
//     len = v;
//     return;
//   }

//   void clusterizer::set_vel(double v){
//     vel = v;
//     return;
//   }

//   void clusterizer::set_rad(double v){
//     rad = v;
//     return;
//   }

//   void clusterizer::set_GG_CELL_pitch (double p){
//     GG_CELL_pitch = p;
//     set_rad (GG_CELL_pitch / cos(M_PI/8.));
//     set_GG_CELL_diam (rad);
//     set_CellDistance (rad);
//     return;
//   }

//   void clusterizer::set_CellDistance(double v){
//     CellDistance = v;
//     return;
//   }

//   void clusterizer::set_SuperNemo(bool v){
//     SuperNemo = v;
//     return;
//   }

//   void clusterizer::set_SuperNemoChannel(bool v){
//     if (v)
//       {
//         set_SuperNemo (true);
//         SuperNemoChannel = true;
//         set_NemoraOutput (false);
//         set_N3_MC (false);
//         setDoDriftWires (false);
//         set_MaxBlockSize (1);
//       }
//     else
//       {
//         SuperNemoChannel = false;
//       }
//     return;
//   }

//   void clusterizer::set_NemoraOutput(bool no){
//     NemoraOutput = no;
//     return;
//   }

//   void clusterizer::set_N3_MC(bool v){
//     N3_MC = v;
//     return;
//   }

//   void clusterizer::set_FoilRadius(double v){
//     FoilRadius = v;
//     return;
//   }

//   void clusterizer::set_xsize(double v){
//     xsize = v;
//     return;
//   }

//   void clusterizer::set_ysize(double v){
//     ysize = v;
//     return;
//   }

//   void clusterizer::set_zsize(double v){
//     zsize = v;
//     return;
//   }

//   void clusterizer::set_bfield(double v){
//     bfield = v;
//     return;
//   }

}
