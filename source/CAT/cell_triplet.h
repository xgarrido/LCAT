/* -*- mode: c++ -*- */
#ifndef FALAISE_CAT_CELL_TRIPLET_H
#define FALAISE_CAT_CELL_TRIPLET_H 1

#include <string>
#include <iostream>
#include <vector>
#include <CAT/cell.h>
#include <CAT/cell_couplet.h>
#include <CAT/joint.h>


namespace CAT{

  /// \brief A cell_triplet is composed of three cells and a list of joints
  class cell_triplet
  {


  protected:

    // first cell
    cell ca_;

    // second cell
    cell cb_;

    // third cell
    cell cc_;

    // list of chi2 values
    std::vector<double> chi2s_;

    // list of prob values
    std::vector<double> probs_;

  public:

    // list of joints
    std::vector<joint> joints_;

    // status of cell triplet
    bool free_;

    // begun cell triplet
    bool begun_;

    //!Default constructor
    cell_triplet();

    //!Default destructor
    virtual ~cell_triplet();

    //! constructor
    cell_triplet(const cell &ca, const cell &cb, const cell &cc);

    /*** dump ***/
    virtual void dump (std::ostream & a_out         = std::clog,
                       const std::string & a_title  = "",
                       const std::string & a_indent = "",
                       bool a_inherit          = false) const;

    // /*** dump ***/
    // virtual void dump_joint (joint j,
    //   		       std::ostream & a_out         = std::clog,
    //   		       const std::string & a_title  = "",
    //   		       const std::string & a_indent = "",
    //   		       bool a_inherit          = false) const;


    //! set cells
    void set(const cell_couplet &cca, const cell_couplet &ccb);

    //! set cells
    void set(const cell &ca, const cell &cb, const cell &cc);

    //! set free level
    void set_free(bool free);
    //! set begun level

    void set_begun(bool begun);

    //! set joints
    void set_joints(const std::vector<joint> & joints);

    //! set chi2 list
    void set_chi2s(const std::vector<double> & chi2s);

    //! set prob list
    void set_probs(const std::vector<double> & probs);

    //! get first cell couplet
    cell_couplet cca();

    //! get second cell couplet
    cell_couplet ccb();

    //! get joints
    const std::vector<joint>& joints() const;

    //! get first cell
    const cell& ca()const;

    //! get second cell
    const cell& cb()const;

    //! get third cell
    const cell& cc()const;

    //! get list of chi2
    const std::vector<double>& chi2s() const;

    //! get list of prob
    const std::vector<double>& probs() const;

    //! get free level
    bool free()const;

    //! get begun level
    bool begun()const;

  public:

    void calculate_joints(double ratio_, double phi_limit_);

    std::vector<joint> refine(const std::vector<joint> & joints, double Ratio, size_t max_njoints=4);

    size_t iteration()const;

    cell_triplet invert();

    void set_all_used();

    friend bool operator==(const cell_triplet& left,
                           const cell_triplet& right);

    bool same_last_cell(cell c)const;

  };

}

#endif // FALAISE_CAT_CELL_TRIPLET_H
