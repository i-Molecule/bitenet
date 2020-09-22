#pragma once

#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cstdlib>


#include "cPDBfile.hpp"
#include "cVoxelGrid.hpp"

class Descriptor
{
  public:
    // MEMBERS

    PDBfile *pdb_;
    VoxelGrid<Atom> grid_;

    std::vector<std::string> channel_names_;


    double* export_res_ = NULL;
    // bool* export_res_bool_ = NULL;

    // CTORS

    Descriptor() = delete;
    Descriptor(PDBfile &pdb)
    {
        pdb_ = &pdb;
        //init_grid();
    };

    ~Descriptor(){
      if (export_res_) delete [] export_res_; 
      // if (export_res_bool_) delete [] export_res_bool_;
    };

    // FUNCTIONS

    void determine_interface(double r_cutoff=10.0);
    // SINGLE MOLECULE DESCRIPTORS

    // voxel grid
    void calculate_grid_atomic_density(float r_cutoff, std::set<std::string>* atom_types=NULL, size_t threads=1);

    double calculate_atomic_density(const Atom& a, const Position<float>& pos, double r_cutoff);

    void write_grid_csv(const std::string& file_name);

    double* export_grid(const size_t threads=1);
    // bool* export_grid_bool();
    // int* export_grid_bool_sparse();

    //calculate_grid_vdw();
    //calculate_grid_el();
    //calculate_grid_rdf();

    // vectors

    // INTERACTION DESCRIPTORS

    // voxel grid
    //calculate_grid_rdf();

    // vectors
    //calculate_sif();

    // THIRD-PARTY-DESCRIPTORS
    //calculate_splif();  // chemaxon
    //calculate_ecfp();   // chemaxon
    //calculate_charge(); // terachem
    //calculate_rdkit();  // rdkit
    //calculate_evolution(); // hhblits
    // ...

    void init_grid(const float &voxel_size,
                   const float &spacing,
                   const std::vector<std::string> *channel_names = NULL,
                   bool flag_init_protein = true,
                   bool flag_init_ligand = true,
                   const size_t threads=1);

  private:

    bool flag_grid_initialized_ = false;
};