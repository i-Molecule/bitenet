#pragma once

#include <vector>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
//#include <iomanip>
//#include <thread>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <utility>
//#include <tgmath.h>

#include "cAtom.hpp"

// Defining tuples for ligands and residues names for easier work with it
//#define LigandName std::tuple<std::string, char, int>
//#define ResidueName std::tuple<std::string, char, int>
struct Ligand {
    std::string name_="-";
    char chain_='-';
    int id_=0;
    std::vector<std::vector<int> > conects_; // TODO: add relevant conections from pdb file
};
bool operator==(const Ligand& l1, const Ligand& l2);
bool operator<(const Ligand& l1, const Ligand& l2);
std::ostream& operator<<(std::ostream& os, const Ligand& l);

struct Residue {
    //int id_=0;
    //std::string name_="RES";
    std::vector<Atom*> atoms_;
};

template <typename T>
struct Box {
    Position<T> center;
    Position<T> size;
};

class PDBfile {

private:
    // Dictionary for short amino acid names
    static const std::map<std::string, char> amino_acids_3_to_1_;
    // List of aminoacids
    static const std::set<std::string> amino_acids_;

    // Dictionary for replacing of residue names, e.g. HSD -> HIS
    static const std::map<std::string, std::string> amino_replace_;
    
    // set of names that are not considered as ligands
    static const std::set<std::string> unknown_ligands_;
    static const std::set<std::string> detergents_;
    static const std::set<std::string> list_of_modified_residues_;
    static const std::set<std::string> list_of_cofactors_;
    static const std::set<std::string> list_of_carbohydrates_;
    
    std::set<std::string> not_a_ligand_;    
    // ligand names that are discarded (union of unknown_ligands_ and detergents_)
    // not a static const because it varies depending on CLR // TODO? FIX IT?
    
    std::set<Ligand> modified_residues_;

    std::string file_name_; // name of the pdb file
    char chain_id_;         // protein chain to consider

    const int max_number_of_atoms = 100000; // to allocate memory for atoms

    unsigned int atoms_in_ligand_cutoff_;    // Compound is considered to be ligands if number of its heavy atoms is more than cutoff and its name is not in list below. (for binding sites)
    unsigned int atoms_in_interface_cutoff_; // Interface is considered as relevant, if the number of atoms in the interface >= this cutoff (for binding sites)
    unsigned int atoms_in_protein_cutoff_;   // Proteins with smaller number of atoms are disregarded (for binding sites)

    float interface_radius_;        // Two atoms interacts if the distance between them <= this cutoff (Angstroms) (for binding sites)
    float interface_squared_radius_;

    // borders of the protein
    float x_min_ = 0.0;
    float y_min_ = 0.0;
    float z_min_ = 0.0;
    float x_max_ = 0.0;
    float y_max_ = 0.0;
    float z_max_ = 0.0;

    void init(); // basic initializer in each constructor
    void init(const std::string &file_name,
              const char chain_id,
              const std::string &typization = "-",
              const unsigned int atoms_in_ligand_cutoff = 14,
              const unsigned int atoms_in_interface_cutoff = 20,
              const float &interface_radius = 4.0,
              const unsigned int atoms_in_protein_cutoff = 400,
              const bool flag_filter_binding_sites = true,
              const bool only_heavy_atoms = true,
              const std::set<std::string> *atom_types = NULL);

  public:
    
    std::vector<Atom> protein_;                         // vector of atoms of the protein
    std::map<Ligand, std::vector<Atom> > ligands_map_;  // atoms in a ligand
    std::map<Ligand, std::vector<Atom*> > interfaces_; // one interface per ligand, stores pointers to atoms in protein_

    std::map<Ligand, Box<float>> interface_boxes_;

    std::map<int, Residue> residues_;                    // vector of residues, [residue_id, Residue]
    
    std::string typization_ = "-";        // key word for the type assignemt
    
    bool silent_ = false;
    std::string add_header = "";
    float* export_res_ = NULL;
    bool filter_ligands_ = true;

    PDBfile();
    PDBfile(    const std::string& file_name, 
                const char chain_id, 
                const std::string& typization="-",
                const unsigned int atoms_in_ligand_cutoff=14,
                const unsigned int atoms_in_interface_cutoff = 20,
                const float& interface_radius = 4.0,
                const unsigned int atoms_in_protein_cutoff = 400,
                const bool flag_filter_binding_sites=true,
                const bool only_heavy_atoms=true,
                const std::set<std::string>* atom_types = NULL,
                const bool silent=false,
                const bool filter_ligands=true
            );
    ~PDBfile();
    
    void readPDBFile(const bool only_heavy_atoms=true, const std::set<std::string>* atom_types=NULL);
    Atom readPDBlineAtom(const std::string& pdb_line) const;
    std::string writePDBlineAtom(const Atom& a) const;
    void print_interface(const Ligand& l, const std::string& fn, const bool capped=false, const bool flat_split_protein_ligand=false); // print .pdb file with interface residues and corresponding ligand
    
    void writePDBfile(const std::string& filename, 
        const bool& flag_write_protein=true, const bool& flag_write_ligands=true) const;

    void filter_binding_sites();  // filter non-relevant ligands, determine interface residues
    void center();       // move COM to O(0,0,0) and rotate with respect to x,y,z (can look to modeller function)
    void rotate(const double& theta, const double& phi, const double& psi);
    void compute_bounding_box(const float& buffer=0.0);  // compute size of the box that fit protein with ligands

    void compute_interface_boxes();
    float* export_interfaces();

    // GETTERS
    std::string get_file_name()             const {return file_name_;};
    char  get_chain_id()                    const {return chain_id_;};
    float get_interface_radius()            const {return interface_radius_;};
    float get_interface_squared_radius_()   const {return interface_squared_radius_;};
    float get_atoms_in_ligand_cutoff()      const {return atoms_in_ligand_cutoff_;};
    float get_atoms_in_interface_cutoff()   const {return atoms_in_interface_cutoff_;};
    float get_atoms_in_protein_cutoff()     const {return atoms_in_protein_cutoff_;};
    std::vector<float> get_borders()        const {return std::vector<float> {x_min_, y_min_, z_min_, x_max_, y_max_, z_max_};};
};