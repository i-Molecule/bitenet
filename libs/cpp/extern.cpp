#include "cDescriptor.hpp"

extern "C"{

    PDBfile* pdb_new(const char* file_name, 
                // const char* chain_id,
                const char* typization="ITSCORE",
                const unsigned int atoms_in_ligand_cutoff=14,
                const unsigned int atoms_in_interface_cutoff = 20,
                const float interface_radius = 4.0,
                const unsigned int atoms_in_protein_cutoff = 400,
                const bool filter_ligands = true);
                // const bool flag_filter_binding_sites = true,
                // const bool only_heavy_atoms = true,
                // const std::set<std::string> *atom_types = NULL);

    void rotate(PDBfile* pdb, const double theta, const double phi, const double psi);
    void center(PDBfile* pdb);
    void compute_boxes(PDBfile* pdb, const float buffer=0.);
    void write(PDBfile* pdb, const char* filename, 
        const bool flag_write_protein=true, const bool flag_write_ligands=true);
    int get_interface_number(PDBfile* pdb);
    float* export_interfaces(PDBfile* pdb);

    void pdb_del(PDBfile *pdb);

    Descriptor* descriptor_new(PDBfile* pdb);
                
    void init_grid(Descriptor* descriptor,
        const float voxel_size = 1.,
        const float spacing = 0.0,
        const float r_cutoff=4.,
        const int threads=1);
        // const std::vector<std::string> *channel_names = NULL,
        // bool flag_init_protein = true,
        // bool flag_init_ligand = true,
        //std::set<std::string>* atom_types=NULL);
    void calculate_density(Descriptor*, const float r_cutoff=4., const int threads=1);
    int* get_grid_size(Descriptor*);
    float* get_grid_floor(Descriptor*);
    // int get_interface_number(Descriptor*);
    double* export_grid(Descriptor*, const int threads=1);
    // bool* export_grid_bool(Descriptor*);
    // float* export_interfaces(Descriptor*);
    void descriptor_del(Descriptor*);

}

extern "C"{

    PDBfile* pdb_new(const char* file_name, 
        // const char* chain_id,
        const char* typization,
        const unsigned int atoms_in_ligand_cutoff,
        const unsigned int atoms_in_interface_cutoff,
        const float interface_radius,
        const unsigned int atoms_in_protein_cutoff,
        const bool filter_ligands){

        const bool flag_filter_binding_sites = true;
        const bool only_heavy_atoms = true;
        const std::set<std::string> *atom_types = NULL;

        const bool silent = true;

        const std::string file_name_ = file_name;
        // const char chain_id_ = chain_id[0];
        const char chain_id_ = 'A';

        // const std::string typization_ = typization;

        return new PDBfile(
            file_name_,
            chain_id_,
            typization,
            atoms_in_ligand_cutoff,
            atoms_in_interface_cutoff,
            interface_radius,
            atoms_in_protein_cutoff,
            flag_filter_binding_sites,
            only_heavy_atoms,
            atom_types,
            silent,
            filter_ligands);
    }

    void rotate(PDBfile* pdb, const double theta, const double phi, const double psi){
        pdb->rotate(theta, phi, psi);
        return;
    }
    void center(PDBfile* pdb){
        pdb->center();
        return;
    }
    void compute_boxes(PDBfile* pdb, const float buffer){
        pdb->compute_bounding_box(buffer);
        pdb->compute_interface_boxes();
        return;
    }
    void write(PDBfile* pdb, const char* filename, const bool flag_write_protein, const bool flag_write_ligands){
        pdb->writePDBfile(filename, flag_write_protein, flag_write_ligands);
        return;
    }
    int get_interface_number(PDBfile* pdb){
        return int(pdb->interfaces_.size());
    }
    float* export_interfaces(PDBfile* pdb){
        return pdb->export_interfaces();
    }
    void pdb_del(PDBfile* pdb){
        delete pdb;
        return;
    }

    Descriptor* descriptor_new(PDBfile* pdb){
        return new Descriptor(*pdb);
    
    }

    void init_grid(Descriptor* descriptor,
        const float voxel_size,
        const float spacing,
        const float r_cutoff,
        const int threads){

        bool flag_init_protein = false;
        bool flag_init_ligand = false;
        std::set<std::string>* atom_types=NULL;
        
        std::vector<std::string> channel_names = Atom::get_channel_names(descriptor->pdb_->typization_);
        descriptor->init_grid(voxel_size, spacing, &channel_names, 
            flag_init_protein, flag_init_ligand);
        descriptor->calculate_grid_atomic_density(r_cutoff, NULL, threads);
        return;
    }

    void calculate_density(Descriptor* descriptor, const float r_cutoff, const int threads){
        descriptor->calculate_grid_atomic_density(r_cutoff, NULL, threads);
        return;
    }

    int* get_grid_size(Descriptor* descriptor){
        auto sizes = descriptor->grid_.get_n_voxels();
        unsigned int n_channels = descriptor->channel_names_.size();
        static int res[4];
        res[0] = sizes.x_;
        res[1] = sizes.y_;
        res[2] = sizes.z_;
        res[3] = n_channels;
        return res;
    }
    float* get_grid_floor(Descriptor* descriptor){
        auto grid_floor = descriptor->grid_.get_floor();
        static float res[3];
        res[0] = grid_floor.x_;
        res[1] = grid_floor.y_;
        res[2] = grid_floor.z_;
        return res;
    }
    // int get_interface_number(Descriptor* descriptor){
    //     return int(descriptor->pdb_->interfaces_.size());
    // }

    double* export_grid(Descriptor* descriptor, const int threads){
        return descriptor->export_grid(threads);
    }
    // bool* export_grid_bool(Descriptor*  descriptor){
    //     return descriptor->export_grid_bool();
    // }
    // float* export_interfaces(Descriptor* descriptor){
    //     return descriptor->export_interfaces();
    // }
    void descriptor_del(Descriptor* descriptor){
        delete descriptor;
        return;
    }
}

