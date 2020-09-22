#include "cDescriptor.hpp"
#include <future>

double exponent_approximation(const double& y);


void Descriptor::init_grid( const float& voxel_size, 
                            const float& spacing, 
                            const std::vector<std::string>* channel_names, 
                            bool flag_init_protein, 
                            bool flag_init_ligand,
                            const size_t threads) {

    channel_names_ = *channel_names;
    
    if (pdb_==NULL) {
        std::cerr<<"cannot initialize grid - no pdb file"<<std::endl;
        return;
    }

    // if (!flag_init_protein && !flag_init_ligand) {
    //     std::cerr << "pass at least one flag_init_protein/ligand as true" << std::endl;
    //     return;
    // }


    // set up voxel grid parameters
    std::vector<float> borders = pdb_->get_borders();

    float size_x = borders[3] - borders[0];
    float size_y = borders[4] - borders[1];
    float size_z = borders[5] - borders[2];

    float center_x = 0.5*(borders[3] + borders[0]);
    float center_y = 0.5*(borders[4] + borders[1]);
    float center_z = 0.5*(borders[5] + borders[2]);

    Position<float> center(center_x, center_y, center_z);
    grid_.silent_ = pdb_->silent_;

    // initializing grid
    flag_grid_initialized_ = grid_.init(size_x, size_y, size_z, voxel_size, spacing, center, channel_names, threads);

    //std::cout << size_x << " " << size_y << " " << size_z << " " <<std::endl;
    //std::cout << center << std::endl;
    //std::cout   << borders[0] << " : " << borders[3] << "\n" 
    //            << borders[1] << " : " << borders[4] << "\n" 
    //            << borders[2] << " : " << borders[5] << "\n" 
    //            << std::endl;
    
    if (!flag_grid_initialized_) return;

    // local vars    
    Position<float> pos;    // local position
    Atom* a;                // local atom
    std::string t;          // local type
    float x,y,z;            // local coordinates
    Voxel<Atom>* v;         // local pointer to voxel
    Position<int> index;     // local voxel indexes

    if (flag_init_protein)
    {
        for (size_t p = 0; p < pdb_->protein_.size(); ++p)
        {
            // fill voxel
            a = &pdb_->protein_[p];
            pos = a->get_position();
            index = grid_.get_voxel_indexes(pos);

            if (index.x_ == -1 || index.y_ == -1 || index.z_ == -1) {
                continue;
            }

            //std::cout<<index<<std::endl;
            v = &grid_.voxels_[index.x_][index.y_][index.z_];
            v->add_object(a);
        }
    }

    if (flag_init_ligand)
    {
        for (auto &[ligand, vector] : pdb_->ligands_map_)
        {
            for (size_t p = 0; p < vector.size(); ++p)
            {
                // fill voxel
                a = &vector[p];
                pos = a->get_position();
                index = grid_.get_voxel_indexes(pos);
                if (index.x_ == -1 || index.y_ == -1 || index.z_ == -1) {
                    continue;
                }
                v = &grid_.voxels_[index.x_][index.y_][index.z_];
                v->add_object(a);
            }
        }
    }

    return;

}


void Descriptor::calculate_grid_atomic_density(float r_cutoff, std::set<std::string>* atom_types, size_t threads) {
    
    if (threads > 1){
        std::vector<std::future<void>> futures;
        size_t chunk_size = std::ceil(static_cast<float>(grid_.get_n_voxels()[0]) / threads);
        for (size_t thread=0; thread < threads; ++thread){
            futures.push_back(
                std::async([this, thread, chunk_size, r_cutoff]{
                    std::vector<Voxel<Atom>* > neighbors;
                    for (const Atom& a : pdb_->protein_){
                        neighbors = grid_.get_neighbours(a.get_position(), r_cutoff, 
                            {thread*chunk_size, 0, 0}, 
                            {(thread+1)*chunk_size, grid_.get_n_voxels()[1], grid_.get_n_voxels()[2]});
                        for (auto v: neighbors){
                            v->channels_[a.get_type()] += calculate_atomic_density(a, v->p_, r_cutoff);
                        }
                    }
                }
                )
            );
        }
    } else {
        // local vars    
        Position<float> pos;    // local position
        const Atom* a;                // local atom
        std::string t;          // local type
        float x,y,z;            // local coordinates
        Voxel<Atom>* v;         // local pointer to voxel
        Position<int> index;     // local voxel indexes
        // TODO: Maybe store pointers to relevant atoms upon initialization and run over them?
        std::vector<Voxel<Atom>* > neighbors;
        for (size_t p=0; p<pdb_->protein_.size(); ++p) {
            
            // fill voxel
            a = &pdb_->protein_[p];
            t = a->get_type();

            pos = a->get_position();
            // calculate density in neibouring voxels
            neighbors = grid_.get_neighbours(pos, r_cutoff);
            for (auto v: neighbors) {
                v->channels_[t] += calculate_atomic_density((*a), v->p_, r_cutoff);
            }
        }
    }
}


double Descriptor::calculate_atomic_density(const Atom& a, const Position<float>& pos, double r_cutoff) {
    double d = (a.get_position() - pos).get_squared_length();
    if (d<r_cutoff*r_cutoff) {
        return std::exp(-0.5*d);
        //return exponent_approximation(-0.5*d);
    }
    return 0.0;
}


void Descriptor::write_grid_csv(const std::string& file_name) {


    if (grid_.voxels_ == NULL) {
        std::cerr   << "cannot write empty grid (grid_.voxels_=NULL)" 
                    << std::endl;
        return;
    }

    std::ofstream os;
    os.open(file_name);

    // writing header
    os << "i,j,k";
    const auto channels = grid_.voxels_[0][0][0].channels_;
    for (const auto& [channel, value] : channels) {
        os << "," << channel;
    }
    os << '\n';

    Position<unsigned int> n_voxels = grid_.get_n_voxels();

    Voxel<Atom>* v;
    for (size_t i=0; i<n_voxels.x_; ++i) {
        for (size_t j=0; j<n_voxels.y_; ++j) {
            for (size_t k=0; k<n_voxels.z_; ++k) {
                v = &grid_.voxels_[i][j][k];
                os  << i << "," << j << "," << k;
                for (const auto& [channel, value] : v->channels_) {
                    os << "," << value;
                }
                os << '\n';
            }
        }
    }
    os.close();
}

double* Descriptor::export_grid(const size_t threads){
    auto n_voxels = grid_.get_n_voxels();
    if (export_res_){
        delete [] export_res_;
        export_res_ = NULL;
    }
    export_res_ = new double[n_voxels[0] * n_voxels[1] * n_voxels[2] * int(channel_names_.size())];
    if (threads > 1){
        std::vector<std::future<void>> futures;
        size_t chunk_size = std::ceil(static_cast<float>(n_voxels.x_) / threads);
        for (size_t thread=0; thread < threads; ++thread){
            futures.push_back(
                std::async([this, thread, chunk_size, n_voxels]{
                    unsigned int index;
                    std::map<std::string, double>* channels;
                    for (size_t i=thread*chunk_size; i<n_voxels.x_ && i<(thread+1)*chunk_size; ++i) {
                        for (size_t j=0; j<n_voxels.y_; ++j) {
                            for (size_t k=0; k<n_voxels.z_; ++k) {
                                channels = &grid_.voxels_[i][j][k].channels_;
                                for (size_t l=0; l<channel_names_.size(); ++l){
                                    index = i*n_voxels[1]*n_voxels[2]*int(channel_names_.size()) + j*n_voxels[2]*int(channel_names_.size()) +
                                        k*int(channel_names_.size()) + l;
                                    if (channels->find(channel_names_[l]) != channels->end()) 
                                        export_res_[index] = (*channels)[channel_names_[l]];
                                    else
                                        export_res_[index] = 0.;
                                }
                            }
                        }
                    }
                }
            ));
        }
    } else {
        unsigned int index;
        std::map<std::string, double>* channels;
        for (size_t i=0; i<n_voxels.x_; ++i) {
            for (size_t j=0; j<n_voxels.y_; ++j) {
                for (size_t k=0; k<n_voxels.z_; ++k) {
                    channels = &grid_.voxels_[i][j][k].channels_;
                    for (size_t l=0; l<channel_names_.size(); ++l){
                        index = i*n_voxels[1]*n_voxels[2]*int(channel_names_.size()) + j*n_voxels[2]*int(channel_names_.size()) +
                            k*int(channel_names_.size()) + l;
                        if (channels->find(channel_names_[l]) != channels->end()) 
                            export_res_[index] = (*channels)[channel_names_[l]];
                        else
                            export_res_[index] = 0.;
                    }
                }
            }
        }
    }
    return export_res_;
}

void Descriptor::determine_interface(double r_cutoff) {

    std::unordered_set< Position<int>,  Position_Hash_Function<int> > neighbours;
    std::vector<Position<int> > vs;
    Voxel<Atom> *voxel;

    for (const auto& [ligand, atoms] : pdb_->ligands_map_) {
        
        auto key = ligand;

        for (size_t i=0; i<atoms.size(); ++i) {
            vs = grid_.get_neighbours_index(atoms[i].get_position(), r_cutoff);
            for (const auto& v : vs) {
                neighbours.insert(std::move(v));
            }
        }

        int counter = 0;

        for (const auto& n : neighbours) {
            voxel = &grid_.voxels_[n.x_][n.y_][n.z_];
            for (Atom *obj : voxel->objects_) {
                pdb_->interfaces_[key].push_back(obj);
                counter+=1;
            }
        }
    }

    std::cout   << "There are " 
                << pdb_->interfaces_.size() 
                << " ligand-binding interfaces for " 
                << pdb_->get_file_name() 
                << " chain " 
                << pdb_->get_chain_id() 
                << std::endl;

    for (const auto& interface: pdb_->interfaces_) {
        printf("\t-- %lu atoms on interface with ligand %s\n", (interface.second).size(), interface.first.name_.c_str() );
    }

    return;
}

double exponent_approximation(const double& y) {
    double x = y;
    x = 1.0 + x / 256.0;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    return x;
}