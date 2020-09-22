/**********************************/
/*
This class aims to create 3d-cube of voxels given the size along the x-, y-, z- dimensions and number of voxels, respecitvely. 
*/
/**********************************/


#pragma once

#include <map>
#include <string>
#include <tgmath.h>
#include <tuple>
#include <future>

#include "cVoxel.hpp"


template <typename T>
class VoxelGrid {

private:
    
    // TODO: why do we need it here? it should be static member of voxel!
    bool init(const std::vector<std::string>* channel_names=NULL, const size_t threads=1); // names of the channels in a voxel

    // initial parameters of the grid
    //float size_;                              // Grid size
    float size_x_, size_y_, size_z_;   // length of each edge of the cube
    float voxel_size_;                  // length of a voxel edge (same in each direction)
    float spacing_;                     // distance between the voxel (same in each direction)
    Position<float> center_;            // coordinates of the center of the grid

    // other parameters
    float cell_size_;
    //unsigned int n_voxels_;             // Number of voxels along x, y, z axes
    unsigned int n_voxels_x_, n_voxels_y_, n_voxels_z_;             // Number of voxels along x, y, z axes
    float x_min_, x_max_, y_min_, y_max_, z_min_, z_max_; // grid borders

    void fit_grid_size(unsigned int& n_voxels, float& size, const char direction); // fit size to integer number of voxels

public:
    
    Voxel<T>*** voxels_ = NULL;

    bool silent_ = true;
    
    VoxelGrid();
    VoxelGrid(  const float& size, 
                const float& voxel_size, 
                const float& spacing,
                const Position<float>& center,
                const std::vector<std::string>* channel_names=NULL
             );

    VoxelGrid(  const float& size_x,
                const float& size_y,
                const float& size_z, 
                const float& voxel_size, 
                const float& spacing, 
                const Position<float>& center,
                const std::vector<std::string>* channel_names=NULL
             );
    ~VoxelGrid();
    
    Position<int> get_voxel_indexes(const Position<float>& p) const;  // check if coordinates belongs to voxel
    Position<float> get_voxel_center(const Position<int>& p) const;    // get x-, y-, z- coordinates of a voxel
    
    std::vector<Voxel<T>* > get_neighbours(const Position<int>& pos, int n_order) const; // get neigbours of n-order from the current voxel(i,j,k), including this voxel
    std::vector<Voxel<T>* > get_neighbours(const Position<float>& pos, int n_order) const; // get neigbours of n-order from the current point, including voxel containing this point
    std::vector<Voxel<T>* > get_neighbours(const Position<float>& p, const double& r_cutoff) const; // get neigbours within r_cutoff (converted to n_order) from the current point, including voxel containing this point
    
    std::vector<Voxel<T>* > get_neighbours(const Position<int>& pos, int n_order, 
        const Position<size_t>& start_indexes, const Position<size_t>& end_indexes) const;
    std::vector<Voxel<T>* > get_neighbours(const Position<float>& pos, const double& r_cutoff, 
        const Position<size_t>& start_indexes, const Position<size_t>& end_indexes) const;

    std::vector<Position<int> > get_neighbours_index(const Position<int>& pos, int n_order) const; // get neigbours of n-order from the current voxel(i,j,k), including this voxel
    std::vector<Position<int> > get_neighbours_index(const Position<float>& pos, int n_order) const; // get neigbours of n-order from the current point, including voxel containing this point
    std::vector<Position<int> > get_neighbours_index(const Position<float>& p, const double& r_cutoff) const; // get neigbours within r_cutoff (converted to n_order) from the current point, including voxel containing this point

    // GETTERS
    float get_voxel_size() const {return voxel_size_;};
    float get_spacing() const {return spacing_;};
    Position<float> get_sizes() const {return Position<float>{size_x_, size_y_, size_z_}; };
    Position<unsigned int> get_n_voxels() const {return Position<unsigned int>{n_voxels_x_, n_voxels_y_, n_voxels_z_}; };
    Position<float> get_floor() const {return Position<float>{x_min_, y_min_, z_min_};};

    bool init(  const float& size_x,
                const float& size_y,
                const float& size_z, 
                const float& voxel_size, 
                const float& spacing, 
                const Position<float>& center,
                const std::vector<std::string>* channel_names=NULL,
                const size_t threads=1);

};

template <typename T>
VoxelGrid<T>::VoxelGrid() 
    :   size_x_(0.0), size_y_(0.0), size_z_(0.0),
        spacing_(0.0),
        voxel_size_(1.0),
        center_({0.0, 0.0, 0.0})
{
    //init();
}

template <typename T>
VoxelGrid<T>::VoxelGrid(const float& size, const float& voxel_size, const float& spacing, const Position<float>& center, const std::vector<std::string>* channel_names) 
    :   size_x_(size), size_y_(size), size_z_(size),
        spacing_(spacing),
        voxel_size_(voxel_size), center_(center)
{
    init(channel_names);
}

template <typename T>
VoxelGrid<T>::VoxelGrid(const float& size_x, const float& size_y, const float& size_z, const float& voxel_size, const float& spacing, const Position<float>& center, const std::vector<std::string>* channel_names) 
    :   size_x_(size_x), size_y_(size_y), size_z_(size_z),
        spacing_(spacing), voxel_size_(voxel_size), center_(center)
{
    init(channel_names);
}

template <typename T>
VoxelGrid<T>::~VoxelGrid() {
    
    if (voxels_) {

        for (size_t i=0; i<n_voxels_x_; ++i) {

            if (voxels_[i]) {
                for (size_t j=0; j<n_voxels_y_; ++j) {
                    if (voxels_[i][j]) {
                        delete [] voxels_[i][j];
                    }
                }
                delete [] voxels_[i];
            }
        }
    } 
    delete [] voxels_;
}


template <typename T>
void VoxelGrid<T>::fit_grid_size(unsigned int& n_voxels, float& size, const char direction) {
    n_voxels = std::floor(size / cell_size_); // voxels+spacing within the size

    if ( float(n_voxels) < size / cell_size_) {
        if (!silent_)
        std::cerr   << "warning: resizing bounding box to fit integer number of voxels in " 
                    << direction 
                    << " direction" 
                    << std::endl;
        size = n_voxels * cell_size_ + voxel_size_; // add the last voxel
        n_voxels += 1;
    }

    return;
}

template <typename T>
bool VoxelGrid<T>::init(const std::vector<std::string>* channel_names, const size_t threads) {

    // check sizes
    if (size_x_<=0.0 || size_y_<=0.0 || size_z_<=0.0) {
        std::cerr<<"warning: cannot initialize voxel grid of non-positive size : " << size_x_ << " " << size_y_ << " " << size_z_ << std::endl;
        return false;
    }

    // check voxel size
    if (voxel_size_<=0.0) {
        std::cerr<<"warning: cannot initialize voxel grid of non-positive voxel size : "<< voxel_size_ << std::endl;
        return false;
    } 
    
    // check spacing
    if (spacing_<0.0) {
        std::cerr<<"warning: cannot initialize voxel grid of negative spacing: " << spacing_ << std::endl;
        return false;
    }

    // check pointer to voxels
    if (voxels_) {
        std::cerr<<"warning: voxel grid was initialized before - clearning memory"<<std::endl;
        delete [] voxels_;
    }

    voxels_ = NULL; // no content yet // TODO: need to protect this pointer from change

    cell_size_ = voxel_size_ + spacing_;

    // fit borders 
    fit_grid_size(n_voxels_x_, size_x_, 'x');
    fit_grid_size(n_voxels_y_, size_y_, 'y');
    fit_grid_size(n_voxels_z_, size_z_, 'z');
    
    // compute borders
    x_min_ = center_.x_ - 0.5*size_x_;
    x_max_ = center_.x_ + 0.5*size_x_;

    y_min_ = center_.y_ - 0.5*size_y_;
    y_max_ = center_.y_ + 0.5*size_y_;

    z_min_ = center_.z_ - 0.5*size_z_;
    z_max_ = center_.z_ + 0.5*size_z_;


    // local vars
    float margin_x = voxel_size_/2.0 + x_min_;
    float margin_y = voxel_size_/2.0 + y_min_;
    float margin_z = voxel_size_/2.0 + z_min_;

    voxels_ = new Voxel<T>** [n_voxels_x_];

    if (threads > 1){
        std::vector<std::future<void>> futures;
        size_t chunk_size = std::ceil(static_cast<float>(n_voxels_x_) / threads);
        for (size_t thread=0; thread < threads; ++thread){
            futures.push_back(
                std::async([this, thread, chunk_size, margin_x, margin_y, margin_z]{
                    Voxel<T>* pv;
                    for (size_t i=thread*chunk_size; i<n_voxels_x_ && i<(thread+1)*chunk_size; ++i) {
                        
                        voxels_[i] = new Voxel<T>* [n_voxels_y_];
                        for (size_t j=0; j<n_voxels_y_; ++j) {

                            voxels_[i][j] = new Voxel<T> [n_voxels_z_];
                            for (size_t k=0; k<n_voxels_z_; ++k) {
                                // initializing voxel coordinates
                                pv = &(voxels_[i][j][k]);
                                pv->p_.x_ = i*cell_size_ + margin_x;
                                pv->p_.y_ = j*cell_size_ + margin_y;
                                pv->p_.z_ = k*cell_size_ + margin_z;
                            }
                        }
                    }
                }
                )
            );
        }
    } else {
        Voxel<T>* pv;
        for (size_t i=0; i<n_voxels_x_; ++i) {
            
            voxels_[i] = new Voxel<T>* [n_voxels_y_];
            for (size_t j=0; j<n_voxels_y_; ++j) {

                voxels_[i][j] = new Voxel<T> [n_voxels_z_];
                for (size_t k=0; k<n_voxels_z_; ++k) {

                    // initializing voxel coordinates
                    pv = &(voxels_[i][j][k]);
                    pv->p_.x_ = i*cell_size_ + margin_x;
                    pv->p_.y_ = j*cell_size_ + margin_y;
                    pv->p_.z_ = k*cell_size_ + margin_z;
                    
                    // initializing voxel channels
                    // if (channel_names) pv->init_channels(*channel_names);
                }
            }
        }
    }

    unsigned int volume = n_voxels_x_ * n_voxels_y_ * n_voxels_z_;  
    if (!silent_)  
    printf("Creating voxel grid of %.3f*%.3f*%.3f units, containing %d*%d*%d=%d voxels of size %.3f with %.3f spacing\n", size_x_, size_y_, size_z_, n_voxels_x_, n_voxels_y_, n_voxels_z_, volume, voxel_size_, spacing_);
    //printf("X : %8.3f - %8.3f\n", x_min_, x_max_);
    //printf("Y : %8.3f - %8.3f\n", y_min_, y_max_);
    //printf("Z : %8.3f - %8.3f\n", z_min_, z_max_);

    return true;
}


template <typename T>
bool VoxelGrid<T>::init(  
                const float& size_x,
                const float& size_y,
                const float& size_z, 
                const float& voxel_size, 
                const float& spacing, 
                const Position<float>& center,
                const std::vector<std::string>* channel_names,
                const size_t threads) {

    size_x_ = size_x;
    size_y_ = size_y;
    size_z_ = size_z;
    voxel_size_ = voxel_size;
    spacing_ = spacing;
    center_ = center;

    return init(channel_names, threads);
}

template <typename T>
Position<int> VoxelGrid<T>::get_voxel_indexes(const Position<float>& p) const {
    
    // given 3d coordinate return determines the corresponding voxels and returns its position

    // std::cout << p << "  " << Position<float>{x_min_, y_min_, z_min_} << std::endl;

    int ind_i = std::floor((p.x_ - x_min_)/cell_size_);
    int ind_j = std::floor((p.y_ - y_min_)/cell_size_);
    int ind_k = std::floor((p.z_ - z_min_)/cell_size_);

    // std::cout << ind_i << " " << ind_j << " " << ind_k << std::endl;

    if (ind_i<0 || ind_j<0 || ind_k<0 || 
        ind_i>=n_voxels_x_ || ind_j>=n_voxels_y_ || ind_k>=n_voxels_z_) {
        if (!silent_)
        printf("Point (%8.3f,%8.3f,%8.3f) is not in grid \n", p.x_, p.y_, p.z_);
        ind_i = -1;
        ind_j = -1;
        ind_k = -1;
        return Position<int> {ind_i, ind_j, ind_k};
    }
    
    if ( ((p.x_ - x_min_) -ind_i*cell_size_ > voxel_size_) || ((p.y_ - y_min_) - ind_j*cell_size_ > voxel_size_) || 
        ((p.z_ - z_min_) - ind_k*cell_size_ > voxel_size_) ) {
        if (!silent_)
        printf("Point (%8.3f,%8.3f,%8.3f) is not in voxel \n", p.x_, p.y_, p.z_);
        ind_i = -1;
        ind_j = -1;
        ind_k = -1;
        return Position<int> {ind_i, ind_j, ind_k};
    }
    
    return Position<int> {ind_i, ind_j, ind_k};
}


template <typename T>
std::vector<Voxel<T>* > VoxelGrid<T>::get_neighbours(const Position<int>& pos, int n_order) const
{
    
    // stores pointers to n-layer neibourhood of given voxel, including itself
    std::vector<Voxel<T>* > vs;
    int i = pos.x_;
    int j = pos.y_;
    int k = pos.z_;
    
    for (int a=-n_order; a<=n_order; ++a) {
        if (i+a<0 || i+a>=n_voxels_x_) continue;

        for (int b=-n_order; b<=n_order; ++b) {
            if (j+b<0 || j+b>=n_voxels_y_) continue;

            for (int c=-n_order; c<=n_order; ++c) {
                if (k+c<0 || k+c>=n_voxels_z_) continue;

                vs.push_back(&voxels_[i+a][j+b][k+c]);
            }
        }
    }
    return vs;
}

template <typename T>
std::vector<Position<int> > VoxelGrid<T>::get_neighbours_index(const Position<int>& pos, int n_order) const
{
    
    // stores pointers to n-layer neibourhood of given voxel, including itself
    std::vector<Position<int> > vs;
    int i = pos.x_;
    int j = pos.y_;
    int k = pos.z_;
    
    for (int a=-n_order; a<=n_order; ++a) {
        if (i+a<0 || i+a>=n_voxels_x_) continue;

        for (int b=-n_order; b<=n_order; ++b) {
            if (j+b<0 || j+b>=n_voxels_y_) continue;

            for (int c=-n_order; c<=n_order; ++c) {
                if (k+c<0 || k+c>=n_voxels_z_) continue;

                vs.push_back(Position<int>{i+a,j+b,k+c});
            }
        }
    }
    return vs;
}

template <typename T>
std::vector<Voxel<T>* > VoxelGrid<T>::get_neighbours(const Position<float>& p, int n_order) const{

    // get voxel that contains the point
    Position<int> pos = get_voxel_indexes(p);
    return get_neighbours(pos, n_order);
}

template <typename T>
std::vector<Position<int> > VoxelGrid<T>::get_neighbours_index(const Position<float>& p, int n_order) const{

    // get voxel that contains the point
    Position<int> pos = get_voxel_indexes(p);
    return get_neighbours_index(pos, n_order);
}

template <typename T>
std::vector<Voxel<T>* > VoxelGrid<T>::get_neighbours(const Position<float>& p, const double& r_cutoff) const{

    // get voxel that contains the point
    Position<int> pos = get_voxel_indexes(p);

    // get n_order that covers r_cutoff
    int n_order = std::ceil(r_cutoff/cell_size_);

    return get_neighbours(pos, n_order);
}



template <typename T>
std::vector<Position<int> > VoxelGrid<T>::get_neighbours_index(const Position<float>& p, const double& r_cutoff) const{

    // get voxel that contains the point
    Position<int> pos = get_voxel_indexes(p);

    // get n_order that covers r_cutoff
    int n_order = std::ceil(r_cutoff/cell_size_);

    return get_neighbours_index(pos, n_order);
}


template <typename T>
std::vector<Voxel<T>* > VoxelGrid<T>::get_neighbours(const Position<int>& pos, int n_order, 
        const Position<size_t>& start_indexes, const Position<size_t>& end_indexes) const{

    std::vector<Voxel<T>* > vs;
    int i = pos.x_;
    int j = pos.y_;
    int k = pos.z_;
    
    for (int a=-n_order; a<=n_order; ++a) {
        if (i+a<0 || i+a>=n_voxels_x_ || i+a<start_indexes[0] || i+a>=end_indexes[0]) continue;

        for (int b=-n_order; b<=n_order; ++b) {
            if (j+b<0 || j+b>=n_voxels_y_ || j+b<start_indexes[1] || j+b>=end_indexes[1]) continue;

            for (int c=-n_order; c<=n_order; ++c) {
                if (k+c<0 || k+c>=n_voxels_z_ || k+c<start_indexes[2] || k+c>=end_indexes[2]) continue;

                vs.push_back(&voxels_[i+a][j+b][k+c]);
            }
        }
    }
    return vs;
}
template <typename T>
std::vector<Voxel<T>* > VoxelGrid<T>::get_neighbours(const Position<float>& p, const double& r_cutoff, 
        const Position<size_t>& start_indexes, const Position<size_t>& end_indexes) const{

    Position<int> pos = get_voxel_indexes(p);
    int n_order = std::ceil(r_cutoff/cell_size_);
    return get_neighbours(pos, n_order, start_indexes, end_indexes);
    
}


