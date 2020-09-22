#pragma once
// Voxel is a 3d cube with pointers to some object and feature channels
// Definitions are here to avoid declaration of all possible typenames in template

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>

#include "cPosition.hpp"

// Class for voxels for voxel grid
template <typename T>
class Voxel {
    
public:

    // MEMBERS
    std::vector<T*> objects_;       // Vector of pointers to objects

    std::map<std::string, double> channels_; // Vector of channels for descriptors, key of map corresponds to the code of the channel
    
    Position<float> p_; // voxel coordinates

    // CTORS
    Voxel() : p_({0.0, 0.0, 0.0}){};
    Voxel(const Position<float>& p) : p_(p) {};
    Voxel(const size_t& number_of_channels) : p_({0.0, 0.0, 0.0}) {objects_.reserve(number_of_channels);};
    Voxel(const Position<float>& p, const size_t& number_of_channels) : p_(p) {
        objects_.reserve(number_of_channels);
    };
    ~Voxel(){};

    
    void add_object(T* obj);        // add pointer to object to voxel
    void del_object(T* obj);        // del pointer to object from voxel
    
    void init_channels(const std::vector<std::string>& channel_names); // initialize channels with 0.0 values
    
    void clear();                   // clear voxel

    //friend std::ostream& operator<<(std::ostream& os, const Voxel<T>& vox);
};
 
template <typename T>
void Voxel<T>::clear() {
    
    objects_.clear();
    channels_.clear();
    p_=Position<float>{0.0, 0.0, 0.0};

    return;
}

template <typename T>
void Voxel<T>::add_object(T* obj) {
    
    if (obj!=NULL) {
        objects_.push_back(obj);
    }
    else {
        std::cerr<<"warning: passing NULL pointer"<<std::endl;
    }
    
    return;
}

template <typename T>
void Voxel<T>::del_object(T* obj) {
    if (obj!=NULL) {
        objects_.erase(std::remove(objects_.begin(), objects_.end(), obj), objects_.end()); // remove all elements with value obj
    }
    else {
        std::cerr<<"warning: passing NULL pointer"<<std::endl;
    }
    
    return;
}

template <typename T>
void Voxel<T>::init_channels(const std::vector<std::string>& channel_names) {
    for (auto s: channel_names) {
        channels_[s];
    }
    return;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Voxel<T>& v) {
    os << "VOXEL POSITION : " << v.p_ << 
    " VOXEL SIZE : " << v.objects_.size() <<
    " VOXEL CHANNELS :";
    for (auto [key, value] : v.channels_) {
        os << " (" << key << "," << value << ")";
    }
    return os;
}