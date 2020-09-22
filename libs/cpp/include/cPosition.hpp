#pragma once
#include <algorithm>
#include <tgmath.h>
#include <iostream>


// Classes and functions for better and easier dealing with dimensional(Angstroms) quantities (atoms' and voxels positions)

// Class for position, three dimensional orthogonal coordinates {x, y, z}


// T is supposed to be number type, e.g. float/double/int
template <typename T>
class Position {

public:

	T x_, y_, z_;								// {X, Y, Z}

	Position();
	Position(const T& x, const T& y, const T& z);

	void set_value(const T& x, const T& y, const T& z);	// Sets position values x, y, z
    void set_value(const Position<T>& pos);              // Sets position values to pos.x, pos.y, pos.z

	T get_length() const;					// Returns position vector length ( Sqrt(x^2 + y^2 + z^2) )
    T get_squared_length() const;             // Returns position vector length x^2 + y^2 + z^2 )

	Position<T> operator+(const Position<T>& pos) const;
	Position<T> operator-(const Position<T>& pos) const;
	Position<T> operator*(const T& t) const;
	Position<T> operator/(const T& t) const;
	Position<T>& operator+=(const Position<T>& pos);
    Position<T>& operator*=(const T& t);
	Position<T>& operator/=(const T& t);
    Position<T>& operator=(const Position<T>& p);
    //friend std::ostream& operator<<(std::ostream& os, const Position<T>&);
    T operator[](size_t i) const;
	T & operator[](size_t i);
};


// Matrix 3x3 class used for calculation rotation matrix from inertia matrix for protein rotation
template <typename T>
class Matrix3 {
	
public:
	
	Position<T> rows_[3];
	
	Matrix3();
	Matrix3(const Position<T> rows[3]);
	Matrix3(const Position<T> &pos1, const Position<T> &pos2, const Position<T> &pos3);
	
	T det() const;
	Matrix3<T> transposed() const;
	Matrix3<T> inv() const;
	Position<T> eigenvalues() const;
	Matrix3<T> eigenvectors() const;
	
	Matrix3<T> operator*(const Matrix3<T>& matrix) const;
	Position<T> operator*(const Position<T>& pos) const;
    Matrix3<T>& operator*=(const T& t);
	Matrix3<T>& operator/=(const T& t);
	Position<T> operator[](size_t i) const;
	Position<T>& operator[](size_t i);
};

template <typename T>
Position<T>::Position() {
    set_value(0, 0, 0);
}

template <typename T>
Position<T>::Position(const T& val1, const T& val2, const T& val3) {
    set_value(val1, val2, val3);
}

template <typename T>
void Position<T>::set_value(const Position<T>& pos) {
    x_ = pos.x_;
    y_ = pos.y_;
    z_ = pos.z_;
}

template <typename T>
void Position<T>::set_value(const T& val1, const T& val2, const T& val3) {
    x_ = val1;
    y_ = val2;
    z_ = val3;
}

template <typename T>
T Position<T>::get_length() const {
    return sqrt(get_squared_length());
}

template <typename T>
T Position<T>::get_squared_length() const {
    return x_ * x_ + y_ * y_ + z_ * z_;
}

template <typename T>
Position<T> Position<T>::operator+(const Position<T>& pos2) const {
    
    Position<T> pos;
    pos.x_ = this->x_ + pos2.x_;
    pos.y_ = this->y_ + pos2.y_;
    pos.z_ = this->z_ + pos2.z_;
    return pos;
}

template <typename T>
Position<T> Position<T>::operator-(const Position<T>& pos2) const {
    
    Position<T> pos;
    pos.x_ = this->x_ - pos2.x_;
    pos.y_ = this->y_ - pos2.y_;
    pos.z_ = this->z_ - pos2.z_;
    return pos;
}


template <typename T>
Position<T> Position<T>::operator*(const T& value) const {
    
    Position<T> pos;
    pos.x_ = this->x_ * value;
    pos.y_ = this->y_ * value;
    pos.z_ = this->z_ * value;
    return pos;
}

template <typename T>
Position<T> Position<T>::operator/(const T& value) const {
    
    Position<T> pos;
    pos.x_ = this->x_ / value;
    pos.y_ = this->y_ / value;
    pos.z_ = this->z_ / value;
    return pos;
}

template <typename T>
Position<T>& Position<T>::operator+=(const Position<T>& pos2){
    this->x_ += pos2.x_;
    this->y_ += pos2.y_;
    this->z_ += pos2.z_;
    return *this;
}

template <typename T>
Position<T>& Position<T>::operator*=(const T& value){
    this->x_ *= value;
    this->y_ *= value;
    this->z_ *= value;
    return *this;
}

template <typename T>
Position<T>& Position<T>::operator/=(const T& value){
    this->x_ /= value;
    this->y_ /= value;
    this->z_ /= value;
    return *this;
}

template <typename T>
Position<T>& Position<T>::operator=(const Position<T>& p){
    this->x_ = p.x_;
    this->y_ = p.y_;
    this->z_ = p.z_;
    return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Position<T>& pos) {
    os << pos.x_ << " " << pos.y_ << " " << pos.z_;
    return os;
}

template <typename T>
bool operator==(const Position<T>& p1, const Position<T>& p2) {
    return (p1.x_==p2.x_ && p1.y_==p2.y_ && p1.z_==p2.z_ );
}

template <typename T>
T Position<T>::operator[](size_t i) const{
	if (i == 0)
		return this->x_;
	else if (i == 1)
		return this->y_;
	else if (i == 2) 
		return this->z_;
	else
		throw std::runtime_error("Position operator[]: index " + std::to_string(i) + " is not valid");
}

template <typename T>
T& Position<T>::operator[] (size_t i){
	if (i == 0)
		return this->x_;
	else if (i == 1)
		return this->y_;
	else if (i == 2) 
		return this->z_;
	else
		throw std::runtime_error("Position operator[]: index " + std::to_string(i) + " is not valid");
}

template <typename T1, typename T2>
Position<T1> operator*(const Position<T1>& pos, const T2& v){
    Position<T1> p (0, 0, 0);
    p.x_ = pos.x_ * v;
    p.y_ = pos.y_ * v;
    p.z_ = pos.z_ * v;
    return p;
}

// class for hash function 
template <typename T>
class Position_Hash_Function { 
public: 
    int operator()(const Position<T>& p) const
    { 
        return 100*p.x_ + 10*p.y_ + p.z_; 
    } 
}; 

template <typename T>
T dot_scalar(const Position<T>& pos1, const Position<T>& pos2){
    return pos1.x_*pos2.x_ + pos1.y_*pos2.y_ + pos1.z_*pos2.z_;
}

template <typename T>
Position<T> dot_vector(const Position<T>& pos1, const Position<T>& pos2){
    Position<T> pos;
    pos.x_ = pos1.y_*pos2.z_ - pos1.z_*pos2.y_;
    pos.y_ = pos1.z_*pos2.x_ - pos1.x_*pos2.z_;
    pos.z_ = pos1.x_*pos2.y_ - pos1.y_*pos2.x_;
    return pos;
}

template <typename T>
Matrix3<T>::Matrix3(){
	rows_[0] = Position<T>(0, 0, 0);
	rows_[1] = Position<T>(0, 0, 0);
	rows_[2] = Position<T>(0, 0, 0);
}

template <typename T>
Matrix3<T>::Matrix3(const Position<T> rows[3]){
	rows_[0] = rows[0];
	rows_[1] = rows[1];
	rows_[2] = rows[2];
}

template <typename T>
Matrix3<T>::Matrix3(const Position<T> &pos1, const Position<T> &pos2, const Position<T> &pos3){
	rows_[0] = pos1;
	rows_[1] = pos2;
	rows_[2] = pos3;
}

template <typename T>
T Matrix3<T>::det() const{
	return ((*this)[0][0]*((*this)[1][1]*(*this)[2][2] - (*this)[1][2]*(*this)[2][1]) - 
				(*this)[0][1]*((*this)[1][0]*(*this)[2][2] - (*this)[2][0]*(*this)[1][2]) + 
				(*this)[0][2]*((*this)[1][0]*(*this)[2][1] - (*this)[2][0]*(*this)[1][1]));
}

template <typename T>
Matrix3<T> Matrix3<T>::transposed() const{
	Matrix3<T> matrix({Position<T>((*this)[0][0], (*this)[1][0], (*this)[2][0]),
		Position<T>((*this)[0][1], (*this)[1][1], (*this)[2][1]),
		Position<T>((*this)[0][2], (*this)[1][2], (*this)[2][2])});
	return matrix;
}

template <typename T>
Matrix3<T> Matrix3<T>::inv() const{
	Matrix3<T> matrix({Position<T>((*this)[1][1]*(*this)[2][2] - (*this)[2][1]*(*this)[1][2], 
								(*this)[0][2]*(*this)[2][1] - (*this)[2][2]*(*this)[0][1],
								(*this)[0][1]*(*this)[1][2] - (*this)[1][1]*(*this)[0][2]),
					Position<T>((*this)[1][2]*(*this)[2][0] - (*this)[2][2]*(*this)[1][0],
							(*this)[0][0]*(*this)[2][2] - (*this)[2][0]*(*this)[0][2],
							(*this)[0][2]*(*this)[1][0] - (*this)[1][2]*(*this)[0][0]),
					Position<T>((*this)[1][0]*(*this)[2][1] - (*this)[2][0]*(*this)[1][1],
							(*this)[0][1]*(*this)[2][0] - (*this)[2][1]*(*this)[0][0],
							(*this)[0][0]*(*this)[1][1] - (*this)[1][0]*(*this)[0][1])});
	matrix /= this->det();
	return matrix;
}

template <typename T>
Position<T> Matrix3<T>::eigenvalues() const{
    // calculating eigenvalues for real symmetrical matrix 3x3
    // from https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
	T p1 = (*this)[0][1]*(*this)[0][1] + (*this)[0][2]*(*this)[0][2] + (*this)[1][2]*(*this)[1][2];
	Position<T> eig(0., 0., 0.);
	if (p1 == 0){
		eig.x_ = (*this)[0][0];
		eig.y_ = (*this)[1][1];
		eig.z_ = (*this)[2][2];
	} else {
		T q = ((*this)[0][0] + (*this)[1][1] + (*this)[2][2]) / 3;
		T p2 = ((*this)[0][0] - q)*((*this)[0][0] - q) + ((*this)[1][1] - q)*((*this)[1][1] - q) + \
			((*this)[2][2] - q)*((*this)[2][2] - q) + 2 * p1;
		T p = sqrt(p2 / 6);
		Matrix3<T> B(Position<T>((1 / p) * ((*this)[0][0] - q), (1 / p) * (*this)[0][1], (1 / p) * (*this)[0][2]),
					Position<T>((1 / p) * (*this)[1][0], (1 / p) * ((*this)[1][1] - q), (1 / p) * (*this)[1][2]),
					Position<T>((1 / p) * (*this)[2][0], (1 / p) * (*this)[2][1], (1 / p) * ((*this)[2][2] - q)));
		T r = B.det() / 2;
		const double pi = std::atan(1)*4;
		double phi;
		if (r <= -1){
			phi = pi / 3;
		} else if (r >= 1) {
			phi = 0;
		} else {
			phi = std::acos(r) / 3;
		}
		
		//eig.x_ = q + 2 * p * std::cos(phi);
		//eig.z_ = q + 2 * p * std::cos(phi + (2*pi/3));
        eig.z_ = q + 2 * p * std::cos(phi);
		eig.x_ = q + 2 * p * std::cos(phi + (2*pi/3));
		eig.y_ = 3*q - eig.x_ - eig.z_;
        // eig.x_ <= eig.y_ <= eig.z_ 
	}
	return eig;
}

template <typename T>
Matrix3<T> Matrix3<T>::eigenvectors() const{
    // calculating eigenvectors 
	Matrix3<T> v(Position<T>(1., 1., 1.), Position<T>(1., 1., 1.), Position<T>(1., 1., 1.));
	Position<T> eig = this->eigenvalues();
	for (size_t i=0; i<3; ++i){
		v.rows_[i].y_ = - v[i][2] * ((*this)[1][2] - (*this)[0][2] * (*this)[1][0] / ((*this)[0][0] - eig[i] + 1e-20)) / 
			((*this)[1][1] - eig[i] - (*this)[0][1] * (*this)[1][0] / ((*this)[0][0] - eig[i] + 1e-20) + 1e-20);
		v.rows_[i].x_ = - 1. / ((*this)[0][0] - eig[i] + 1e-20) * ((*this)[0][1] * v[i][1] + (*this)[0][2] * v[i][2]);
		T norm = std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
		v.rows_[i] /= norm;
	}
    // check if it is right or left handed basis, inverse the last vector
    if (dot_scalar(v[0], dot_vector(v[1], v[2])) < 0){
        v[2] *= -1.;
    }
	v = v.transposed();
	return v;
}

template <typename T>
Matrix3<T> Matrix3<T>::operator*(const Matrix3<T>& matrix) const{
	Matrix3<T> C({Position<T>(0., 0., 0.), Position<T>(0., 0., 0.), Position<T>(0., 0., 0.)});
	for (size_t i=0; i<3; ++i){
		for (size_t j=0; j<3; ++j){
			C[i][j] = (*this)[i][0]*matrix[0][j] + (*this)[i][1]*matrix[1][j] + (*this)[i][2]*matrix[2][j];
		}
	}
	return C;
}

template <typename T>
Position<T> Matrix3<T>::operator*(const Position<T>& pos) const{
	Position<T> w(0., 0., 0.);
	for (size_t i=0; i<3; ++i){
		w[i] = (*this)[i][0]*pos[0] + (*this)[i][1]*pos[1] + (*this)[i][2]*pos[2];
	}
	return w;
}

template <typename T>
Matrix3<T>& Matrix3<T>::operator/=(const T& t){
	this->rows_[0] /= t;
	this->rows_[1] /= t;
	this->rows_[2] /= t;
	return *this;
}

template <typename T>
Matrix3<T>& Matrix3<T>::operator*=(const T& t){
	this->rows_[0] *= t;
	this->rows_[1] *= t;
	this->rows_[2] *= t;
	return *this;
}

template <typename T>
Position<T> Matrix3<T>::operator[](size_t i) const{
	if ((i >= 0) && (i <= 2))
		return this->rows_[i];
	else
		throw std::runtime_error("Matrix3 operator[]: index " + std::to_string(i) + " is not valid");
}

template <typename T>
Position<T>& Matrix3<T>::operator[] (size_t i){
	if ((i >= 0) && (i <= 2))
		return this->rows_[i];
	else
		throw std::runtime_error("Matrix3 operator[]: index " + std::to_string(i) + " is not valid");
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix3<T>& matrix) {
	os << matrix[0] << std::endl;
	os << matrix[1] << std::endl;
	os << matrix[2];
    return os;
}

template <typename T>
Matrix3<T> rotation_matrix(const T& theta, const T& phi, const T& psi){
    Position<T> u(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    T c = cos(psi), s = sin(psi), C = 1 - cos(psi);
    Matrix3<T> m;
    m[0][0] = u[0]*u[0]*C + c;
    m[0][1] = u[0]*u[1]*C - u[2] * s;
    m[0][2] = u[0]*u[2]*C + u[1] * s;
    m[1][0] = u[1]*u[0]*C + u[2] * s;
    m[1][1] = u[1]*u[1]*C + c;
    m[1][2] = u[1]*u[2]*C - u[0] * s;
    m[2][0] = u[2]*u[0]*C - u[1] * s;
    m[2][1] = u[2]*u[1]*C + u[0] * s;
    m[2][2] = u[2]*u[2]*C + c;
    return m;
}

template <typename T1, typename T2>
Position<T2> operator*(const Matrix3<T1>& matrix, const Position<T2>& pos){
    Position<T2> w(0., 0., 0.);
	for (size_t i=0; i<3; ++i){
		w[i] = matrix[i][0]*pos[0] + matrix[i][1]*pos[1] + matrix[i][2]*pos[2];
	}
	return w;
}