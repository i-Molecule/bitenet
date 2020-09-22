#include "cAtom.hpp"

void Atom::init() {
	serial_id_ = 0;
	residue_id_ = 0;
	alternative_location_ = '-';
	name_ = "";
	element_ = "-";
	residue_name_ = "";
	position_ = { 0, 0, 0 };
	occupancy_ = 0;
	temperature_factor_ = 0;
	chain_id_ = '-';
	record_name_ = "";
    type_ = "-";
	//charge_ = 0;
	type_ = "";
	//sasa_ = 0;
	//energy_electrostatic_ = 0;
	//energy_vdw_ = 0;
	//descriptor_evolution_ = {};
    return;
}

// Constructors

// Constructs Atom class with zeros variables
Atom::Atom() {
	init();
    return;
}

// element, position
Atom::Atom(const std::string& el, const Position<float>& p) {
	init();
	set_element(el);
	set_position(p);
    return;
}

// element, name, position
Atom::Atom(const std::string& el, const std::string& n, const Position<float>& p) {
	init();
	set_element(el);
	set_name(n);
	set_position(p);
    return;
}

// element, name, serial, position
Atom::Atom(const std::string& el, const std::string& n, const int i, const Position<float>& p) {
	init();
	set_element(el);
	set_name(n);
	set_serial_id(i);
	set_position(p);
    return;
}

Atom::~Atom() {}

// Setting class fields

void Atom::set_element(const std::string& el) {
	element_ = el;
    return;
}
const std::string& Atom::get_element() const{
	return element_;
}

void Atom::set_name(const std::string& n) {
	name_ = n;
    return;
}
const std::string& Atom::get_name() const{
	return name_;
}

void Atom::set_serial_id(const int i) {
	serial_id_ = i;
    return;
}
const int& Atom::get_serial_id() const {
	return serial_id_;
}

void Atom::set_residue_name(const std::string& r) {
	residue_name_ = r;
    return;
}
const std::string& Atom::get_residue_name() const {
	return residue_name_;
}

void Atom::set_residue_id(const int i) {
	residue_id_ = i;
    return;
}
const int& Atom::get_residue_id() const {
	return residue_id_;
}


void Atom::set_position(const Position<float>& p) {
	position_ = p;
    return;
}
const Position<float>& Atom::get_position() const {
	return position_;
}

void Atom::set_occupancy(const float& o) {
	occupancy_ = o;
}
const float& Atom::get_occupancy() const {
	return occupancy_;
}

void Atom::set_temperature_factor(const float& t) {
	temperature_factor_ = t;
    return;
}
const float& Atom::get_temperature_factor() const {
	return temperature_factor_;
}

void Atom::set_chain_id(const char c) {
	chain_id_ = c;
    return;
}
const char& Atom::get_chain_id() const {
	return chain_id_;
}

void Atom::set_record_name(const std::string& s) {
	record_name_ = s;
    return;
}
const std::string& Atom::get_record_name() const {
	return record_name_;
}

void Atom::set_alternative_location(const char a) {
	alternative_location_ = a;
    return;
}
const char& Atom::get_alternative_location() const {
	return alternative_location_;
}

void Atom::set_charge(const std::string& c) {
	charge_ = c;
    return;
}
const std::string& Atom::get_charge() const {
	return charge_;
}
/*
void Atom::set_sasa(double s) {
	sasa_ = s;
}
const double& Atom::sasa() const {
	return sasa_;
}

void Atom::set_energy_electrostatic(float e) {
	energy_electrostatic_ = e;
}
const float& Atom::energy_electrostatic() const {
	return energy_electrostatic_;
}

void Atom::set_energy_vdw(float e) {
	energy_vdw_ = e;
}
const float& Atom::energy_vdw() const {
	return energy_vdw_;
}
*/
void Atom::set_type(const std::string& t) {
	type_ = t;
    return;
}
const std::string& Atom::get_type() const {
	return type_;
}

float Atom::get_radius_vdw() const{
    return radii_vdw_.at(element_);
}


void Atom::assign_typization(const std::string& typization) {
    std::string key;
    if (typization=="ITSCORE") {
        
        if (name_=="N") key = "backbone N";
        else if (name_=="O") key = "backbone O";
        else if (name_=="C") key = "backbone C";
        else if (name_=="CA") key = "backbone CA";
        else if (name_=="OXT") key = "C-terminal OXT";
        else key = residue_name_+":"+name_;
        
        if (atom_types_.at(typization).find(key)!=atom_types_.at(typization).end()) {
            type_ = (atom_types_.at(typization)).at(key);
        }
        else {
            if (false)
            std::cerr << "warning: unknown key: " << key << " for typization: " << typization << std::endl;
            type_ = "-";
        }
    }
    else {
        if (false)
        std::cerr << "warning: unknown typization: " << typization << std::endl;
        type_ = "-";
    }
    
    return;
}

double Atom::calculate_atomic_density(const Position<float>& p, const float& squared_cutoff) const {
    double d = (p-position_).get_squared_length();
    if (d<squared_cutoff) return exponent_approximation(-d/2.0);
    return 0.0;
}

double Atom::calculate_atomic_density(const float& x, const float& y, const float& z, const float& squared_cutoff) const {
    
    Position<float> p;
    p.x_= x; p.y_ = y; p.z_ = z;
    return calculate_atomic_density(p, squared_cutoff);
}



double Atom::exponent_approximation(const double& y) const {
    double x = y;
    x = 1.0 + x / 256.0;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    return x;
}

/*
void Atom::set_descriptor_evolution(std::vector<float> evolDescriptors) {
	descriptor_evolution_ = evolDescriptors;
}

std::vector<float> Atom::descriptor_evolution() const {
	return descriptor_evolution_;
}

std::vector<Atom*> atomIndexer(std::vector<Atom> &atomVector) {
	// Returns vector of pointers to Atoms from vector of Atom
	std::vector<Atom*> atomIndexer;
	for (int i = 0; i < atomVector.size(); i++) {
		atomIndexer.push_back(&atomVector[i]);
	}
	return atomIndexer;
}
*/

// this is the same for all instances
const std::map< std::string, std::map<std::string, std::string > > Atom::atom_types_= Atom::set_atom_types();

const std::map<std::string, std::map<std::string, std::string> > Atom::set_atom_types() {
    std::map<std::string, std::map<std::string, std::string> > m;
    m["ITSCORE"] = {
        {"CYS:SG", "SULFUR"},
        {"MET:SD", "SULFUR"},
    
        {"ASN:ND2", "AMIDE"},
        {"GLN:NE2", "AMIDE"},
        {"backbone N", "AMIDE"},
    
        {"HIS:ND1", "AROMATIC"},
        {"HIS:NE1", "AROMATIC"},
        {"HIS:NE2", "AROMATIC"},
        {"TRP:NE1", "AROMATIC"},
    
        {"ARG:NE", "GUANIDINIUM"},
        {"ARG:NH1", "GUANIDINIUM"},
        {"ARG:NH2", "GUANIDINIUM"},
        {"ARG:NH3", "GUANIDINIUM"},
    
        {"LYS:NZ",  "AMMONIUM"},
    
        {"ASN:OD1", "CARBONYL"},
        {"GLN:OE1", "CARBONYL"},
        {"backbone O", "CARBONYL"},
    
        {"SER:OG", "HYDROXYL"},
        {"THR:OG1", "HYDROXYL"},
        {"TYR:OH", "HYDROXYL"},
    
        {"ASP:OD1", "CARBOXYL"},
        {"ASP:OD2", "CARBOXYL"},
        {"ASP:OD3", "CARBOXYL"},
        {"GLU:OE1", "CARBOXYL"},
        {"GLU:OE2", "CARBOXYL"},
        {"GLU:OE3", "CARBOXYL"},
        {"C-terminal O", "CARBOXYL"},
        {"C-terminal OXT", "CARBOXYL"},
    
        {"ARG:CZ", "CARBON_SP2"},
        {"ASN:CG", "CARBON_SP2"},
        {"ASP:CG", "CARBON_SP2"},
        {"GLN:CD", "CARBON_SP2"},
        {"GLU:CD", "CARBON_SP2"},
        {"backbone C", "CARBON_SP2"},
    
        {"HIS:CG", "CARBON_AROMATIC"},
        {"HIS:CD2", "CARBON_AROMATIC"},
        {"HIS:CE1", "CARBON_AROMATIC"},
        {"PHE:CG", "CARBON_AROMATIC"},
        {"PHE:CD1", "CARBON_AROMATIC"},
        {"PHE:CD2", "CARBON_AROMATIC"},
        {"PHE:CD3", "CARBON_AROMATIC"},
        {"PHE:CE1", "CARBON_AROMATIC"},
        {"PHE:CE2", "CARBON_AROMATIC"},
        {"PHE:CE3", "CARBON_AROMATIC"},
        {"PHE:CZ", "CARBON_AROMATIC"},
        {"TRP:CG", "CARBON_AROMATIC"},
        {"TRP:CD1", "CARBON_AROMATIC"},
        {"TRP:CD2", "CARBON_AROMATIC"},
        {"TRP:CD3", "CARBON_AROMATIC"},
        {"TRP:CE1", "CARBON_AROMATIC"},
        {"TRP:CE2", "CARBON_AROMATIC"},
        {"TRP:CE3", "CARBON_AROMATIC"},
        {"TRP:CZ1", "CARBON_AROMATIC"},
        {"TRP:CZ2", "CARBON_AROMATIC"},
        {"TRP:CZ3", "CARBON_AROMATIC"},
        {"TRP:CH2", "CARBON_AROMATIC"},
        {"TYR:CG", "CARBON_AROMATIC"},
        {"TYR:CD1", "CARBON_AROMATIC"},
        {"TYR:CD2", "CARBON_AROMATIC"},
        {"TYR:CD3", "CARBON_AROMATIC"},
        {"TYR:CE1", "CARBON_AROMATIC"},
        {"TYR:CE2", "CARBON_AROMATIC"},
        {"TYR:CE3", "CARBON_AROMATIC"},
        {"TYR:CZ", "CARBON_AROMATIC"},
    
        {"ALA:CB", "CARBON_SP3"},
        {"ARG:CB", "CARBON_SP3"},
        {"ARG:CD", "CARBON_SP3"},
        {"ARG:CG", "CARBON_SP3"},
        {"ASN:CB", "CARBON_SP3"},
        {"ASP:CB", "CARBON_SP3"},
        {"CYS:CB", "CARBON_SP3"},
        {"GLN:CB", "CARBON_SP3"},
        {"GLN:CG", "CARBON_SP3"},
        {"GLU:CB", "CARBON_SP3"},
        {"GLU:CG", "CARBON_SP3"},
        {"HIS:CB", "CARBON_SP3"},
        {"ILE:CB", "CARBON_SP3"},
        {"ILE:CG1", "CARBON_SP3"},
        {"ILE:CG2", "CARBON_SP3"},
        {"ILE:CG3", "CARBON_SP3"},
        {"ILE:CD1", "CARBON_SP3"},
        {"ILE:CD", "CARBON_SP3"},// added as trajectory pdb files has replacings ILE:CD1 -> ILE:CD        
        {"LEU:CB", "CARBON_SP3"},
        {"LEU:CG", "CARBON_SP3"},
        {"LEU:CD1", "CARBON_SP3"},
        {"LEU:CD2", "CARBON_SP3"},
        {"LEU:CD3", "CARBON_SP3"},
        {"LYS:CB", "CARBON_SP3"},
        {"LYS:CG", "CARBON_SP3"},
        {"LYS:CD", "CARBON_SP3"},
        {"LYS:CE", "CARBON_SP3"},
        {"MET:CB", "CARBON_SP3"},
        {"MET:CG", "CARBON_SP3"},
        {"MET:CE", "CARBON_SP3"},
        {"PHE:CB", "CARBON_SP3"},
        {"PRO:CB", "CARBON_SP3"},
        {"PRO:CG", "CARBON_SP3"},
        {"PRO:CD", "CARBON_SP3"},
        {"SER:CB", "CARBON_SP3"},
        {"THR:CB", "CARBON_SP3"},
        {"THR:CG2", "CARBON_SP3"},
        {"TRP:CB", "CARBON_SP3"},
        {"TYR:CB", "CARBON_SP3"},
        {"VAL:CB", "CARBON_SP3"},
        {"VAL:CG1", "CARBON_SP3"},
        {"VAL:CG2", "CARBON_SP3"},
        {"VAL:CG3", "CARBON_SP3"},
        {"backbone CA", "CARBON_SP3"},
    };
    
    m["-"] = {};
    
    return m;
}

const std::map< std::string, std::vector<std::string> > Atom::channel_names_ = Atom::set_channel_names();

const std::map< std::string, std::vector<std::string> > Atom::set_channel_names() {
    
    std::map< std::string, std::vector<std::string> > m;
    m["ITSCORE"] = {"SULFUR", "AMIDE", "AROMATIC", "GUANIDINIUM", "AMMONIUM", "CARBONYL", "HYDROXYL", "CARBOXYL", "CARBON_SP2", "CARBON_AROMATIC", "CARBON_SP3"};
    m["-"] = {};

    return m;
}

const std::vector<std::string> Atom::get_channel_names(const std::string& typization) {
    
    std::vector<std::string> s;
    
    if (typization=="ITSCORE") {
        s = Atom::channel_names_.at("ITSCORE");
    }
    else {
        std::cerr << "warning: unknown typization: " << typization << std::endl;
        s = Atom::channel_names_.at("-");
    }
    return s;
}

const std::map<std::string, float> Atom::radii_vdw_ = {
    { "H", 1.20 },
    { "C", 1.70 },
    { "N", 1.55 },
    { "O", 1.52 },
    { "S", 1.80 }
};

std::ostream& operator <<(std::ostream& os, const Atom& a) {
    os << "SERIAL ID: " << a.serial_id_ <<
    ", NAME: " << a.name_ <<
    ", ELEMENT: " << a.element_ <<
    ", TYPE: " << a.type_ <<
    ", RESIDUE ID: " << a.residue_id_ <<
    ", RESIDUE NAME: " << a.residue_name_ <<
    ", CHAIN ID: " << a.chain_id_ <<
    ", POSITION: " << a.position_;
    return os;
}




/*
constexpr std::vector<std::vector<float> > precompute_lj_grid() {

    std::vector<std::vector<float> > lj_grid;

  
    const float step = 0.2;
    const float r_min = 0.0;
    const float r_max = 10.0+step;

    const int ind_min = int(r_min);
    const int ind_max = int((r_max-r_min)/step);

    const std::vector<const std::string> atoms = {"H", "C", "N", "O", "S"};


    for (int ind = ind_min; ind < ind_max; ++ind) {
        for (const std::string& s1 : atoms) {
            for (const std::string& s2 : atoms) {
                lj_grid[s1][s2][ind] = 0.0; // formula goes here
            }
        }
    }
  
    return lj_grid;
}

const std::vector<std::vector<float> > Atom::lj_grid_ = precompute_lj_grid();
*/