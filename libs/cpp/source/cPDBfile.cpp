#include "cPDBfile.hpp"


bool operator==(const Ligand& l1, const Ligand& l2) {
    return (l1.name_==l2.name_ && l1.chain_==l2.chain_ && l1.id_==l2.id_);
}
bool operator<(const Ligand& l1, const Ligand& l2) {
    if (l1.name_!=l2.name_) return (l1.name_<l2.name_);
    else if (l1.chain_!=l2.chain_) return (l1.chain_<l2.chain_);
    else return (l1.id_<l2.id_);
    // return (l1.name_<l2.name_);
}
std::ostream& operator<<(std::ostream& os, const Ligand& l) {
    os << l.name_ << " " << l.chain_ << " " << l.id_;
    return os;
}

const std::map<std::string, char> PDBfile::amino_acids_3_to_1_ = {
    { "GLY", 'G' },{ "ALA", 'A' },{ "VAL", 'V' },{ "ILE", 'I' },{ "LEU", 'L' },
    { "PRO", 'P' },{ "SER", 'S' },{ "THR", 'T' },{ "CYS", 'C' },{ "MET", 'M' },
    { "ASP", 'D' },{ "ASN", 'N' },{ "GLU", 'E' },{ "GLN", 'Q' },{ "LYS", 'K' },
    { "ARG", 'R' },{ "HIS", 'H' },{ "PHE", 'F' },{ "TYR", 'Y' },{ "TRP", 'W' }
};
const std::set<std::string> PDBfile::amino_acids_ = {
    "GLY", "ALA", "VAL", "ILE", "LEU", "PRO", "SER", "THR", "CYS", "MET",
    "ASP", "ASN", "GLU", "GLN", "LYS", "ARG", "HIS", "PHE", "TYR", "TRP" };
const std::map<std::string, std::string> PDBfile::amino_replace_ = {
    {"HSD", "HIS"}, {"HSE", "HIS"}, {"HSP", "HIS"}
};
const std::set<std::string> PDBfile::unknown_ligands_ = { "UNL" };
const std::set<std::string> PDBfile::detergents_ = { 
    "12P", "15P", "1PE", "2CV", "2PE", "78M", "78N", "BCR", "BNG", "BOG",
    "BTB", "C14", "C8E", "CDL", "CLR", "CM5", "CPS", "DAO", "DGA", "DGD",
    "DMU", "DXC", "L2P", "LDA", "LFA", "LHG", "LI1", "LMG", "LMT", "LMU",
    "MPG", "MQ7", "MTN", "MYR", "MYS", "NS5", "OLA", "OLB", "OLC", "P4C",
    "P6G", "PC1", "PCW", "PE4", "PE5", "PE8", "PEE", "PG6", "PGV", "PL9",
    "PLC", "OLM", "PTY", "PX4", "SPN", "SPO", "SQD", "STE", "TGL", "U10", 
    "Y01", "HOH"
};
// from website : https://www.globalphasing.com/buster/manual/maketnt/manual/lib_val/library_validation.html
const std::set<std::string> PDBfile::list_of_modified_residues_ = {
    "CSD", "HYP", "BMT", "5HP", "ABA", "AIB", "CSW", "OCS", "DAL", "DAR",
    "DSG", "DSP", "DCY", "CRO", "DGL", "DGN", "DHI", "DIL", "DIV", "DLE",
    "DLY", "DPN", "DPR", "DSN", "DTH", "DTR", "DTY", "DVA", "CGU", "KCX",
    "LLP", "CXM", "FME", "MLE", "MVA", "NLE", "PTR", "ORN", "SEP", "TPO",
    "PCA", "PVL", "SAR", "CEA", "CSO", "CSS", "CSX", "CME", "TYS", "TPQ",
    "STY",
    "NH2", "CBX", "ACE", "FOR", "IVA", "BOC", // does not correspond to any residue
    "LYR", "MSE", // others
};
const std::set<std::string> PDBfile::list_of_cofactors_ = {
    "ADP", "AMP", "ATP", "CMP", "COA", "FAD", "FMN", "NAP", "NDP"
};
const std::set<std::string> PDBfile::list_of_carbohydrates_ = {
    "BGC", "GLC", "MAN", "BMA", "FUC", "GAL", "GLA", "NAG", "NGA", "SIA",
    "XYS"
};

PDBfile::PDBfile(){
    init();
}

PDBfile::PDBfile(   const std::string& file_name, 
                    const char chain_id, 
                    const std::string& typization,
                    const unsigned int atoms_in_ligand_cutoff,
                    const unsigned int atoms_in_interface_cutoff,
                    const float& interface_radius,
                    const unsigned int atoms_in_protein_cutoff,
                    const bool flag_filter_binding_sites,
                    const bool only_heavy_atoms,
                    const std::set<std::string>* atom_types,
                    const bool silent,
                    const bool filter_ligands
                    ) 
    :   file_name_(file_name), 
        chain_id_(chain_id),
        typization_ (typization),
        atoms_in_ligand_cutoff_ (atoms_in_ligand_cutoff),
        atoms_in_interface_cutoff_ (atoms_in_interface_cutoff),
        interface_radius_ (interface_radius), // Angstroms
        interface_squared_radius_ (interface_radius_*interface_radius_),
        atoms_in_protein_cutoff_ (atoms_in_protein_cutoff), // heavy atoms
        silent_ (silent),
        filter_ligands_ (filter_ligands)
    {
    init();
    readPDBFile(only_heavy_atoms, atom_types);
    if (flag_filter_binding_sites) {
        filter_binding_sites();
    }
    compute_bounding_box();
    //std::cout << typization_ << " " << typization << std::endl;
}

PDBfile::~PDBfile() {if (export_res_) delete [] export_res_;}

void PDBfile::init() {
    not_a_ligand_ = detergents_;
    not_a_ligand_.insert(unknown_ligands_.begin(), unknown_ligands_.end());
    not_a_ligand_.insert(list_of_modified_residues_.begin(), list_of_modified_residues_.end());
    not_a_ligand_.insert(list_of_cofactors_.begin(), list_of_cofactors_.end());
    not_a_ligand_.insert(list_of_carbohydrates_.begin(), list_of_carbohydrates_.end());


    protein_.clear();                         // vector of atoms of the protein
    ligands_map_.clear();  // atoms in a ligand
    interfaces_.clear();
    protein_.reserve(max_number_of_atoms);
}


void PDBfile::init( const std::string &file_name,
              const char chain_id,
              const std::string &typization,
              const unsigned int atoms_in_ligand_cutoff,
              const unsigned int atoms_in_interface_cutoff,
              const float &interface_radius,
              const unsigned int atoms_in_protein_cutoff,
              const bool flag_filter_binding_sites,
              const bool only_heavy_atoms,
              const std::set<std::string> *atom_types) 

    { 

        file_name_ = file_name; 
        chain_id_ = chain_id;
        typization_ = typization;
        atoms_in_ligand_cutoff_ = atoms_in_ligand_cutoff;
        atoms_in_interface_cutoff_ = atoms_in_interface_cutoff;
        interface_radius_ = interface_radius; // Angstroms
        interface_squared_radius_  = interface_radius_*interface_radius_;
        atoms_in_protein_cutoff_  = atoms_in_protein_cutoff; // heavy atoms

        init();
        readPDBFile(only_heavy_atoms, atom_types);
        if (flag_filter_binding_sites) {
            filter_binding_sites();
        }
        
}


void PDBfile::readPDBFile(const bool only_heavy_atoms, const std::set<std::string>* atom_types) {
    
    // Reading PDB file
    // Compound is considered to be a ligand if number of its heavy atoms is more than cutoff and its name is not in list below.
    // Returns list of atoms arrays, array[0] is protein, others are ligands
    
    std::vector<Atom> protein;
    char chain = '-';             // to track protein chains
    bool flag_protein_chain = false;    // to track protein chains
    bool flag_protein_ter = false;      // to track protein chains
    
    std::ifstream pdbfile(file_name_);
    
    if (pdbfile.is_open()) {
        
        if (!silent_)
        printf("Reading %s ...\n", file_name_.c_str());
        
        while (!pdbfile.eof()) {
            
            std::string line;
            getline(pdbfile, line);
            
            std::string record = line.substr(0, 6);
            
            //if cholesterol exist in header or title of pdb (it is complex with CLR) then it will be a ligand
            if (((record == "HEADER") || (record == "TITLE ")) && (line.find("CHOLESTEROL") != std::string::npos)) {
                not_a_ligand_.erase("CLR");
                add_header = "HEADER + CHOLESTEROL";
            }
            
            //check the list of modified residues not to be ligand
            /*
             Record Format
             
             COLUMNS        DATA TYPE     FIELD       DEFINITION
             --------------------------------------------------------------------------------
             1 -  6        Record name   "MODRES"
             8 - 11        IDcode        idCode      ID code of this entry.
             13 - 15        Residue name  resName     Residue name used in this entry.
             17             Character     chainID     Chain identifier.
             19 - 22        Integer       seqNum      Sequence number.
             23             AChar         iCode       Insertion code.
             25 - 27        Residue name  stdRes      Standard residue name.
             30 - 70        String        comment     Description of the residue modification.
             
            */
            if (record == "MODRES") {
                
                std::string resname = line.substr(12, 3);
                resname.erase(std::remove(resname.begin(), resname.end(), ' '), resname.end());
                char chain = line[16];
                int seqNum = std::stoi(line.substr(18, 4));
                Ligand l{resname, chain, seqNum};
                
                modified_residues_.insert(l);
                continue;
            }
            
            // check protein atoms
            // TODO: optimize copying local vars
            if (record == "ATOM  ") {
                Atom atom = readPDBlineAtom(line);
                // if (atom.get_chain_id() != chain_id_) continue;
                
                if ( amino_replace_.find(atom.get_residue_name()) != amino_replace_.end()){
                    atom.set_residue_name(amino_replace_.at(atom.get_residue_name()));
                    atom.assign_typization(typization_);
                }

                if ( std::find(amino_acids_.begin(), amino_acids_.end(), atom.get_residue_name() ) != amino_acids_.end() ) {

                    // check atom types
                    if ( atom_types!=NULL && atom_types->find(atom.get_element())==atom_types->end() ) continue;
                    
                    // check heavy atoms
                    if ( only_heavy_atoms && ((atom.get_element() == "H") || (atom.get_element() == "D"))) continue;
                        
                    if (flag_protein_ter && !silent_) {
                        printf("\t-- protein %c\n", chain_id_);
                        std::cerr<<"warning: there might be several protein chains with the same chain id" << std::endl;
                    }
                    protein_.push_back(atom);
                    residues_[protein_.back().get_residue_id()].atoms_.push_back(&protein_.back());
                }

                else {
                    if (!silent_)
                    std::cerr << "warning: unknown residue name: " << atom.get_residue_name() << std::endl;
                }

                continue;
            }
            
            // check ligand atoms
            if (record == "HETATM") {
                Atom atom = readPDBlineAtom(line);
                chain = atom.get_chain_id();

                    // check atom types
                    if ( atom_types!=NULL && atom_types->find(atom.get_element())==atom_types->end() ) continue;
                    
                    // check heavy atoms
                    if ( only_heavy_atoms && ((atom.get_element() == "H") || (atom.get_element() == "D"))) continue;
                    
                Ligand l{atom.get_residue_name(), atom.get_chain_id(), atom.get_residue_id()};

                if (ligands_map_.find(l) == ligands_map_.end()) {
                    if (!silent_)
                    printf("\t-- ligand %s %c %d\n", atom.get_residue_name().c_str(), atom.get_chain_id(), atom.get_residue_id());
                    ligands_map_.insert({ l,{ atom } });
                }
                else {
                    ligands_map_[l].push_back(atom);
                }
                
                continue;
            }
            if (record.substr(0, 3)=="TER"){
                flag_protein_ter = true;
                if (!silent_)
                if (chain==chain_id_) printf("\t-- protein %c\n", chain_id_);
            }
            
            if (record.substr(0, 3) == "END") break;
        }
    }
    else {
        std::cerr << "cannot open file " + file_name_ << std::endl;
    }
    pdbfile.close();
    return;
}

Atom PDBfile::readPDBlineAtom(const std::string& line) const {
    
    /*
     ATOM records in pdb file:
     
     COLUMNS        DATA  TYPE    FIELD        DEFINITION
     -------------------------------------------------------------------------------------
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
     13 - 16        Atom          name         Atom name.
     17             Character     altLoc       Alternate location indicator.
     18 - 20        Residue name  resName      Residue name.
     22             Character     chainID      Chain identifier.
     23 - 26        Integer       resSeq       Residue sequence number.
     27             AChar         iCode        Code for insertion of residues.
     31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
     39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
     47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
     55 - 60        Real(6.2)     occupancy    Occupancy.
     61 - 66        Real(6.2)     tempFactor   Temperature  factor.
     77 - 78        LString(2)    element      Element symbol, right-justified.
     79 - 80        LString(2)    charge       Charge  on the atom.
     
     In pdbqt file added:
     71 - 76                                      Partial charge
     78 - 79                                      AutoDock atom type                    */
    
    /* HETATM record type:
     
     COLUMNS       DATA  TYPE     FIELD         DEFINITION
     -----------------------------------------------------------------------
     1 - 6        Record name    "HETATM"
     7 - 11       Integer        serial        Atom serial number.
     13 - 16       Atom           name          Atom name.
     17            Character      altLoc        Alternate location indicator.
     18 - 20       Residue name   resName       Residue name.
     22            Character      chainID       Chain identifier.
     23 - 26       Integer        resSeq        Residue sequence number.
     27            AChar          iCode         Code for insertion of residues.
     31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
     39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
     47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
     55 - 60       Real(6.2)      occupancy     Occupancy.
     61 - 66       Real(6.2)      tempFactor    Temperature factor.
     77 - 78       LString(2)     element       Element symbol; right-justified.
     79 - 80       LString(2)     charge        Charge on the atom.
     */
    
    std::string record_name = line.substr(0, 6);

    int serial_id = stoi(line.substr(6, 5));

    std::string name = line.substr(12, 4);
    name.erase(std::remove(name.begin(), name.end(), ' '), name.end());

    char alternative_location = line[16];

    std::string resname = line.substr(17, 3);
    resname.erase(std::remove(resname.begin(), resname.end(), ' '), resname.end());

    char chain_id = line[21];

    int residue_id = stoi(line.substr(22, 4));

    Position<float> pos (stof(line.substr(30, 8)), stof(line.substr(38, 8)), stof(line.substr(46, 8)));
    
    float occupancy = stof(line.substr(54, 6));
    float temperature_factor = stof(line.substr(60, 6));
    std::string element = line.substr(76, 2);
    std::string charge = line.substr(78, 2);

    // float occupancy;
    // try 
    //     occupancy = stof(line.substr(54, 6));
    // catch (const std::exception& e)
    //     occupancy = 0.;
    // float temperature_factor;
    // try 
    //     temperature_factor = stof(line.substr(60, 6));
    // catch (const std::exception& e)
    //     temperature_factor = 0.;
    // std::string element;
    // if (line.length() >= 77)
    //     element = line.substr(76, 2);
    // else element = "";
    // std::string charge;
    // if (line.length() >= 78)
    //     charge = line.substr(78, 2);
    // else charge = "";

    element.erase(std::remove(element.begin(), element.end(), ' '), element.end());
    
    Atom atom(element, name, serial_id, pos);
    atom.set_record_name(record_name);
    atom.set_alternative_location(alternative_location);
    atom.set_residue_name(resname);
    atom.set_chain_id(chain_id);
    atom.set_residue_id(residue_id);
    atom.set_charge(charge);
    if ( occupancy > 0.) {
        atom.set_occupancy(occupancy);
    }
    else {
        occupancy = 1.0;
        atom.set_occupancy(occupancy);
    }
    atom.set_temperature_factor(temperature_factor);
    
    // assign atom type
    if (typization_!="-") {
        atom.assign_typization(typization_);
        //std::cout << "AAA " << typization_ << " " << atom.get_type() <<  std::endl;
    }
    
    //std::cout << line << std::endl;
    //std::cout << atom << std::endl;
    //exit(0);

    return atom;
}

std::string PDBfile::writePDBlineAtom(const Atom& a) const{
    
    char s[100];
    std::string tmp="--"; // TODO: fix charges...
    sprintf(s, "%-6s%5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
    a.get_record_name().c_str(),
    a.get_serial_id(),
    a.get_name().c_str(),
    a.get_alternative_location(),
    a.get_residue_name().c_str(),
    a.get_chain_id(),
    a.get_residue_id(),
    a.get_position().x_,
    a.get_position().y_,
    a.get_position().z_,
    a.get_occupancy(),
    a.get_temperature_factor(),
    a.get_element().c_str(),
    tmp.c_str()
    );
    //std::cout << "AAA : " << a.get_element() << std::endl;
    return s;
/*
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
     13 - 16        Atom          name         Atom name.
     17             Character     altLoc       Alternate location indicator.
     18 - 20        Residue name  resName      Residue name.
     22             Character     chainID      Chain identifier.
     23 - 26        Integer       resSeq       Residue sequence number.
     27             AChar         iCode        Code for insertion of residues.
     31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
     39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
     47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
     55 - 60        Real(6.2)     occupancy    Occupancy.
     61 - 66        Real(6.2)     tempFactor   Temperature  factor.
     77 - 78        LString(2)    element      Element symbol, right-justified.
     79 - 80        LString(2)    charge       Charge  on the atom.
     */

}

/*
void PDBfile::readPDBFile(bool flag_ter, bool flag_endmdl, bool only_heavy_atoms) {
    
    // Reading pdb file "ATOM" and "HETATM" records only.
    // If only_heavy_atoms is True, then hydrogens will be ignored.
    
    std::ifstream pdbfile(file_name_);
    if (pdbfile.is_open())
    {
        std::string line;
        
        while (!pdbfile.eof())
        {
            getline(pdbfile, line);
            if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")
            {
                
                Atom atom = readPDBlineAtom(line);
                
                if (!(only_heavy_atoms && ((atom.element() == "H") || (atom.element() == "D")))) {
                    atoms_.push_back(atom);
                }
            }
            if (line.substr(0, 6) == "ENDMDL" && flag_endmdl) break;
            if (line.substr(0, 3) == "TER" && flag_endmdl) break;
            if (line.substr(0, 3) == "END") break;
        }
        pdbfile.close();
    }
    else {
        std::cerr << "cannot to open file " + file_name_ << std::cout;
    }
    return;
}
 */



void PDBfile::filter_binding_sites() {
    
    if (ligands_map_.size()==0) {
        if (!silent_)
        std::cerr<<"warning: there are no ligands in the pdb file" << std::endl;
        return;
    }
    
    /*
    if (proteins_map_.size()==0) {
        std::cerr<<"warning: there are no proteins in the pdb file" << std::endl;
        return;
    }
    */
    if (protein_.size()<atoms_in_protein_cutoff_) {
        std::cerr << "warning: protein does not pass protein size threshold: " << file_name_ << " " << chain_id_ << " " <<  protein_.size() << std::endl;
        return;
    }
    
    std::vector<Ligand> ligands_to_delete;
    // filtering ligands by
    // - modres
    // - blacklisted
    // - size of the ligand
    unsigned int atom_count;
    for (const auto& [ligand, atoms] : ligands_map_) {

        if (std::find(modified_residues_.begin(), modified_residues_.end(), ligand) != modified_residues_.end()) {
            if (!silent_)
            std::cerr << "warning: ligand corresponds to the modified amino acid residue: " << ligand.name_ << std::endl;
            ligands_to_delete.push_back(ligand);
            continue;
        }
        if (filter_ligands_ && (std::find(not_a_ligand_.begin(), not_a_ligand_.end(), ligand.name_) != not_a_ligand_.end())) {
            if (!silent_)
            std::cerr << "warning: ligand corresponds to the disregarded ligands: " << ligand.name_ << std::endl;
            ligands_to_delete.push_back(ligand);
            continue;
        }
        atom_count = 0;
        for (const auto& atom: atoms){
            if (atom.get_alternative_location() == 'A' || atom.get_alternative_location() == ' ')
                atom_count++;
        }
        if (atom_count < atoms_in_ligand_cutoff_) {
            if (!silent_)
            std::cerr << "warning: ligand does not pass ligand size threshold: " << ligand.name_ << " " << atoms.size() << std::endl;
            ligands_to_delete.push_back(ligand);
            continue;
        }
        
    }
    
    // if all ligands should be discarded there are no binding sites
    if (ligands_map_.size()==ligands_to_delete.size()) {
        if (!silent_)
        std::cerr<<"warning: there are no relevant ligands in the pdb file" << std::endl;
        return;
    }
    else {
        for (const auto& ligand : ligands_to_delete ) {
            ligands_map_.erase(ligand);
        }
    }
    
    /*
    std::vector<char> proteins_to_delete;
    // filtering proteins py
    // - size
    for (auto protein: proteins_map_) {
        if (protein.second.size()<atoms_in_protein_cutoff_) {
            std::cerr << "warning: protein does not pass protein size threshold: " << protein.first << " " << protein.second.size() << std::endl;
            proteins_to_delete.push_back(protein.first);
            continue;
        }
    }
    
    // if all proteins should be discarded there are no binding sites
    if (proteins_map_.size()==proteins_to_delete.size()) {
        std::cerr<<"warning: there are no relevant proteins in the pdb file" << std::endl;
        return;
    }
    else {
        for (auto protein : proteins_to_delete ) {
            proteins_map_.erase(protein);
        }
    }
    */

    // determine interfaces
    Position<float> p;
    int counter = 0;
    for (const auto& [ligand, atoms] : ligands_map_) {
        auto key = ligand;
        interfaces_[key];
        for (size_t i=0; i<protein_.size(); ++i) {
            for (size_t j=0; j<atoms.size(); ++j) {
                p = protein_[i].get_position()-atoms[j].get_position();
                if ( p.get_squared_length() <= interface_squared_radius_) {
                    interfaces_[key].push_back(&protein_[i]);
                    break;
                }
            }
        }
        counter+=1;
    }
    
    ligands_to_delete.clear();
    // delete non-relevant interfaces and corresponding ligands
    for (const auto& [ligand, patoms]: interfaces_) {
        if (patoms.size() < atoms_in_interface_cutoff_) {
            ligands_to_delete.push_back(ligand);
            if (!silent_)
                std::cerr<<"warning: interface with ligand " << ligand.name_ << " has too few atoms: " << patoms.size() << std::endl;
        }
    }

    for (const auto& ligand: ligands_to_delete) {
        ligands_map_.erase(ligand);
        interfaces_.erase(ligand);
    }
    
    if (!silent_){
        std::cout << "There are " << interfaces_.size() << " ligand-binding interfaces for " << file_name_ << " chain " << chain_id_ << std::endl;
        for (const auto& interface: interfaces_) {
            printf("\t-- %lu atoms on interface with ligand %s\n", (interface.second).size(), interface.first.name_.c_str() );
        }
    }
    
    return;
}

void PDBfile::center() {
    
    // translate to COM
    Position<float> com(0,0,0);
    Position<float> p;
    for (size_t i=0; i<protein_.size(); ++i) {
        com+=protein_[i].get_position();
    }
    com/=int(protein_.size());
    
    for (size_t i=0; i<protein_.size(); ++i) {
        p = protein_[i].get_position() - com;
        protein_[i].set_position(p);
    }
    
    for (auto& ligand: ligands_map_) {
        for (size_t i=0; i<ligand.second.size(); ++i) {
            p = ligand.second[i].get_position() - com;
            ligand.second[i].set_position(p);
        }
    }
    
    // rotate w.r.t. principle axes

    // Inertia matrix
    Matrix3<double> I;
    for (size_t i=0; i<protein_.size(); ++i) {
        p = protein_[i].get_position();
        I[0][0] += (p.y_ * p.y_ + p.z_ * p.z_);
        I[1][1] += (p.x_ * p.x_ + p.z_ * p.z_);
        I[2][2] += (p.y_ * p.y_ + p.x_ * p.x_);
        I[0][1] -= p.x_ * p.y_;
        I[0][2] -= p.x_ * p.z_;
        I[1][2] -= p.y_ * p.z_;
    }
    I[1][0] = I[0][1];
    I[2][0] = I[0][2];
    I[2][1] = I[1][2];
    I /= int(protein_.size());

    // eigenvectors
    auto v = I.eigenvectors();
    // inverse (transpose because orthogonal)
    // v = v.inv();
    v = v.transposed();
    
    // rotate protein
    for (size_t i=0; i<protein_.size(); ++i){
        p = v * protein_[i].get_position();
        protein_[i].set_position(p);
    }
    // rotate ligands
    for (auto& ligand: ligands_map_){
        for (size_t i=0; i<ligand.second.size(); ++i){
            p = v * ligand.second[i].get_position();
            ligand.second[i].set_position(p);
        }
    }

    return;
}

void PDBfile::rotate(const double& theta, const double& phi, const double& psi){
    auto Q = rotation_matrix(theta, phi, psi);
    // rotate protein
    Position<float> p;
    for (size_t i=0; i<protein_.size(); ++i){
        p = Q * protein_[i].get_position();
        protein_[i].set_position(p);
    }
    // rotate ligands
    for (auto& ligand: ligands_map_){
        for (size_t i=0; i<ligand.second.size(); ++i){
            p = Q * ligand.second[i].get_position();
            ligand.second[i].set_position(p);
        }
    }
    return;
}


void PDBfile::compute_bounding_box(const float& buffer) {
    
    //printf("ZZZ: %f %f %f %f %f %f\n", x_min_,y_min_,z_min_,x_max_,y_max_,z_max_);

    x_min_ = protein_[0].get_position().x_;
    y_min_ = protein_[0].get_position().y_;
    z_min_ = protein_[0].get_position().z_;
    x_max_ = protein_[0].get_position().x_;
    y_max_ = protein_[0].get_position().y_;
    z_max_ = protein_[0].get_position().z_;

    for (int i=0; i<protein_.size(); ++i) {
        if (x_min_ > protein_[i].get_position().x_) x_min_ = protein_[i].get_position().x_;
        if (y_min_ > protein_[i].get_position().y_) y_min_ = protein_[i].get_position().y_;
        if (z_min_ > protein_[i].get_position().z_) z_min_ = protein_[i].get_position().z_;
        if (x_max_ < protein_[i].get_position().x_) x_max_ = protein_[i].get_position().x_;
        if (y_max_ < protein_[i].get_position().y_) y_max_ = protein_[i].get_position().y_;
        if (z_max_ < protein_[i].get_position().z_) z_max_ = protein_[i].get_position().z_;
    }
    
    //printf("ZZZ: %f %f %f %f %f %f\n", x_min_,y_min_,z_min_,x_max_,y_max_,z_max_);

    
    if (buffer!=0.0) {
        x_min_-=buffer;
        y_min_-=buffer;
        z_min_-=buffer;
        x_max_+=buffer;
        y_max_+=buffer;
        z_max_+=buffer;
    }
    
    //printf("ZZZ: %f %f %f %f %f %f\n", x_min_,y_min_,z_min_,x_max_,y_max_,z_max_);

    
    return;
}

void PDBfile::compute_interface_boxes(){
    Position<float> p(0, 0, 0);
    Box<float> box;
    for (const auto &[ligand, atoms] : interfaces_){
        box.center = Position<float>(0, 0, 0);
        box.size = Position<float>(0, 0, 0);
        for (size_t i=0; i<atoms.size(); ++i){
            p = atoms[i]->get_position();
            box.center += p;
        }
        box.center /= int(atoms.size());
        for (size_t i=0; i<atoms.size(); ++i){
            p = atoms[i]->get_position() - box.center;
            box.size.x_ += p.x_ * p.x_;
            box.size.y_ += p.y_ * p.y_;
            box.size.z_ += p.z_ * p.z_;
        }
        box.size /= int(atoms.size());
        box.size.x_ = 2 * std::sqrt(box.size.x_);
        box.size.y_ = 2 * std::sqrt(box.size.y_);
        box.size.z_ = 2 * std::sqrt(box.size.z_);
        interface_boxes_[ligand] = box;
    }
    return;
}

void PDBfile::writePDBfile(const std::string& filename, 
    const bool& flag_write_protein, const bool& flag_write_ligands) const{
    std::ofstream f;
    if (!flag_write_protein && !flag_write_ligands){
        std::cerr << "pass at least one flag to write pdb" << std::endl;
        return;
    }
    f.open(filename);
    if (add_header.size() > 0){
        f << add_header << std::endl;
    }
    f << "REMARK   FILE GENERATED BY  ?descriptors_" << std::endl;
    if (flag_write_protein){
        for (size_t i=0; i<protein_.size(); ++i){
            f << writePDBlineAtom(protein_[i]) << std::endl;
        }
    }
    if (flag_write_ligands){
        for (auto& ligand: ligands_map_){
            if (interfaces_.find(ligand.first) != interfaces_.end()){
                for (size_t i=0; i<ligand.second.size(); ++i){
                    f << writePDBlineAtom(ligand.second[i]) << std::endl;
                }
            }
        }
    }
    f.close();
    return;
}


void PDBfile::print_interface(const Ligand& l, const std::string& fn, const bool capped, const bool flag_split_protein_ligand) {
    
    // Print pdb file with ligand and neigbouring residues
    // capped option mimics ACE and CT3 caps by taking atoms from previous and next residues, respectively 

    std::set<int> residue_ids;

    std::ofstream f, fp, fl;
    if (flag_split_protein_ligand) {
        fp.open(fn+"_protein.pdb");
        fl.open(fn+"_ligand.pdb");
    }
    else {
        f.open(fn+".pdb");
    }
    

    std::string res_name; // tmp name to use in cycle
    std::string atm_name;
    int min_atom_number, max_atom_number; // min and max atom numbers, need for cap
    
    for (const auto& pa : interfaces_[l]) {
        residue_ids.insert(pa->get_residue_id());
    }
    
    // i-th and i+2th residues will have overlapped caps, e.g. CAT=CAY
    // so add i+1th residue
    std::set<int> additional_residues;
    for (const auto ind : residue_ids) {
        if (residue_ids.find(ind+2)!=residue_ids.end()) {
            if (residues_.find(ind+1)!=residues_.end()) {
                additional_residues.insert(ind+1);
            }
        }
    }
    
    residue_ids.insert(additional_residues.begin(), additional_residues.end());
    //std::cout << "INTERFACE : ";
    //for (const auto ind : residue_ids) {
    //    std::cout << ind << " ";
    //}
    //std::cout << std::endl;


    for (const auto& id : residue_ids) {

        res_name = residues_[id].atoms_[0]->get_residue_name();
        std::set<int> ids;
        for (const auto& a : residues_[id].atoms_) {
            ids.insert(a->get_serial_id());
        }
        min_atom_number = (*ids.begin());
        max_atom_number = (*prev(ids.end()));

        // if left residue is not on the interface skip
        if  ( (capped==false || residue_ids.find(id-1)!=residue_ids.end()) == false ) {

            // if left residue exists
            if (residues_.find(id-1)!=residues_.end()) {

                for (const auto& pa: residues_[id-1].atoms_) {
                    atm_name = pa->get_name();
                    if (atm_name=="CA" || atm_name=="C" || atm_name=="O") {
                        Atom a_tmp = (*pa);
                        a_tmp.set_residue_name(res_name);
                        a_tmp.set_residue_id(id);
                        a_tmp.set_name(atm_name+"Y");
                        if (atm_name=="CA") {
                            a_tmp.set_serial_id(min_atom_number-3);
                        }
                        else if (atm_name=="O") {
                            a_tmp.set_serial_id(min_atom_number-1);
                        }
                        else if (atm_name=="C") {
                            a_tmp.set_serial_id(min_atom_number-2);
                        }


                        if (flag_split_protein_ligand) {
                            fp << writePDBlineAtom(a_tmp) << '\n';
                        }
                        else {
                            f << writePDBlineAtom(a_tmp) << '\n';
                        }
                    }
                }
            }
        }

        // write residue itself
        for (const auto& pa : residues_[id].atoms_) {
            if (flag_split_protein_ligand) {
                fp << writePDBlineAtom((*pa)) << '\n';
            }
            else {
                f << writePDBlineAtom((*pa)) << '\n';
            }
            
        }

        // if right residue is not on the interface skip
        if ( (capped==false || residue_ids.find(id+1)!=residue_ids.end()) == false ) {

            // if right residue exists
            if (residues_.find(id+1)!=residues_.end()) {

                for (const auto& pa: residues_[id+1].atoms_) {
                    atm_name = pa->get_name();
                    if (atm_name=="N" || atm_name=="CA") {
                        Atom a_tmp = (*pa);
                        a_tmp.set_residue_name(res_name);
                        a_tmp.set_residue_id(id);
                        a_tmp.set_name(atm_name+"T");
                        if (atm_name=="CA") {
                            a_tmp.set_serial_id(max_atom_number+2);
                        }
                        else if (atm_name=="N") {
                            a_tmp.set_serial_id(max_atom_number+1);
                        }


                        if (flag_split_protein_ligand) {
                            fp << writePDBlineAtom(a_tmp) << '\n';
                        }
                        else {
                            f << writePDBlineAtom(a_tmp) << '\n';
                        }
                    }
                }
            }
        }
    }

    for (const auto& a : ligands_map_[l]) {

        if (flag_split_protein_ligand) {
            fl << writePDBlineAtom(a) << '\n';
        }
        else {
            f << writePDBlineAtom(a) << '\n';
        }
    }

    /*
    // TODO: manage connect records
    std::string s;
    for (const auto& v : l.conects_) {
        f<<"CONECT";        
        for (const auto& i : v) {
            char stmp[6];
            sprintf(stmp, "%5d", i);
            f << stmp;
        }
        f<<'\n';
    }
    */
    
    if (flag_split_protein_ligand) {
        fp.close();
        fl.close();
    }
    else {
        f.close();
    }

    return;
}

float* PDBfile::export_interfaces() {
    if (export_res_){
        delete [] export_res_;
        export_res_ = NULL;
    }
    export_res_ = new float[interface_boxes_.size() * 6];
    unsigned int index = 0;
    for (const auto &[ligand, box] : interface_boxes_){
        export_res_[index * 6 + 0] = box.center.x_;
        export_res_[index * 6 + 1] = box.center.y_;
        export_res_[index * 6 + 2] = box.center.z_;
        export_res_[index * 6 + 3] = box.size.x_;
        export_res_[index * 6 + 4] = box.size.y_;
        export_res_[index * 6 + 5] = box.size.z_;
        index++;
    }
    return export_res_;
}

std::ostream& operator<< (std::ostream& os, const PDBfile& pdb) {
    os << "file : " << pdb.get_file_name() << " " 
    << "chain : " << pdb.get_chain_id() << " "
    << "# protein atoms : " << pdb.protein_.size() << " "
    << "ligands :" ;
    for (const auto& [ligand, atoms] : pdb.ligands_map_) {
        os << "\n\tname:" << ligand << " atoms: " << atoms.size();
    }
    return os;
}