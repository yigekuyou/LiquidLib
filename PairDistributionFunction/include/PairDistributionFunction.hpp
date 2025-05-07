//
//  PairDistributionFunction.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_PairDistributionFunction_hpp
#define LiquidLib_PairDistributionFunction_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class PairDistributionFunction : public Trajectory {
public:
    PairDistributionFunction();
    virtual ~PairDistributionFunction();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_g_r();
    void write_g_r();

protected:
    
private:
// private member function
    void check_parameters() throw();
    inline void histogram_g_r(size_t const & frame_number,
                              size_t const & atom1_index, size_t const & atom2_index,
                              size_t const & atom1_type, size_t const & atom2_type,
                              double const & delta_r);
    
    void determine_atom_indexes(vector < string > const & atom_types,
                                vector < double > const & scattering_lengths,
                                string            const & atom_group,
                                vector < vector < unsigned int > > & atom_types_indexes,
                                double & average_scattering_length,
                                size_t & number_of_atoms);
    
    bool check_atom_types_equivalent();
    inline void print_status(size_t & status);
    
// private member variables
    string input_file_name_;
    string output_file_name_;
    
    vector < string > atom_types1_;
    vector < double > scattering_lengths1_;
    string atom_group1_;

    vector < string > atom_types2_;
    vector < double > scattering_lengths2_;
    string atom_group2_;

    unsigned int number_of_bins_;
    unsigned int number_of_frames_to_average_;
	
	double max_cutoff_length_;
	
    vector< double > r_values_;
    vector< double > g_r_;
};

#endif // defined (LiquidLib_PairDistributionFunction_hpp)
