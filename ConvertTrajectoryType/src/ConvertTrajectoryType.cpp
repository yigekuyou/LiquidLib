//
//  ConvertTrajectoryType.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//

#include "ConvertTrajectoryType.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <iomanip>
#include <cstring>
#include <sstream>

#include "Trajectory.hpp"
#ifdef GROMACS
extern "C" {
    #include "xdrfile.h"
    #include "xdrfile_trr.h"
    #include "xdrfile_xtc.h"
}
#endif

using namespace std;

ConvertTrajectoryType::ConvertTrajectoryType() :
    outfile_type_(""),
    input_file_name_("convert_trajectory.in"),
    output_file_name_("")
{
}


ConvertTrajectoryType::~ConvertTrajectoryType()
{
}


void ConvertTrajectoryType::read_command_inputs(int argc, char * argv[])
{
    for (int input = 1; input < argc; ++input) {
        if (strcmp(argv[input], "-i") == 0) {
            input_file_name_ = argv[++input];
            continue;
        }
        if (strcmp(argv[input], "-o") == 0) {
            output_file_name_ = argv[++input];
            continue;
        }
        if (strcmp(argv[input], "-t") == 0) {
            trajectory_file_name_ = argv[++input];
            continue;
        }
        if (strcmp(argv[input], "-v") == 0) {
            is_run_mode_verbose_ = 1;
            continue;
        }
        
        cerr << "\nERROR: Unrecognized flag '" << argv[input] << "' from command inputs.\n";
        exit(1);
    }
}


void ConvertTrajectoryType::read_input_file()
{
    ifstream input_file(input_file_name_);
    
    if (!input_file) {
        cerr << "ERROR: Input file location, ";
        cerr << "\033[1;25m";
        cerr << input_file_name_;
        cerr << "\033[0m";
        cerr << ", does not exist, please check input";
        cerr << endl;
        exit(1);
    }
    
    string input_word;
    
    while (input_file >> input_word) {
        // check for comment
        if (input_word[0] == '#') {
            getline(input_file, input_word);
            continue;
        }
        
        // check for memeber bools
        if (input_word == "is_run_mode_verbose") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            if (input_word == "true" || input_word == "yes") {
                is_run_mode_verbose_ = true;
            }
            else if(input_word == "false" || input_word == "no") {
                is_run_mode_verbose_ = false;
            }
            else {
                is_run_mode_verbose_ = stoi(input_word);
            }
            continue;
        }
        if (input_word == "is_wrapped") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            if (input_word == "true" || input_word == "yes") {
                is_wrapped_ = true;
            }
            else if(input_word == "false" || input_word == "no") {
                is_wrapped_ = false;
            }
            else {
                is_wrapped_ = stoi(input_word);
            }
            continue;
        }
        
        // check if equal to member ints
        if (input_word == "start_frame") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            start_frame_ = stoi(input_word);
            continue;
        }
        if (input_word == "end_frame") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            end_frame_ = stoi(input_word);
            continue;
        }
        if (input_word == "dimension") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            dimension_ = stoi(input_word);
            continue;
        }
        if (input_word == "frame_interval") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            frame_interval = stoi(input_word);
            continue;
        }
        
        // check if equal to member doubles
        if (input_word == "trajectory_delta_time") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_delta_time_ = stod(input_word);
            continue;
        }
        
        //check if equal to member strings
        if (input_word == "output_file_name") {
            if (output_file_name_ != "") {
                cerr << "ERROR: Please do not set output file by command line and input file,\n";
                cerr << "     : we are unsure on which to prioritize.";
                cerr << endl;
                exit(1);
            }
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            output_file_name_ = input_word;
            continue;
        }
        if (input_word == "trajectory_file_name") {
            if (trajectory_file_name_ != "") {
                cerr << "ERROR: Please do not set trajectory file by command line and input file,\n";
                cerr << "     : we are unsure on which to prioritize.";
                cerr << endl;
                exit(1);
            }
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_file_name_ = input_word;
            continue;
        }
        if (input_word == "trajectory_file_type") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_file_type_ = input_word;
            continue;
        }
#ifdef GROMACS
        if (input_word == "gro_file_name") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            gro_file_name_ = input_word;
            continue;
        }
#else
        if (input_word == "gro_file_name") {
            cerr << "gro files cannot be used in non gromacs";
            cerr << "compatible version of LiquidLib\n";
            cerr << endl;
        }
#endif
        if (input_word == "outfile_type") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            outfile_type_ = input_word;
            continue;
        }
        if (input_word == "output_precision") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            output_precision_ = stod(input_word);
            continue;
        }
        if (input_word == "boxlength") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                average_box_length_.push_back(stod(input_word));
            }
            continue;
        }
        
        //check for everything else
        cerr << "WARNING: no matching input type for: ";
        cerr << "\033[1;33m";
        cerr << input_word;
        cerr << "\033[0m";
        cerr << " disregarding this variable and continueing to next line\n";
        getline(input_file, input_word);
    }
    check_parameters();
    
    input_file.close();
}


void ConvertTrajectoryType::convert_trajectory()
{
    if (outfile_type_ == "xtc") {
        convert_to_xtc();
    }
    else if (outfile_type_ == "xyz") {
        convert_to_xyz();
    }
    else {
        cerr << "ERROR:  No matching file type to convert to";
    }
}

void ConvertTrajectoryType::convert_to_xyz()
{
    ofstream xyz_file(output_file_name_);
    
    size_t status = 0;
    
    for (int i_frame = 0; i_frame < (end_frame_ - start_frame_)/frame_interval; ++i_frame) {
        xyz_file << number_of_system_atoms_ << "\n";
        xyz_file << "unused line" << "\n";
        for (int i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            xyz_file << "  ";
            xyz_file << system_atom_types_[i_atom] << "\t";
            
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                xyz_file << trajectory_[i_frame][i_atom][i_dimension] << "\t";
            }
            xyz_file << "\n";
        }
        
        if (is_run_mode_verbose_) {
            print_status(status);
        }
    }
    xyz_file.close();
}


void ConvertTrajectoryType::convert_to_xtc()
{
#ifdef GROMACS
    float precision = 8;
    
    string gro_file_path = output_file_name_.substr(0, output_file_name_.find_last_of(".")) + ".gro";
    
    rvec *configuration;
    configuration = new rvec[number_of_system_atoms_];
    
    matrix box;
    box[0][0] = average_box_length_[0];
    box[1][1] = average_box_length_[1];
    box[2][2] = average_box_length_[2];
    
    // Create gro file
    ofstream gro_file(gro_file_path);
    gro_file << "Generic GRO File" << "\n";
    gro_file << " " << number_of_system_atoms_ << "\n";
    
    for (int i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
        // print molecuele number
        gro_file << setw(5) << setfill(' ') << 1;
        
        // print molecule name
        gro_file << setw(5) << setfill(' ') << left << "UNK" << right;
        
        // print atom type
        gro_file << setw(5) << setfill(' ') << system_atom_types_[i_atom];
        
        // print atom number
        gro_file << setw(5) << setfill(' ') << i_atom+1;
        
        
        // print coordinates
        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            gro_file << setw(8) << setfill(' ') << setprecision(3) << fixed;
            gro_file << trajectory_[0][i_atom][i_dimension];
        }
        gro_file << endl;
    }
    
    gro_file.close();
    
    XDRFILE *output_file = xdrfile_open(output_file_name_.c_str(), "w");
    
    size_t status = 0;
    
    for (int i_frame = 0; i_frame < (end_frame_ - start_frame_)/frame_interval; ++i_frame) {
        for (int i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                configuration[i_atom][i_dimension] = trajectory_[i_frame][i_atom][i_dimension];
            }
        }
        
        write_xtc(output_file, number_of_system_atoms_, i_frame, i_frame*trajectory_delta_time_, box, configuration, precision);
        
        if (is_run_mode_verbose_) {
            print_status(status);
        }
    }
    xdrfile_close(output_file);
#else
    cerr << "cannot convert to xtc without GROMACS xdrfile library" << endl;
    exit(0);
#endif
}


// Function to check that all the parameters provided by the user for
// mean squared displacement are usable.  Ensures that the code will
// not have with needing enough data points
void ConvertTrajectoryType::check_parameters() throw()
{
    if (end_frame_ == 0) {
        cerr << "ERROR: We require more information to proceed, either end_frame or number_of_time_points\n";
        cerr << "       must be povided for us to continue.";
        cerr << endl;
        exit(1);
    }
    
    if (is_wrapped_) {
        cerr << "WARNING: the trajectory provided is not unwrapped\n";
        cerr << "       : We will unwrapp it for you, but user discretion\n";
        cerr << "       : is advised";
        cerr << endl;
    }
    
    if (output_file_name_ == "") {
        output_file_name_ = trajectory_file_name_.substr(0, output_file_name_.find_last_of(".")) + outfile_type_;
    }
    
    // check output file can be opened
    ofstream output_r2_t_file(output_file_name_);
    if (!output_r2_t_file) {
        cerr << "ERROR: Output file: "
        << "\033[1;25m"
        << output_file_name_
        << "\033[0m"
        << ", could not be opened"
        << endl;
        exit(1);
    }
    output_r2_t_file.close();
}


void ConvertTrajectoryType::print_status(size_t & status)
{
    ++status;
    cout << "\rcurrent progress of calculating the mean squared displacement is: ";
    cout << status * 100.0/end_frame_;
    cout << " \%";
    cout << flush;
}