//
//  ConvertTrajectoryType.hpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_ConvertTrajectoryType_hpp
#define LiquidLib_ConvertTrajectoryType_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class ConvertTrajectoryType : public Trajectory {
public:
    ConvertTrajectoryType();
    virtual ~ConvertTrajectoryType();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void convert_trajectory();
    
    
protected:
    //protected member functions
    void check_parameters() throw();
    inline void print_status(size_t & status);
    void convert_to_xtc();
    void convert_to_xyz();
    // TODO
    //void convert_to_XDATCAR();
    //void convert_to_trr();
    
    string outfile_type_;
    string input_file_name_;
    string output_file_name_;
    
private:
};

#endif // defined (LiquidLib_MeanSquaredDisplacement_hpp)
