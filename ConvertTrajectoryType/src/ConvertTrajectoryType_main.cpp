//
//  ConvertTrajectoryType_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "ConvertTrajectoryType_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "ConvertTrajectoryType.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
    
    ConvertTrajectoryType convert_trajectory;
    
    convert_trajectory.read_command_inputs(argc, argv);
    convert_trajectory.read_input_file();
    convert_trajectory.read_trajectory();
    convert_trajectory.convert_trajectory();
    
    cout << "Successfully computed mean squared displacement\n";
    cout << "Rolling Out.";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--          Convert Trajectroy Type         --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: convert_trajectory.in)\n";
    cout << "\n";
}