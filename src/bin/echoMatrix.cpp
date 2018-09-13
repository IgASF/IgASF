/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <exception>
#include <iostream>

#include <tools/matrixio.h>


void print_help(char** args)
{
    std::cout<<"\n"<<args[0]
    <<" mat\n\n"
    <<"  mat       name of the matrix file to be echoed to the console; use\n"
    <<"            \"stdin\" to read from console.\n"
    <<std::endl;
}

int main (int argn, char** args)
{
    if(argn!=2)
        {
        if(argn<2)
            std::cout << "\nThe following error occured: Input missing."<< std::endl;
        else
            std::cout << "\nThe following error occured: Too much input."<< std::endl;            
        print_help(args);
        return 1;
        }

    inData matrix;

    try {
        matrix = readMatrix(args[1]);
        }
    catch (std::exception& e) {
        std::cout << "There was a problem in reading the file: " << e.what() << "\n" << std::endl;
        return 3;
    }

    switch(matrix.rowMajor)
    {
        case 0: std::cout << FMatrix(matrix.asColMatrix()) << std::endl; break;
        case 1: std::cout << FMatrix(matrix.asRowMatrix()) << std::endl; break;
    }

    return 0;
}
