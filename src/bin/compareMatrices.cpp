/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <exception>
#include <iostream>
#include <array>

#include <tools/matrixio.h>

void print_help(char** args)
{
    std::cout<<"\n"<<args[0]
    <<" mat1 mat2\n\n"
    <<"  mat1 mat2    name of the matrixs file to be compared.\n"
    <<std::endl;
}

int main (int argn, char** args)
{
    if(argn!=3)
        {
        if(argn<3)
            std::cout << "\nThe following error occured: Input missing."<< std::endl;
        else
            std::cout << "\nThe following error occured: Too much input."<< std::endl;            
        print_help(args);
        return 1;
        }

    std::array<inData,2> matrix;

    for (int i=0; i<2; ++i)
    {
        try {
            matrix[i] = readMatrix(args[i+1]);
            }
        catch (std::exception& e) {
            std::cout << "There was a problem in reading the "
                << ( i==0 ? "first" : "second" ) << " file: " << e.what() << "\n" << std::endl;
            return 3;
        }
    }

    if (matrix[0].rows!=matrix[1].rows
          || matrix[0].cols!=matrix[1].cols
          || matrix[0].rows<=0
          || matrix[0].cols<=0
    )
    {
        std::cout << "The sizes of the matrices do not agree: "
            << matrix[0].rows << "x" << matrix[0].cols << " vs. "
            << matrix[1].rows << "x" << matrix[1].cols << "\n" << std::endl;
        return 1;
    }
    
    double err;
    switch(matrix[0].rowMajor+matrix[1].rowMajor*2)
    {
        case 0: err = (matrix[0].asColMatrix()-matrix[1].asColMatrix()).norm(); break;
        case 1: err = (matrix[0].asRowMatrix()-RowM(matrix[1].asColMatrix())).norm(); break;
        case 2: err = (matrix[0].asColMatrix()-ColM(matrix[1].asRowMatrix())).norm(); break;
        case 3: err = (matrix[0].asRowMatrix()-matrix[1].asRowMatrix()).norm(); break;
    }

    std::cout << err << std::endl;

    return (err < 1.e-7) ? 0 : 1;
}
