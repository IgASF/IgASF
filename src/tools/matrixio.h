/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include<vector>
#include<unistd.h>
#include<algebra/eigen.h>

typedef Eigen::SparseMatrix<double,Eigen::RowMajor> RowM;
typedef Eigen::SparseMatrix<double,Eigen::ColMajor> ColM;
typedef Eigen::Map<const RowM>                      RMap;
typedef Eigen::Map<const ColM>                      CMap;

struct inData
{
    int   rowMajor;
    int   rows;
    int   cols;
    std::vector<int>    innerBeg;
    std::vector<int>    innerPos;
    std::vector<double> values;

    RMap asRowMatrix() const;
    CMap asColMatrix() const;

    bool    valid() const;
    double  norm()  const;
};

struct outData
{
    outData(const RowM &mat );
    outData(RMap mat );
    outData(const ColM &mat );
    outData(CMap mat );

    int     rowMajor;
    int     rows;
    int     cols;
    int     nnzs;
    const int*    innerBeg;
    const int*    innerPos;
    const double* values;
};

void     writeMatrix(const outData &data, const std::string& fn);
inData   readMatrix (const std::string& fn);






