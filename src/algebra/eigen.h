/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

// Eigen
#define EIGEN_DEFAULT_INDEX_TYPE int
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wdocumentation"
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wnonportable-include-path"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#pragma GCC diagnostic pop
typedef double real_t;
typedef EIGEN_DEFAULT_INDEX_TYPE index_t;

// shortnames

typedef Eigen::Matrix<index_t,-1,-1> IMatrix;
typedef Eigen::Map<IMatrix>          IMap;
typedef Eigen::Map<const IMatrix>    CIMap;

typedef Eigen::Matrix<real_t,-1,-1>  FMatrix;
typedef Eigen::Map<FMatrix>          FMap;
typedef Eigen::Map<const FMatrix>    CFMap;

typedef Eigen::Array<real_t,-1,-1>   FArray;
typedef Eigen::Map<FArray>           FAMap;
typedef Eigen::Map<const FArray>     CFAMap;

typedef Eigen::Matrix<real_t,-1,1>   FVector;
typedef Eigen::Map<FVector>          FVMap;
typedef Eigen::Map<const FVector>    CFVMap;

typedef Eigen::SparseMatrix<real_t,Eigen::RowMajor,index_t> SMatrix;
typedef Eigen::Map<SMatrix>         SMap;
typedef Eigen::Map<const SMatrix>   CSMap;


typedef Eigen::PermutationMatrix<-1,-1, index_t> Perm;

typedef Eigen::Triplet<real_t,index_t> Triplet;
typedef std::vector<Triplet>           TripletVec;


// input and output
#include <json.hpp>
using Json = nlohmann::json;

namespace Eigen {

void to_json  (Json& j, const FMatrix& p);
void to_json  (Json& j, const FMap& p);
void to_json  (Json& j, const CFMap& p);

void to_json  (Json& j, const IMatrix& p);
void to_json  (Json& j, const IMap& p);
void to_json  (Json& j, const CIMap& p);

void to_json  (Json& j, const SMatrix& p);

void from_json(const Json& j, FMatrix& p);
void from_json(const Json& j, IMatrix& p);
void from_json(const Json& j, SMatrix& p);
}
