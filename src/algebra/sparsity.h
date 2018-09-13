/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/eigen.h>
#include <tools/view.hpp>

// describes the position of the nonZeros inside the matrix
// M_{i,j} = Sum_k Ind_{suppB_i}(x_k)Ind_{supp C_j}(x_k)
// This is used for mapping the computed vector of coefficients
// to a sparse matrix.
struct Sparsity
{
    OwningView<index_t> m_start;
    OwningView<index_t> m_idxs;
    index_t             m_cols;
    
    Sparsity () = default;
    Sparsity (const Sparsity &)= default;
    Sparsity (Sparsity &&)= default;
    Sparsity & operator=(const Sparsity &)= default;
    Sparsity & operator=( Sparsity &&)= default;


    Sparsity (index_t rows, index_t cols, index_t nnzs=0)
    : m_start(rows+1), m_idxs(nnzs), m_cols(cols)
    {}
    Sparsity (index_t rows, index_t cols, const std::vector<index_t> & starts, const std::vector<index_t> & idxs )
    : m_start(starts.begin(), starts.end()), m_idxs(idxs.begin(),idxs.end()), m_cols(cols)
    {
    (void)rows;
    assert(static_cast<index_t>(starts.size())==rows+1);
    }
    Sparsity (index_t rows, index_t cols, std::vector<index_t> &&starts, std::vector<index_t> &&idxs )
    : m_start(starts.begin(), starts.end()), m_idxs( idxs.begin(), idxs.end()), m_cols( cols)
    {
    (void)rows;
    assert(static_cast<index_t>(starts.size())==rows+1);
    }
    Sparsity (index_t rows, index_t cols, OwningView<index_t> &&starts, OwningView<index_t> &&idxs )
    : m_start(std::move(starts)), m_idxs(std::move(idxs)), m_cols( cols)
    {
    (void)rows;
    assert(static_cast<index_t>(starts.size())==rows+1);
    }

    
    // global properties
    index_t cols() const ;
    index_t rows() const ;
    index_t nnzs() const ;
    // per row properties
    index_t nnzs (index_t r) const ; // non zeros in a row
    index_t start(index_t r) const;
    // access the position of the nnz in row r
    view<const index_t> rowPos  (index_t r) const;
          view<index_t> rowPos  (index_t r);
    
    SMap   matrix (real_t* data)        ;
    CSMap  matrix (const real_t* data) const ;
    
    void   matrix (real_t* data, SMap&dst)        ;
    void   matrix (const real_t* data, CSMap&dst) const ;
    
    // returns the row and column index
    // of a nonZero position and viceversa
    std::array<index_t, 2> posToRC (index_t pos) const;
    index_t   RCToPos (index_t row, index_t col) const;
    IMatrix   posToRCMap () const;
    

};

// returns the sparsity of test^T trial if both are row majors
Sparsity bilinearSparsity (const Sparsity &test,  const Sparsity &trial);
// returns the sparsity of op1 Kronecker op2 trial if both are row majors
Sparsity kroneckerSparsity(const Sparsity &op1, const Sparsity &op2);


