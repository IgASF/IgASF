/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/sparsity.h>
#include <tools/view.hpp>
#include <map>

struct PartialDerivative
{
    PartialDerivative(uint32_t dat) : raw(dat) {}
    PartialDerivative() : raw(0) {}
    
    struct DirectionDer
    {
    DirectionDer( uint32_t &raw, int dir) : m_data(raw), m_dir(dir) {}
    void operator=(int val)
    { m_data =  (m_data &~(0xF<<(4*m_dir))) | ((val&0xF)<<(4*m_dir));}
    operator PartialDerivative() const
    {return PartialDerivative((m_data>>(4*m_dir)) &0xF);}
    explicit operator int() const {return (m_data>>(4*m_dir)) &0xF;}
    
    private:
    uint32_t &m_data;
    int       m_dir;
    };
    
    DirectionDer operator[] (int dir)       { return DirectionDer(raw,dir);}
    const DirectionDer operator[] (int dir) const
    { return DirectionDer(const_cast<uint32_t &>(raw),dir); }
    bool operator <  (const PartialDerivative &other) const
    {return raw<other.raw;}
    bool operator == (const PartialDerivative &other) const
    {return raw==other.raw;}
    
    uint32_t raw;
};
void to_json(Json& j, const PartialDerivative & pd);

std::vector<PartialDerivative> gradient(int domain);
std::vector<PartialDerivative> hessian (int domain);


struct BasisValues : public Sparsity// collections of the data of a 1D basis at all quad points
{
public:
    BasisValues() {}
    
    index_t numData  () const;
    index_t getIndex (PartialDerivative der) const ;
    SMap  operator[] (index_t der);
    CSMap operator[] (index_t der) const ;
    SMap  operator[] (PartialDerivative der);
    CSMap operator[] (PartialDerivative der) const ;
    SMap  operator[] (PartialDerivative::DirectionDer der);
    CSMap operator[] (PartialDerivative::DirectionDer der) const ;
    
    void  matrix (index_t der, SMap&dst);
    void  matrix (index_t der, CSMap&dst) const ;
    void  matrix (PartialDerivative der, SMap&dst);
    void  matrix (PartialDerivative der, CSMap&dst) const ;
    void  matrix (PartialDerivative::DirectionDer der, SMap&dst);
    void  matrix (PartialDerivative::DirectionDer der, CSMap&dst) const ;
    
    const real_t *data(index_t der) const;
          real_t *data(index_t der);
    const real_t *data(index_t der, index_t row) const;
          real_t *data(index_t der, index_t row);
private:
    friend void initBasisValues(BasisValues&, int, size_t, size_t, index_t);
    std::vector<real_t>            m_data;
    std::vector<view<real_t>>      valueVec; // coefficients
public:
    std::vector<PartialDerivative> derVec;   // index of computed values
};

// init a basis values with constant number of nnzs per row
void initBasisValues(BasisValues& val,
                       int num_derivs,
                       size_t rows,
                       size_t cols,
                       index_t nnzs);


