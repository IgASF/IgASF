/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/eigen.h>

template <int DIM>
struct MultiIndex
{
    Eigen::Matrix<index_t,DIM,1> m_data;
    
    MultiIndex(IMatrix m)
    {
        m_data.resize(m.rows());
        m_data(0,0)=1;
        for (int r=1;r<m_data.rows();++r )
            m_data(r)=m_data(r-1)*m(r-1);
    }
    
    template <typename T>
    index_t multi2flat(T multi)
    { return m_data.dot(multi);}
    
    Eigen::Matrix<index_t,DIM,1> flat2multi(index_t flat)
    {
        Eigen::Matrix<index_t,DIM,1> res(m_data.rows(),1);
        for (index_t r=m_data.rows()-1; r>=0; --r)
        {
            res(r)=flat/m_data(r);
            flat=flat%m_data(r);
        }
        return res;
    }
    
    struct SubDomain
    {
        index_t shift;
        Eigen::Matrix<index_t,DIM,2> m_data;

        SubDomain(MultiIndex m, IMatrix s)
        {
            shift=m.m_data.dot(s.col(0));
            m_data.col(0)=m.m_data;
            m_data.col(1)=MultiIndex(s.col(1)-s.col(0)).m_data;
        }

        index_t local2global (index_t flat)
        {
            index_t res=0;
            for (index_t r=m_data.rows()-1; r>=0; --r)
            {
                res += (flat/m_data(r,1))*m_data(r,0);
                flat = flat%m_data(r,1);
            }
            return res+shift;
        }
        index_t global2local (index_t flat)
        {
            flat-=shift;
            index_t res=0;
            for (index_t r=m_data.rows()-1; r>=0; --r)
            {
                res += (flat/m_data(r,0))*m_data(r,1);
                flat = flat%m_data(r,0);
            }
            return res;
        }
    };
    SubDomain subdomain(IMatrix m) {return SubDomain(*this,m);}
};
