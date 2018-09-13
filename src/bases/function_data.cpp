/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <sstream>

#include <bases/function_data.h>

// BasisValues

index_t BasisValues::numData  () const {return valueVec.size();}


SMap  BasisValues::operator[] (PartialDerivative::DirectionDer der)
{ return operator[](getIndex(PartialDerivative(der)));}
CSMap BasisValues::operator[] (PartialDerivative::DirectionDer der) const
{ return operator[](getIndex(PartialDerivative(der)));}

SMap  BasisValues::operator[] (PartialDerivative der)
{ return operator[](getIndex(der));}
CSMap BasisValues::operator[] (PartialDerivative der) const
{ return operator[](getIndex(der));}


SMap BasisValues::operator[] (index_t i)
{
    return valueVec.size()>static_cast<size_t>(i) && valueVec[i].size() ? Sparsity::matrix(valueVec[i].begin()):  SMap(0,0, 0,nullptr,nullptr,nullptr);
}

CSMap BasisValues::operator[] (index_t i) const
{
    return valueVec.size()>static_cast<size_t>(i) && valueVec[i].size() ? Sparsity::matrix(valueVec[i].begin()):  CSMap(0,0, 0,nullptr,nullptr,nullptr);
}


void  BasisValues::matrix (index_t der, SMap&dst)
{
    if (valueVec.size()>static_cast<size_t>(der) )
        Sparsity::matrix(valueVec[der].begin(), dst);
    else
        new (&dst)SMap(0,0, 0,nullptr,nullptr,nullptr);
}

void  BasisValues::matrix (index_t der, CSMap&dst) const
{
   if (valueVec.size()>static_cast<size_t>(der) )
        Sparsity::matrix(valueVec[der].begin(), dst);
    else
        new (&dst)CSMap(0,0, 0,nullptr,nullptr,nullptr);
}

void  BasisValues::matrix (PartialDerivative der, SMap&dst)                      {matrix(getIndex(der),dst);}
void  BasisValues::matrix (PartialDerivative der, CSMap&dst) const               {matrix(getIndex(der),dst);}
void  BasisValues::matrix (PartialDerivative::DirectionDer der, SMap&dst)        {matrix(getIndex(PartialDerivative(der)),dst);}
void  BasisValues::matrix (PartialDerivative::DirectionDer der, CSMap&dst) const {matrix(getIndex(PartialDerivative(der)),dst);}


index_t BasisValues::getIndex (PartialDerivative der) const
{ return std::find(derVec.begin(),derVec.end(),der) -derVec.begin(); }

const real_t *BasisValues::data(index_t d) const
{ return valueVec.size()>static_cast<size_t>(d) ? valueVec[d].begin() :nullptr; }

real_t *BasisValues::data(index_t d)
{ return valueVec.size()>static_cast<size_t>(d) ? valueVec[d].begin() :nullptr; }

const real_t *BasisValues::data(index_t d, index_t r) const
{ return valueVec.size()>static_cast<size_t>(d) ? valueVec[d].begin()+start(r) :nullptr; }

real_t *BasisValues::data(index_t d, index_t r)
{ return valueVec.size()>static_cast<size_t>(d) ? valueVec[d].begin()+start(r) :nullptr; }



void initBasisValues(BasisValues& val,
                     int num_derivs,
                     size_t rows,
                     size_t cols,
                     index_t nnzs)
{
    size_t size=rows*nnzs;
    val.m_cols=cols;
    val.valueVec.resize(num_derivs);
    
    val.m_idxs.swap(OwningView<index_t>(size));
    val.m_start.swap(OwningView<index_t>(rows+1));
    for (size_t i=0; i<rows+1;++i)
        val.m_start[i]=(nnzs)*i;
    
    val.m_data.resize(size*num_derivs);
    for(int i=0;i<num_derivs;++i)
    { val.valueVec[i]=view<real_t>(val.m_data.data()+size*i,size); }
}

std::vector<PartialDerivative> gradient(int domain)
{
   std::vector<PartialDerivative> der(domain,0);
    for (int i=0; i<domain ; ++i )
        der[i][i]=1;
    return der;
}

std::vector<PartialDerivative> hessian (int domain)
{
   std::vector<PartialDerivative> der(domain,0);
    for (int i=0; i<domain ; ++i )
        der[i][i]=1;
    return der;
}

void to_json(Json& j, const PartialDerivative & pd)
{
    std::stringstream s;
    s<<std::oct<<pd.raw;
    j=s.str();
}
