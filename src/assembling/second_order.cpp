/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <assembling/second_order.h>
#include <maps/transform_coefs.h>


bool EqCoef::hasA() const {return A.rows()==dim && A.cols()==dim && A.norm()!=0;}
bool EqCoef::hasB() const {return B.rows()==dim && B.cols()==1 && B.norm()!=0; }
bool EqCoef::hasC() const {return C!=0;}


void to_json  (Json& j, const EqCoef& eq)
{
    j["dim"]=eq.dim;
    if (eq.hasC() ) j["C"]=eq.C;
    if (eq.hasA() ) j["A"]=eq.A;
    if (eq.hasB() ) j["B"]=eq.B;
}
void from_json(const Json& j, EqCoef& p)
{
    FMatrix A(0,0),B(0,1); real_t C=0;
    try {A=j.at("A").get<FMatrix>();} catch (...) {}
    try {B=j.at("B").get<FMatrix>();} catch (...) {}
    try {C=j.at("C").get<real_t>();}  catch (...) {}
    if (A.rows()==0 && B.rows()==0 && C==0) C=1;
    p=EqCoef({j["dim"].get<int>(),A,B,C});
}



void to_json  (Json& j, const SecondOrderModel& p)
{
    j["type"]="SecondOrderModel";
    j["coefs"]=p.coefs;
    GeoPtr ptr(p.geo);
    j["geometry"]=ptr;
    ptr.release();
}
void from_json(const Json& j, SecondOrderModel& p)
{
    if (j["type"]!=std::string("SecondOrderModel"))
        throw std::logic_error("different type expected");
    p=SecondOrderModel(j["geometry"].get<EqCoef>(),j["geometry"].get<GeoPtr>().release());
}

thread_local std::vector<real_t> SecondOrderModel::m_data;

std::vector<Model::Part> SecondOrderModel::initParts(const TensorQuadrature  &quad) const
{
    if (geo !=nullptr)
        return initWithGeo(quad);
    else
        return initNoGeo(quad);
}

void SecondOrderModel::done( ) const
{
    m_data.resize(0);
}

#include <maps/bc_map.h>
std::vector<Model::Part> SecondOrderModel::initWithGeo(const TensorQuadrature  &quad) const
{
    std::vector<Part> data;
    const int dim = quad.domainDim();
    const index_t numNodes=quad.size();

    // allocate memory
    int numComp=0;
    if (coefs.hasA())
    {
        numComp+=dim*dim;
    }
    if (coefs.hasB())
    {
        numComp+=dim;
    }
    if (coefs.hasC())
    {
        ++numComp;
    }
    m_data.resize(numComp*numNodes);

    // prepare parts
    PartialDerivative tstD(0);
    PartialDerivative trlD(0);

    std::vector<view<real_t>>        tA(dim*dim);
    std::vector<view<real_t>>        tB(dim);
    view<real_t> tC;

    real_t *next=m_data.data();
    if (coefs.hasA())
    {
        for (int i=0;i<dim;++i)
        {
            tstD[i]=1;
            for (int j=0;j<dim;++j)
                {
                tA[dim*i+j]=view<real_t>(next, numNodes);
                trlD[j]=1;
                data.push_back(Part{tstD,trlD,tA[dim*i+j]} );
                next+=numNodes;
                trlD[j]=0;
                }
            tstD[i]=0;
        }
    }
    if (coefs.hasB())
    {
        for (int i=0;i<dim;++i)
                {
                tB[i]=view<real_t>(next,numNodes);
                trlD[i]=1;
                data.push_back(Part{tstD,trlD,tB[i]} );
                next+=numNodes;
                trlD[i]=0;
                }
    }
    if (coefs.hasC())
    {
        tC=view<real_t>(next,numNodes);
        data.push_back(Part{tstD,trlD,tC} );
    }

    // compute geometry and write coefficients
    GeoPtr ptr(geo);
    transformCoefs(ptr,quad.grid(),coefs,tA,tB,tC);
    ptr.release();
    return data;
}

std::vector<Model::Part> SecondOrderModel::initNoGeo  (const TensorQuadrature  &quad) const
{
    std::vector<Model::Part> data;
    const int dim = quad.domainDim();
    assert(dim==coefs.dim);
    const index_t numNodes=quad.size();
    
    // allocate memory
    int numComp=0;
    if (coefs.hasA())
    {
        for (int i=0;i<dim;++i)
        for (int j=0;j<dim;++j)
            if (coefs.A(i,j)!=0)
                ++numComp;
    }
    if (coefs.hasB())
    {
        for (int i=0;i<dim;++i)
            if (coefs.B(i)!=0)
                ++numComp;
    }
    if (coefs.hasC())
    {
        ++numComp;
    }
    m_data.resize(numComp*numNodes);
    
    // prepare parts
    PartialDerivative tstD(0);
    PartialDerivative trlD(0);
    real_t *next=m_data.data();
    if (coefs.hasA())
    {
        for (int i=0;i<dim;++i)
        {
            tstD[i]=1;
            for (int j=0;j<dim;++j)
            if (coefs.A(i,j)!=0)
                {
                trlD[j]=1;
                view<real_t> v(next, numNodes);
                v.varray()=coefs.A(i,j);
                next+=numNodes;
                data.push_back(Part{tstD,trlD,v} );
                trlD[j]=0;
                }
            tstD[i]=0;
        }
    }
    if (coefs.hasB())
    {
        for (int i=0;i<dim;++i)
            if (coefs.B(i)!=0)
                {
                trlD[i]=1;
                view<real_t> v(next, numNodes);
                v.varray()=coefs.B(i);
                next+=numNodes;
                data.push_back(Part{tstD,trlD,v} );
                trlD[i]=0;
                }
    }
    if (coefs.hasC())
    {
        view<real_t> v(next, numNodes);
        v.varray()=coefs.C;
        data.push_back(Part{tstD,trlD,v} );
    }
    return data;
}
