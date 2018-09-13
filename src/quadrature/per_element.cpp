/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <quadrature/per_element.h>
#include <tools/utils.hpp>

REGISTER_SUBTYPE(Quadrature,ElementQuadrature);

int  ElementQuadrature::domainDim()    const
{
    return m_ns.rows();
}

index_t ElementQuadrature::size() const
{ return m_ns.size(); }

const std::vector<real_t> &ElementQuadrature::weights() const
{ return m_ws; }

const FMatrix& ElementQuadrature::nodes() const
{ return m_ns; }

ElementQuadrature::ElementQuadrature(QuadTemplate qu, std::vector<real_t> els)
: m_qu(qu), m_brks(els)
 {
   std::vector<real_t> &nodes=qu.nodes;
   std::vector<real_t> &weights=qu.weights;
    // NOTE assumes col major maps
    // NOTE assumes quadrature on [-1,1]
    assert (nodes.size()==weights.size() );
    CFMap srcN(nodes.data(),nodes.size(),1);
    CFMap srcW(weights.data(),weights.size(),1);

    m_ns.resize(1,(els.size()-1)*nodes.size());
    m_ws.resize((els.size()-1)*nodes.size());
    FMap dstW(m_ws.data(),nodes.size(),els.size()-1);
    FMap dstN(m_ns.data(),nodes.size(),els.size()-1);
    index_t c=0;
    for_each_pair(els.begin(), els.end(), [srcN,&dstN,srcW,&dstW,&c](real_t a, real_t b)
    {
        dstN.col(c).array()=srcN.array()*(b-a)/2+(a+b)/2;
        dstW.col(c).array()=srcW.array()*(b-a)/2;
        ++c;
    }
    );
}

ElementQuadrature::ElementQuadrature(view<const real_t> qs, view<const real_t> ws)
    :  m_ws(ws.begin(),ws.end()), m_ns(qs.vector().transpose())
{}



void ElementQuadrature::applyToValues  (BasisValues &val) const
{
    for (index_t s=0; s< val.numData();++s)
    {
        if (static_cast<size_t>(val[s].rows())==static_cast<size_t>(m_ns.size()))
        {
            for (index_t r=0;r<val[s].rows();++r)
            {
                index_t start=val.start(r);
                index_t sz=val.nnzs(r);
                FMap block(val.data(s)+start, sz,1 );
                block*= m_ws[r];
            }
        }
    }
}

void to_json  (Json& j, const ElementQuadrature& p)
{
    j={
        {"type", "ElementQuadrature"},
        {"template", p.m_qu},
        {"elements", p.m_brks}
    };
}

void from_json(const Json& j, ElementQuadrature& p)
{
    if (j["type"] !="ElementQuadrature")
        throw std::runtime_error("expected ElementQuadrature, got "+j.at("type").get<std::string>() );
    p=ElementQuadrature(j.at("template").get<QuadTemplate>(), j.at("elements").get<std::vector<real_t>>());
}

std::vector<real_t> knotsToElements (const std::vector<real_t> &kns1, const std::vector<real_t> &kns2)
{
    std::vector<real_t> merged;
    merged.reserve(kns1.size()+kns2.size()-2);
    auto i1=kns1.begin(), i2=kns2.begin(),  e1=kns1.end(), e2=kns2.end();
    assert(*i1==*i2);
    merged.push_back(*i1);
    while (1)
    {
        const real_t a=merged.back();
        i1=std::find_if(i1,e1,[a](real_t c){return c>a;});
        i2=std::find_if(i2,e2,[a](real_t c){return c>a;});
        if (i1<e1 && i2<e2)
            merged.push_back(std::min(*i1,*i2));
        else if (i1<e1)
            merged.push_back(*i1);
        else if (i2<e2)
            merged.push_back(*i2);
        else break;
    }
    return merged;
}

