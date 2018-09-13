/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <bases/bspline.h>
#include <tools/utils.hpp> // factorial

Bspline::Bspline(int d, std::vector<real_t> knots)
: m_kns(knots)
{
    m_deg=d;
    std::unique_copy(m_kns.begin()+d, m_kns.end()-d, std::back_inserter(m_brks));
}

Bspline::Bspline(int d, const view<const real_t> knots)
: m_kns(knots.begin(), knots.end())
{
    m_deg=d;
    std::unique_copy(m_kns.begin()+d, m_kns.end()-d, std::back_inserter(m_brks));
}

namespace {
void bsplines(index_t deg, index_t der, real_t x, typename std::vector<real_t>::const_iterator k, FMap o);
} // anonymous namespace

BasisValues Bspline::evaluate (std::vector<PartialDerivative> derivs,  view<const real_t> xs ) const
{
    BasisValues res;
    initBasisValues(res,derivs.size(), xs.size() ,size(),m_deg+1);
    res.derVec=derivs;
    
    auto k=m_kns.begin()+m_deg;
    
    size_t pos=0;
    FMap val(nullptr,m_deg+1,1);
    
    for (auto x: xs)
        {
        while (k<m_kns.end()-m_deg-1 && *(k+1)<=x ) ++k; // find segment
        
        for (int d=0; d<static_cast<int>(derivs.size()); ++d)
            {
            new (&val) FMap(res.data(d)+pos, m_deg+1, 1);
            bsplines(m_deg, derivs[d].raw, x, k, val );
            }
        for (int i=0;i<m_deg+1;++i)
            res.m_idxs[pos+i]=k-m_kns.begin()+i-m_deg;
        pos+=m_deg+1;
        }
    return res;
}

namespace {

void bsplines(index_t deg, index_t der, real_t x, typename std::vector<real_t>::const_iterator k, FMap o)
{
    if (der>deg) { o.setZero(); return; }
    o(0)=factorial(deg,deg-der);
    // first loop with derivatives coefficients so that the number
    // of possibly unstable operations is minimized
    for (int r=1;r<=der;++r)
        {
        o(r)=o(r-1)/(*(k+r)-*k);
        for(index_t j=r-1;j>0;--j)
            {
            o(j)=o(j-1)/(*(k+j)-*(k+j-r)) -o(j)/(*(k+j+1)-*(k+j+1-r));
            }
        o(0)/=-(*(k+1)-*(k+1-r));
        }
    // then complete the recursion for the non derivative
    for (int r=der+1; r<=deg;++r)
        {
        o(r)=o(r-1)*(x-*k)/(*(k+r)-*k);
        for(index_t j=r-1;j>0;--j)
            {
            o(j)=o(j-1)*(x-*(k+j-r))/(*(k+j)-*(k+j-r))+o(j)*(*(k+j+1)-x)/(*(k+j+1)-*(k+j+1-r));
            }
        o(0)*=(*(k+1)-x)/(*(k+1)-*(k+1-r));
        }
}

} // anonymous namespace

void to_json  (Json& j, const Bspline& p)
{
    j = Json({ {"type", "Bspline"}, {"degree", p.degree()}, {"knots", p.knots()} });
}

void from_json(const Json& j, Bspline& p)
{
    if ( j.at("type").get<std::string>() != "Bspline")
        throw std::logic_error("Expected Bspline, got "+j.at("type").get<std::string>() );
    p=Bspline(j.at("degree").get<int>(), j.at("knots").get<std::vector<real_t>>());
}


