/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <quadrature/tensor_quadrature.h>
#include <quadrature/per_element.h>
#include <algebra/mmatrix.h>

REGISTER_SUBTYPE(Quadrature,TensorQuadrature);

QuadPtr getRecommendedQuadrature(const TensorBasis& tstRef, const TensorBasis& trlRef)
{
    int domDim = tstRef.domainDim();
    assert( domDim==trlRef.domainDim() );
    
    std::vector<QuadPtr>  quads;
    for (int c=0; c<domDim;++c)
        {
        int quad_deg= (tstRef[c].degree()+trlRef[c].degree())/2 +1;
        quads.push_back(
                        QuadPtr(new ElementQuadrature(Rules::Gauss(quad_deg), knotsToElements(tstRef[c].breaks(),trlRef[c].breaks()))));
        }
    return QuadPtr(new TensorQuadrature(quads));
}


int TensorQuadrature::domainDim() const
{
    return m_comps.size();
}

TensorQuadrature::TensorQuadrature(std::vector<QuadPtr> &qu)
{
    for (auto &c:qu)
        m_comps.push_back(std::move(c));
}

void TensorQuadrature::applyToValues  (BasisValues &val) const
{
    (void) val;
    throw std::runtime_error("Not implemented");
}

void TensorQuadrature::applyToValues  (std::vector<BasisValues> &vals) const
{
    for (size_t i=0;i<vals.size();++i)
        m_comps[i]->applyToValues(vals[i]);
}

const FMatrix& TensorQuadrature::nodes()     const
{
    if (! m_compu)
        m_nodes=grid().toPoints();
    return m_nodes;
}

CartesianGrid  TensorQuadrature::grid()     const
{
    std::vector<view<const real_t> > pieces(m_comps.size());
    for (size_t i=0;i<m_comps.size();++i)
        pieces[i]=view<const real_t>(m_comps[i]->nodes());
    return CartesianGrid(pieces);
}

index_t TensorQuadrature::size() const
{
    index_t r=1;
    for (auto &v : m_comps)
        r*=v->size();
    return r;
}

void to_json  (Json& j, const TensorQuadrature& p)
{
    j["type"]="TensorQuadrature";
    j["components"]=p.m_comps;
}

void from_json(const Json& j, TensorQuadrature& p)
{
    if ( j.at("type").get<std::string>() != "TensorQuadrature")
        throw std::runtime_error("expected TensorQuadrature, got "+j.at("type").get<std::string>() );
    p.m_comps=j["components"].get<std::vector<QuadPtr>>();
}

