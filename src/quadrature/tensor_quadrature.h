/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <bases/tensor_basis.h>
#include <quadrature/quadrature.h>


QuadPtr getRecommendedQuadrature(const TensorBasis& tstRef, const TensorBasis& trlRef);

struct TensorQuadrature : public SubType<Quadrature, TensorQuadrature>
{
    TensorQuadrature()=default;
    TensorQuadrature(std::vector<QuadPtr> &qu);

    virtual void applyToValues  (BasisValues &val) const;
    virtual void applyToValues  (std::vector<BasisValues> &vals) const;
    virtual const FMatrix& nodes() const;
    virtual CartesianGrid  grid()  const;
    virtual index_t        size()  const;
    virtual int domainDim()    const;
    
    QuadPtr&       operator[](int c){return m_comps[c];}
    const QuadPtr& operator[](int c) const {return m_comps[c];}

    std::vector<QuadPtr> m_comps;
    mutable FMatrix      m_nodes;
    bool                 m_compu; // are nodes computed

    friend void to_json  (Json& j, const TensorQuadrature& p);
    friend void from_json(const Json& j, TensorQuadrature& p);
};
void to_json  (Json& j, const TensorQuadrature& p);
void from_json(const Json& j, TensorQuadrature& p);
