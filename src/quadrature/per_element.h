/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <quadrature/quadrature.h>
#include <quadrature/quad_templates.h>

// WARNING, think well about what you want to do before putting nodes
// in the breakpoints

struct ElementQuadrature : public SubType<Quadrature,ElementQuadrature>
{
    std::vector<real_t>  m_ws; // weights
    FMatrix  m_ns; // nodes
   
    QuadTemplate         m_qu; // template
    std::vector<real_t>  m_brks;  // elements;
    
    ElementQuadrature()=default;
    virtual ~ElementQuadrature()=default;
    ElementQuadrature(QuadTemplate qu, std::vector<real_t> els);
    ElementQuadrature(view<const real_t> qs, view<const real_t> ws);

    virtual index_t size() const;
    virtual int  domainDim()    const;

    virtual void applyToValues(BasisValues &val) const;
public:
    const std::vector<real_t>& weights() const;
    const FMatrix& nodes()   const;
    friend void to_json  (Json& j, const ElementQuadrature& p);
    friend void from_json(const Json& j, ElementQuadrature& p);
};
void to_json  (Json& j, const ElementQuadrature& p);
void from_json(const Json& j, ElementQuadrature& p);

std::vector<real_t> knotsToElements (const std::vector<real_t> &kns1, const std::vector<real_t> &kns2);
