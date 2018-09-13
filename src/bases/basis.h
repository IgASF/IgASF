/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <tools/registry.h>
#include <tools/view.hpp>
#include <bases/function_data.h>
#include <quadrature/quadrature.h>

// Interface for a real valued generating system
// this has a dynamic hieararchy and should be used through
// shared pointers

struct Basis : public Factory<Basis>
{
    typedef std::unique_ptr<Basis> BasisPtr;
public:
    int domainDim() const {return m_domDim;}
    
    virtual index_t size()    const = 0;
    virtual BasisValues evaluate (std::vector<PartialDerivative> derivs, view<const real_t> xs ) const = 0;
    virtual ~Basis() = default;
protected:
    Basis(int dom) : m_domDim(dom) {}
    Basis() =default;
    int m_domDim;
};
typedef std::unique_ptr<Basis> BasisPtr;


struct Basis1D : public Basis
{
    const std::vector<real_t>& breaks() const {return m_brks;};
    int                        degree() const {return m_deg;}
protected:
    Basis1D() : Basis(1) {}
// this is used to construct default quadrature
    std::vector<real_t> m_brks;
    int                 m_deg;
};
