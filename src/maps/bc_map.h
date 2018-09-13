/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once 

#include <maps/geomap.h>
#include <bases/tensor_basis.h>

struct BasisCoefficientMap : public SubType<GeoMap,BasisCoefficientMap>
{
    BasisCoefficientMap() = default;
    BasisCoefficientMap(BasisPtr &bas, const  FMatrix& coefs);
    BasisCoefficientMap(BasisPtr &bas,  FMatrix &&coefs);
    
    bool                     m_tensor;
    BasisPtr                 m_basis;
    FMatrix                  m_coefs; // one control point per column

    int targetDim() const;
    int domainDim() const;

    FMatrix evaluate(const view<real_t> &points) const;
    FMatrix jacobian(const view<real_t> &points) const;
    
    FMatrix evaluate(const CartesianGrid &points) const;
    FMatrix jacobian(const CartesianGrid &points) const;
};
void to_json  (Json& j, const BasisCoefficientMap& p);
void from_json(const Json& j, BasisCoefficientMap& p);
