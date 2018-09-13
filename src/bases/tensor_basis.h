/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <bases/basis.h>
#include <tools/grid.h>

struct TensorBasis : public SubType<Basis,TensorBasis>
{
    typedef std::vector<std::vector<PartialDerivative>> ComponentRequest;

    TensorBasis() = default;
    TensorBasis(std::vector<BasisPtr> & );
    TensorBasis(std::vector<BasisPtr> && );
    
    index_t size() const;
    const Basis1D& operator[](int b) const {return static_cast<const Basis1D&>(*m_bas[b]);}
          Basis1D& operator[](int b)       {return static_cast<Basis1D&>(*m_bas[b]);}
    virtual BasisValues evaluate (std::vector<PartialDerivative> derivs, view<const real_t> pts ) const;
    virtual std::vector<BasisValues> evaluateComponents (const ComponentRequest derivs, const CartesianGrid &xs ) const;
protected:
    std::vector<BasisPtr> m_bas;
    friend void to_json  (Json& j, const TensorBasis& p);
    friend void from_json(const Json& j, TensorBasis& p);
};
void to_json  (Json& j, const TensorBasis& p);
void from_json(const Json& j, TensorBasis& p);


// conversion to  component request from global request
std::vector<std::vector<PartialDerivative> > compDeriv (const std::vector<PartialDerivative> &der, int mdim, IMatrix &rec);

