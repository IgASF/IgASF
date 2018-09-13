/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <bases/basis.h>

struct Bspline : public SubType<Basis,Bspline,Basis1D>
{
    Bspline()=default;
    Bspline (int d, const view<const real_t> kns);
    Bspline (int d, std::vector<real_t> kns);
    index_t  size()    const {return m_kns.size()-m_deg-1;}
    BasisValues evaluate (std::vector<PartialDerivative> der,  view< const real_t> xs ) const;
    const std::vector<real_t>& knots()  const {return m_kns;}
          std::vector<real_t>& knots()        {return m_kns;}
protected:
    std::vector<real_t> m_kns;
};
void to_json  (Json& j, const Bspline& p);
void from_json(const Json& j, Bspline& p);
// Bspline


