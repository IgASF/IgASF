/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <tools/view.hpp>

struct CartesianGrid
{
    CartesianGrid( view<view<const real_t>> points);
    view<real_t>&       operator[](int i);
    view<const real_t>  operator[](int i) const;

    int     domainDim() const;
    index_t numPoints() const;
    FMatrix toPoints()  const;

    std::vector<real_t>       m_data;
    std::vector<view<real_t>> m_comp;
};
void to_json  (Json& j, const CartesianGrid& p);
void from_json(const Json& j, CartesianGrid& p);
