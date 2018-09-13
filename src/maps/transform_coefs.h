/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <maps/geomap.h>
#include <assembling/second_order.h>

void transformCoefs(
    const GeoPtr& geo,
    const CartesianGrid &nodes,
    view<const real_t> A,
    view<const real_t> b,
    view<const real_t> c,
    std::vector<view<real_t>>  transformedA,
    std::vector<view<real_t>>  transformedB,
    view<real_t>               transformedC
    );


void transformCoefs(
    const GeoPtr& geo,
    const CartesianGrid &nodes,
    EqCoef eq,
    std::vector<view<real_t>> transformedA,
    std::vector<view<real_t>> transformedB,
    view<real_t>              transformedC);
