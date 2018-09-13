/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <maps/geomap.h>

struct RationalMap : public SubType<GeoMap,RationalMap>
{
    RationalMap() = default;
    RationalMap(GeoPtr & geo);

    int targetDim() const;
    int domainDim() const;

    virtual FMatrix evaluate(const view<real_t> &points) const;
    virtual FMatrix jacobian(const view<real_t> &points) const;
    
    virtual FMatrix evaluate(const CartesianGrid &points) const;
    virtual FMatrix jacobian(const CartesianGrid &points) const;

    const GeoPtr&  underlying() const;
protected:
    GeoPtr  m_original;
};
void to_json  (Json& j, const RationalMap& p);
void from_json(const Json& j, RationalMap& p);
