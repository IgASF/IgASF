/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <tools/registry.h>
#include <tools/grid.h>

struct GeoMap : public Factory<GeoMap>
{
    typedef std::unique_ptr<GeoMap> GeoPtr;

    virtual int targetDim() const=0;
    virtual int domainDim() const=0;
    
    virtual FMatrix evaluate(const view<real_t> &points) const=0;
    virtual FMatrix jacobian(const view<real_t> &points) const=0;
    
    virtual FMatrix evaluate(const CartesianGrid &points) const;
    virtual FMatrix jacobian(const CartesianGrid &points) const;
    
    virtual ~GeoMap()=default;
};
typedef std::unique_ptr<GeoMap> GeoPtr;
