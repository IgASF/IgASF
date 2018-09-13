/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <bases/function_data.h>
#include <tools/registry.h>


struct Quadrature :Factory<Quadrature>
{
    virtual void  applyToValues (BasisValues &val) const =0;
    virtual const FMatrix& nodes()     const =0;
    virtual index_t size()             const =0;
    virtual int domainDim()        const =0;

    virtual ~Quadrature() = default;
};
typedef std::unique_ptr<Quadrature> QuadPtr;








