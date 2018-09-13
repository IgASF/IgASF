/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/mmatrix.h>
#include <bases/function_data.h>
#include <tools/view.hpp>

// nD optimized limits
constexpr index_t MAX_TMP=10;

void recursiveAssemble(
                       view<const BasisValues>         tsts,
                       view<const BasisValues>         trls,
                       PartialDerivative               dtst,
                       PartialDerivative               dtrl,
                       view<const real_t>              coefs,
                       view<real_t>                    mOut,
                       view<const Sparsity>            sprs,
                       view<const Sparsity>            krns,
                       rview<2,const index_t>          eles,
                       rview<3,real_t>                 mems
                       );

