/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/mmatrix.h>
#include <bases/function_data.h>
#include <tools/view.hpp>

// coefs gets overwritten by the product of them with the evaluated trial
void KroneckerApply(
                       view<const BasisValues>         tsts,
                       view<const BasisValues>         trls,
                       PartialDerivative               dtst,
                       PartialDerivative               dtrl,
                       view<const real_t>              coefs,
                       view<real_t>                    eval_mem,
                       view<real_t>                    int_mem,
                       view<real_t>                    tmp,
                       view<const real_t>              in_v,
                       view<real_t>                    out_v
                       );

