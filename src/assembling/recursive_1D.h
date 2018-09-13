/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/mmatrix.h>
#include <bases/function_data.h>
#include <tools/view.hpp>

// 1D optimized limits
constexpr index_t MAX_PT =10;
constexpr index_t MAX_TST=10;
constexpr index_t MAX_TRL=10;

// Fast implementation of bilinear sum-factorization that uses element information and provided temporaries
void assemble1D(
                const BasisValues          &mtst,
                const BasisValues          &mtrl,
                PartialDerivative           dtst,
                PartialDerivative           dtrl,
                view<const real_t>          coefs,
                view<real_t>                mOut,
                const Sparsity             &mSpr,
                view<const index_t>         eles
                );
