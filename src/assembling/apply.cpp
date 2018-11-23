/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <assembling/apply.h>
#include <algebra/kronecker.hpp>
#include <tools/timing.h>


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
                       )
{
    const int dim=trls.size();

    // premultiply coefficients by trial evaluation
    {
    typedef decltype(trls[0][0].transpose()) OP_TYPE;
    auto t_start = std::chrono::high_resolution_clock::now();
    
    std::vector<OP_TYPE> ops(dim, trls[0][0].transpose());
    for (int c=0;c<dim;++c)
    { new(&ops[c]) decltype(trls[0][0].transpose())(trls[dim-c-1][dtrl[dim-c-1]].transpose()); }
    
    applyRight(1,  in_v, eval_mem, tmp, view<const OP_TYPE>(ops));
    eval_mem.varray()*=coefs.varray();
    
    auto t_end = std::chrono::high_resolution_clock::now();
    time_apply_trial.fetch_add(std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    }
    
    // integrate and sum to output vector
    {
    typedef decltype(tsts[0][0].transpose().transpose()) OP_TYPE;
    auto t_start = std::chrono::high_resolution_clock::now();
    
    std::vector<OP_TYPE> ops(dim, tsts[0][0].transpose().transpose());
    for (int c=0;c<dim;++c)
    { new(&ops[c]) OP_TYPE (tsts[dim-c-1][dtst[dim-c-1]].transpose().transpose()); }
    
    applyRight(1, eval_mem, int_mem, tmp, view<const OP_TYPE>(ops));
    
    auto t_end = std::chrono::high_resolution_clock::now();
    time_apply_kronecker.fetch_add(std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    }
    out_v.vector()+=int_mem.vector();
}
