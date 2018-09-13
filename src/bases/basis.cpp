/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <bases/basis.h>
#include <bases/bspline.h>
#include <bases/tensor_basis.h>

REGISTER_SUBTYPE(Basis,Bspline);
REGISTER_SUBTYPE(Basis,TensorBasis);


// force instantiaton of Factory<Basis>
