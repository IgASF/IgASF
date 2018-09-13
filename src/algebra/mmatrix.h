/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/sparsity.h>

struct MMatrix : public Sparsity
{
    MMatrix()=default;
    // return 0 matrix
    MMatrix(Sparsity &&sp, std::vector<real_t> &&co) // t
    : Sparsity(std::move(sp)), coefs(std::move(co))
    {}
    MMatrix(Sparsity &&sp)
    : Sparsity(std::move(sp))
    {coefs.resize(nnzs(),0);}
    std::vector<real_t> coefs;
    
    inline const real_t* data(index_t r=0) const {return coefs.data()+Sparsity::start(r);}
    inline       real_t* data(index_t r=0)       {return coefs.data()+Sparsity::start(r);}

    inline CSMap matrix() const  {return Sparsity::matrix(coefs.data());}
    inline SMap  matrix()        {return Sparsity::matrix(coefs.data());}
    
    operator CSMap() const  {return matrix();}
    operator SMap()         {return matrix();}

};
