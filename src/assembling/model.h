/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/mmatrix.h>
#include <quadrature/tensor_quadrature.h>
#include <bases/tensor_basis.h>
#include <tools/registry.h>


struct Model : Factory<Model>
{
public:
    struct Part
    {
    PartialDerivative  test;
    PartialDerivative  trial;
    view<const real_t> coefs;
    };

    std::vector<Part> parts;
public:
    void assemble(const TensorBasis &test, const TensorBasis &trial, const TensorQuadrature  &quad, MMatrix &output) const;
protected:
    virtual std::vector<Part> initParts(const TensorQuadrature  &) const {return parts;}
    virtual void done     ( ) const {}

protected:
    static std::vector<Sparsity> getBilSparsities(const std::vector<BasisValues> &tsts, const std::vector<BasisValues> &trls );
    static std::vector<Sparsity> getKroSparsities(const std::vector<Sparsity> &sprs);
    std::vector<std::vector<index_t>> computeElementSplitting(const TensorBasis &test, const TensorBasis &trial, const TensorQuadrature  &quad)const;
    typedef const PartialDerivative (Part::*Role);
    std::vector<PartialDerivative> getModelRequest(Role role, int component,const std::vector<Part> &data )const;
    TensorBasis::ComponentRequest getTestRequest (int dim, const std::vector<Part> &data)const;
    TensorBasis::ComponentRequest getTrialRequest(int dim, const std::vector<Part> &data)const;
};
typedef std::unique_ptr<Model> ModelPtr;

void to_json  (Json& j, const Model::Part& p);



