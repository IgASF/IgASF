/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <maps/bc_map.h>
#include <algebra/kronecker.hpp>
#include <tools/timing.h>

//REGISTER_SUBTYPE(GeoMap,BasisCoefficientMap);

namespace {
// forward declaration of implementation
FMatrix evaluateTensor(const TensorBasis& , const FMatrix & ,const CartesianGrid &points);
FMatrix jacobianTensor(const TensorBasis&, const FMatrix &, const CartesianGrid &points);
} // anonymous namespace

// public interface
BasisCoefficientMap::BasisCoefficientMap(BasisPtr &bas, const FMatrix &coefs)
: m_basis(std::move(bas)), m_coefs(coefs)
{m_tensor= dynamic_cast<TensorBasis*>(m_basis.get())!=nullptr;}

BasisCoefficientMap::BasisCoefficientMap(BasisPtr &bas,  FMatrix &&coefs)
: m_basis(std::move(bas)), m_coefs(coefs)
{m_tensor= dynamic_cast<TensorBasis*>(m_basis.get())!=nullptr;}


int BasisCoefficientMap::targetDim() const
{return m_coefs.rows();}

int BasisCoefficientMap::domainDim() const
{return m_basis->domainDim();}

FMatrix BasisCoefficientMap::evaluate(const view<real_t> &points) const
{
    BasisValues val=m_basis->evaluate(std::vector<PartialDerivative>(1,PartialDerivative(0)), points);
    return m_coefs*val[0].transpose();
}

FMatrix BasisCoefficientMap::jacobian(const view<real_t> &points) const
{
    FMatrix result( domainDim()*targetDim(),points.size()/domainDim() );
    BasisValues val=m_basis->evaluate(gradient(domainDim()), points);
    for (int d=0; d< domainDim(); ++d)
        result.block(0,d*targetDim(),result.rows(),targetDim())= m_coefs*val[d].transpose();
    return result;
}

FMatrix BasisCoefficientMap::evaluate(const CartesianGrid &points) const
{
    if (m_tensor) return evaluateTensor(*static_cast<TensorBasis*>(m_basis.get()),m_coefs,points);
    else return GeoMap::evaluate(points);
}


FMatrix BasisCoefficientMap::jacobian(const CartesianGrid &points) const
{
    if (m_tensor)
        return jacobianTensor(*static_cast<TensorBasis*>(m_basis.get()),m_coefs,points);
    else return GeoMap::evaluate(points);
}

void to_json(Json& j, const BasisCoefficientMap& p) {
    j = Json{ { "type", "BasisCoefficientMap"}, {"basis", p.m_basis}, {"coefs", p.m_coefs} };
}

void from_json(const Json& j, BasisCoefficientMap& p) {
    if ( j.at("type").get<std::string>() != "BasisCoefficientMap")
        throw std::runtime_error("expected BasisCoefficientMap , got "+j.at("type").get<std::string>() );
    p.m_basis=j.at("basis").get<BasisPtr> ();
    p.m_coefs=j.at("coefs").get<FMatrix>  ();
    p.m_tensor= dynamic_cast<TensorBasis*>(p.m_basis.get())!=nullptr;
}

namespace{
FMatrix evaluateTensor(const TensorBasis& bas, const FMatrix &coefs ,const CartesianGrid &points)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    int tar=coefs.rows();
    int dom=bas.domainDim();
    
    std::vector<std::vector<PartialDerivative>> der(dom,std::vector<PartialDerivative>(1,PartialDerivative(0)));
    auto val=bas.evaluateComponents(der, points);
    
    FMatrix result(tar, points.numPoints());
    switch (dom)
    {
    case 1:
        applyRight(tar, view<const real_t>(coefs), view<real_t>(result), val[0][0].transpose() );
    break;
    case 2:
        applyRight(tar, view<const real_t>(coefs), view<real_t>(result),
        val[0][0].transpose(),val[1][0].transpose() );
    break;
    case 3:
        applyRight(tar, view<const real_t>(coefs), view<real_t>(result), val[0][0].transpose(), val[1][0].transpose(),
        val[2][0].transpose());
    break;
    case 4:
        applyRight(tar, view<const real_t>(coefs), view<real_t>(result),
        val[0][0].transpose(), val[1][0].transpose(), val[2][0].transpose(), val[3][0].transpose());
    break;
    default: throw std::runtime_error("Implemented only for domains of dimension up to 4");
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    time_geo_compute.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    return result;
}

FMatrix jacobianTensor(const TensorBasis& bas, const FMatrix &coefs, const CartesianGrid &points)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    const int dom=bas.domainDim();
    const int tar=coefs.rows();
    const int num=points.numPoints();

    FMatrix jacs(dom*tar,num);

    std::vector<std::vector<PartialDerivative>> der(dom,std::vector<PartialDerivative>({PartialDerivative(0),PartialDerivative(1)}));
    auto val=bas.evaluateComponents(der, points);

    typedef decltype(val[0][0].transpose()) OP_TYPE;
    std::vector<OP_TYPE> ops(dom, val[0][0].transpose());
    for (int c=0;c<dom;++c)
    { new(&ops[c]) decltype(val[0][0].transpose())(val[dom-c-1][0].transpose()); }
    
    index_t max_size=0;
    {
    index_t size=coefs.size();
    for (int i=ops.size()-1; i>=0;--i)
        {
        size/=ops[i].rows();
        size*=ops[i].cols();
        max_size=std::max(max_size, size);
        }
    }
    std::vector<real_t> tmp(num*tar+max_size);
    
    view<real_t> partialDer = view<real_t>(tmp).subView(0,num*tar);
    view<real_t> workMem    = view<real_t>(tmp).subView(num*tar, std::max<index_t>(num*tar,max_size  ));
    for (int c=0;c<dom;++c)
    {
        new(&ops[dom-c-1]) OP_TYPE(val[c][1].transpose());
        applyRight(tar,view<const real_t>(coefs), partialDer,workMem, view<const OP_TYPE>(ops));
        new(&ops[dom-c-1]) OP_TYPE(val[c][0].transpose());
        jacs.middleRows(c*tar,tar)=partialDer.matrix(tar,num);
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    time_geo_compute.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    return jacs;
}

} // anonymous namespace


