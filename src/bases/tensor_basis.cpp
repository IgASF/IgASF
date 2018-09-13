/** This file is part of the IgASF library.
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <bases/tensor_basis.h>
#include <set>
#include <cassert>

namespace  {
    CartesianGrid ptsToGrid(view<const real_t> pts, int mdim, IMatrix &rec);
    BasisValues buildTensorResult(const std::vector<BasisValues> &cs, IMatrix &der_rec, IMatrix pts_rec);
}

TensorBasis::TensorBasis(std::vector<BasisPtr> &&vec )
: m_bas(std::move(vec))
{
    m_domDim=m_bas.size();
    for(auto &b : m_bas )
        {
        if (b->domainDim()!=1) throw std::logic_error("Not Implemented!");
        }
}

TensorBasis::TensorBasis(std::vector<BasisPtr> &vec )
{
    m_domDim=vec.size();
    for(auto &b : vec )
        {
        if (b->domainDim()!=1) throw std::logic_error("Not Implemented!");
        m_bas.push_back(std::move(b));
        }
}


index_t TensorBasis::size() const
{
    index_t r=1;
    for (const auto &a : m_bas)
        {r*=a->size();}
    return r;
}

// the idea of the following function is to translate a request of
// Partial derivatives to a ComponentREQUEST SO THAT Evaluation can
// exploit the tensor product structure
std::vector<std::vector<PartialDerivative> > compDeriv (const std::vector<PartialDerivative> &der, int mdim, IMatrix &rec)
{
    std::vector<std::vector<PartialDerivative> > res(mdim);
    std::vector<std::set<PartialDerivative> >    tmp(mdim);
    
    rec.resize(mdim, der.size());
    
    for (size_t j=0;j<der.size();++j)
        {
        for (int i=0;i<mdim;++i)
            {
            tmp[i].insert(der[j][i]);
            rec(i,j)=static_cast<int>(der[j][i]);
            }
        }
    for (int i=0;i<mdim;++i)
        res[i]=std::vector<PartialDerivative>(tmp[i].begin(), tmp[i].end());
    return res;
}

BasisValues TensorBasis::evaluate (std::vector<PartialDerivative> derivs, view<const real_t> pts ) const
{
    const int mdim=domainDim();
    IMatrix der_rec;
    IMatrix pts_rec;
    auto comp=evaluateComponents(compDeriv(derivs,mdim,der_rec) ,ptsToGrid(pts,mdim,pts_rec) );
    return buildTensorResult(comp, der_rec,pts_rec);
}


std::vector<BasisValues> TensorBasis::evaluateComponents (std::vector<std::vector<PartialDerivative>> derivs, const CartesianGrid &xs ) const
{
    const int mdim=domainDim();
    assert( xs.domainDim()==mdim && derivs.size()==static_cast<size_t>(mdim));// ,"Number of components differ from dimension of grid" );
    std::vector<BasisValues> res(mdim);
    for (int c=0;c<mdim;++c)
        res[c]=m_bas[c]->evaluate(derivs[c], xs[c] );
    return res;
}


void to_json  (Json& j, const TensorBasis& p)
{
    j["type"]="TensorBasis";
    j["components"]=p.m_bas;
}

void from_json(const Json& j, TensorBasis& p)
{
    if ( j.at("type").get<std::string>() != "TensorBasis")
        throw std::logic_error("expected TensorBasis, got "+j.at("type").get<std::string>() );
    p=TensorBasis(j["components"].get<std::vector<BasisPtr>>());
}



namespace {
    std::vector<real_t> sortUniqueWithIndex (  view<const real_t> pts, int dim, int comp, decltype(IMatrix().row(0)) index2 )
    {
    index_t comp_size=pts.size()/dim;
    view<const real_t> cs = pts.subView(comp_size*comp, comp_size*(comp+1));
    
    std::vector<index_t> index1(comp_size);
    for (index_t i=0;i<comp_size;++i)
        index1[i]=i;
    
    std::sort(index1.begin(), index1.end(), [&cs](const index_t &a, const index_t &b){return cs[a]<cs[b];} );
    
    std::vector<real_t>       data;
    data.push_back(pts[0]);
    for (index_t i=1; i<comp_size;++i)
        {
        if( cs[index1[i]]> cs[index1[i-1]] )
            {
            index2(0,index1[i])=data.size();
            data.push_back(cs[index1[i]]);
            }
        else if (cs[index1[i]] == cs[index1[i-1]])
            {
            index2(0,index1[i])=index2(0,index1[i-1]);
            }
        else throw std::logic_error("Logic error in algorithm");
        }
    std::sort(data.begin(), data.end());
    return data;
    }
    
    CartesianGrid ptsToGrid(view<const real_t> pts, int mdim, IMatrix &rec)
    {
    index_t comp_size=pts.size()/mdim;
    
    FMatrix mpts=pts.matrix(mdim, pts.size()/mdim).transpose();
    std::vector<view<const real_t>> pt_cms(mdim);
    
    std::vector<std::vector<real_t>> unique_coords(mdim);
    for (int i=0;i<mdim;++i)
        {
        unique_coords[i]=sortUniqueWithIndex (pts.subView(comp_size*i,comp_size*(i+1)),mdim, i, rec.row(i));
        pt_cms[i]=view<real_t>(unique_coords[i]);
        }
    return CartesianGrid(view<view<const real_t>>(pt_cms) );
    }
    
    BasisValues buildTensorResult(const std::vector<BasisValues> &cs, IMatrix &der_rec, IMatrix pts_rec)
    {
    (void) cs; (void) pts_rec;
    BasisValues res;
    for (index_t d=0; d<der_rec.cols(); ++d)
        for (index_t pt=0; pt<pts_rec.cols(); ++pt)
            {
            
            }
    return res;
    }
    
} // anonymous namespace

