/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <assembling/recursive_1D.h>
#include <tools/utils.hpp>

namespace {

template<int row, int col>
struct accumulateT
{
    void operator()( index_t &pt, index_t el, FMatrix &tmp, int derTst, int derTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs);
};

void accumulate(index_t &pt, index_t el, FMatrix &tmp, int derTst, int derTrl, index_t numTst, index_t numTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs);

} // anonymous namespace


void assemble1D(
                const BasisValues          &tst,
                const BasisValues          &trl,
                PartialDerivative           dtst,
                PartialDerivative           dtrl,
                view<const real_t>          coefs,
                view<real_t>                mOut,
                const Sparsity             &mSpr,
                view<const index_t>         eles
                )
{
    index_t pt_beg=0; // next quad point
    static thread_local FMatrix tmp;
    for (auto pt_end:eles)
        {
        view<const index_t> activeTst = tst.rowPos(pt_beg);
        view<const index_t> activeTrl = trl.rowPos(pt_beg);
        const index_t numTst=activeTst.size();
        const index_t numTrl=activeTrl.size();
        
        tmp.resize(numTst, numTrl);
        tmp.setZero();
        
        const int derTst=tst.getIndex(dtst[0]);
        const int derTrl=trl.getIndex(dtrl[0]);
        
        if (numTst<=MAX_TST && numTrl<=MAX_TRL)
            dispatch2<MAX_TST,MAX_TRL, accumulateT>(numTst,numTrl, pt_beg,pt_end,std::ref(tmp),derTst, derTrl, std::ref(tst), std::ref(trl), coefs);
        else
            accumulate(pt_beg,pt_end,tmp,derTst, derTrl,numTst,numTrl, tst, trl, coefs);
        pt_beg=pt_end;
        
        for (index_t r=0; r<numTst ;++r)
            {
            auto idTst = activeTst[r];
            view<const index_t> idTrl=mSpr.rowPos(idTst);
            view<real_t>        data=mOut.subView(mSpr.start(idTst),mSpr.nnzs(idTst));
            for (index_t c=0; c<numTrl;++c)
                {
                while (idTrl.front()<activeTrl[c] )
                    { idTrl.popFront(1); data.popFront(1); }
                data.front()+=tmp(r,c);
                }
            }
        }
}


// addition without optimization based on numer of actives

namespace {

template <int num>
auto addExpr(index_t pt, int derTst, int derTrl, index_t numTst, index_t numTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs)
{
    auto vTst=view<const real_t>(tst.data(derTst,pt),numTst).vector();
    auto vTrl=view<const real_t>(trl.data(derTrl,pt),numTrl).vector();
    if constexpr (num==1)
        return vTst*vTrl.transpose()*coefs[pt];
    else
        return vTst*vTrl.transpose()*coefs[pt]+addExpr<num-1>(pt+1, derTst, derTrl,numTst,numTrl, tst, trl, coefs);
}

template <int num>
struct addMany
{
    void operator ()(FMatrix &dst, index_t &pt, int derTst, int derTrl, index_t numTst, index_t numTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs)
    {
    dst+=addExpr<num>(pt,derTst, derTrl,numTst,numTrl, tst, trl, coefs);
    }
};

void accumulate(index_t &pt, index_t el, FMatrix &tmp, int derTst, int derTrl, index_t numTst, index_t numTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs)
{
    while (el-pt>= MAX_PT)
    {
        addMany<MAX_PT>()(tmp, pt,derTst, derTrl,numTst,numTrl, tst, trl, coefs);
        pt+=MAX_PT;
    }
    dispatch<MAX_PT-1,addMany>(el-pt,std::ref(tmp), std::ref(pt),derTst, derTrl,numTst,numTrl, std::ref(tst), std::ref(trl), coefs);
}


// addition with optimization based on active tests and trials


template <int row, int col,int num>
struct addManyT
{
    void operator()(FMatrix &dst, index_t pt, int derTst, int derTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs)
    {
    dst+=buildExpression(pt,derTst, derTrl, tst, trl, coefs);
    }
    static auto buildExpression( index_t pt, int derTst, int derTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs)
    {
    auto vTst=view<const real_t>(tst.data(derTst,pt),row).matrix<row,1>();
    auto vTrl=view<const real_t>(trl.data(derTrl,pt),col).matrix<1,col>();
    if constexpr (num==1)
        return vTst*vTrl*coefs[pt];
    else
        return vTst*vTrl*coefs[pt]+addManyT<row,col,num-1>::buildExpression(pt+1, derTst, derTrl, tst, trl, coefs);
    }
};

template <unsigned int R, unsigned int C, template <int, int, int > typename FF>
struct accumulateTWrap
{
    template <int M>
    using F=FF<R,C,M>;
};


template<int row, int col>
    void accumulateT<row,col>::operator()( index_t &pt, index_t el, FMatrix &tmp, int derTst, int derTrl, const BasisValues &tst, const BasisValues &trl, const view<const real_t> &coefs)
    {
    while (el-pt>= MAX_PT)
        {
        addManyT<row,col,MAX_PT>()(tmp,pt,derTst, derTrl, tst, trl, coefs);
        pt+=MAX_PT;
        }
    dispatch<MAX_PT-1, accumulateTWrap<row,col,addManyT>::template F >(el-pt,std::ref(tmp),pt,derTst, derTrl, std::ref(tst), std::ref(trl), coefs);
    pt=el;
    }
    
} // anonymous namespace
