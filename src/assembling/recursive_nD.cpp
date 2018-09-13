/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <assembling/recursive_nD.h>
#include <assembling/recursive_1D.h>
#include <tools/utils.hpp>

namespace {

template <int Pts>
auto SumExpr(
             index_t block, index_t offset,
             index_t TS, index_t tstJ, const real_t *tstV,
             index_t TR, index_t trlJ, const real_t *trlV,
             rview<2,real_t> mMem)
{
    if constexpr(Pts==1) {
        return mMem[0].subView(offset,block).vector()*tstV[TS]*trlV[TR];
    } else {
        return mMem[0].subView(offset,block).vector()*tstV[TS]*trlV[TR]+SumExpr<Pts-1>(block,offset,TS+tstJ,tstJ,tstV,TR+trlJ,trlJ,trlV, mMem.dropFront());
    }
}


template <int pt_blk>
struct writePtBlock
{
    void operator()(
               const BasisValues& mTst, int derTst,
               const BasisValues& mTrl, int derTrl,
               index_t pt,
               const Sparsity  &mSpr,
               const Sparsity  &mKrn,
               const Sparsity  &tKrn,
               rview<2,real_t>  mMem,
               view<real_t>     mOut )
    {
    const real_t  *tstV=mTst.data(derTst)+mTst.start(pt);
    const auto     tstI=mTst.rowPos(pt);
    index_t  tstN=mTst.nnzs(pt);
    
    const real_t  *trlV=mTrl.data(derTrl)+mTrl.start(pt);
    const auto     trlI=mTrl.rowPos(pt);
    index_t  trlN=mTrl.nnzs(pt);
    
    for (index_t TS=0; TS<tstN;++TS)
        {
        const auto pos  = mSpr.rowPos(tstI[TS]);
        for (index_t r=0; r<tKrn.rows(); ++r )
            {
            const auto block=tKrn.nnzs(r);
            const auto curTst=tstI[TS]*tKrn.rows()+r;
            
            auto iter = pos.begin();
            for (index_t TR=0; TR<trlN;++TR)
                {
                iter=std::find(iter,pos.end(), trlI[TR]);
                real_t *curOut= mOut.begin()+mKrn.start(curTst)+block*(iter-pos.begin());
                FMap dst(curOut, block, 1);
                dst+=SumExpr<pt_blk>(block,tKrn.start(r),TS,tstN,tstV,TR,trlN,trlV, mMem);
                }
            }
        }
    }
};

} // anonymous namespace

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
                       )
{
    auto dim = tsts.size();
    if (dim==1)
        {
        assemble1D(tsts.front(),trls.front(), dtst, dtrl, coefs,  mOut, sprs.front(), eles.front());
        return;
        }
    
    // data for the current dimension
    const auto &mTst = tsts.back();
    const auto &mTrl = trls.back();
    const auto &mSpr = sprs.back();
    const auto &mKrn = krns.back();
    const auto &mEle = eles.back();
    auto mMem = mems.back();

    // arguments of the recursive call
    auto rTsts = tsts.dropBack();
    auto rTrls = trls.dropBack();
    auto rMems = mems.dropBack();
    auto rSprs = sprs.dropBack();
    auto rKrns = krns.dropBack();
    auto rEles = eles.dropBack();
    // result of the recursive call
    const auto &tKrn = krns.dropBack().back();
    
    // useful constants
    const auto pt_num    = mTst.rows();
    const auto coefBlock = coefs.size()/pt_num;
    const auto derTst=mTst.getIndex(dtst[dim-1]);
    const auto derTrl=mTrl.getIndex(dtrl[dim-1]);
    
    
    const index_t *pt_end=mEle.begin();
    index_t pt_beg=0;
    while (pt_beg<pt_num)
    {
        const index_t pt_blk=min(MAX_TMP,static_cast<index_t>(mMem.size()),*pt_end-pt_beg);
        // compute lower dimensional integrals
        for (int t=0;t<pt_blk;++t)
            {
            mMem[t].subView(0,tKrn.nnzs()).vector().setZero();
            recursiveAssemble(rTsts, rTrls, dtst, dtrl, coefs.subView((pt_beg+t)*coefBlock,coefBlock),mMem[t], rSprs,rKrns, rEles, rMems );
            }
        dispatch<MAX_TMP, writePtBlock >(pt_blk, std::ref(mTst), derTst, std::ref(mTrl), derTrl, pt_beg, std::ref(mSpr), std::ref(mKrn), std::ref(tKrn), mMem, mOut );
        pt_beg+=pt_blk;
        if (pt_beg==*pt_end) ++pt_end;
    }
}
