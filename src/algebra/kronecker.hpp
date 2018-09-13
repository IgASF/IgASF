/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <tools/view.hpp>

struct MergeAssign
{ template <typename T1, typename T2> void operator() (T1 &dst, T2 &src) {dst.noalias()=src.eval();} };
struct MergeSum
{ template <typename T1, typename T2> void operator() (T1 &dst, T2 &src) {dst.noalias()+=src;} };
struct MergeProduct
{ template <typename T1, typename T2> void operator() (T1 &dst, T2 &src) {dst.noalias()*=src;} };

// apply operator to all blocks of proper size
template <typename OpType, typename Merge=MergeAssign>
void  applyRight(index_t block, view<const real_t> src, view<real_t> dst, const OpType &op)
{
    index_t srcCols =op.rows();
    index_t dstCols =op.cols();
    index_t blkNum  =src.size()/srcCols/block;
    
    CFMap srcMat = const_cast<const view<const real_t>&>(src).matrix(block,blkNum*srcCols);
    FMap  dstMat = dst.matrix(block,blkNum*dstCols);
    
    for (index_t b=0; b<blkNum; ++b)
        {
        auto dstBlk=dstMat.middleCols(b*dstCols,dstCols );
        auto srcBlk=srcMat.middleCols(b*srcCols,srcCols )*op;
        Merge()(dstBlk,srcBlk );
        }
}

// apply many operators one at a time as in Kronecker composition, all of the same type
// applyRight(n, src,dst, ( A , B, ..., C)) computes  src* (A x B x ... x C x In) and stores it in dst
// where In is the n-dimensional identity

template <typename OP>
void applyRight (index_t block, view<const real_t> src, view<real_t> dst, view<real_t> tmp, view<const OP> ops)
{
    view<const real_t>  msrc=src;
    if (ops.size()%2 )
        {
        applyRight( block, msrc, dst, ops.back() );
        block *=ops.back().cols();
        msrc   =dst.subView(0, src.size()/ops.back().rows()*ops.back().cols() );
        ops.popBack();
        }
    for (int i=static_cast<int>(ops.size()-1); i>0; --i)
        {
        applyRight( block, msrc, tmp, ops[i] );
        block *= ops[i].cols();
        msrc = tmp.subView(0, msrc.size()/ops[i].rows()*ops[i].cols());
        
        --i;
        
        applyRight( block, msrc, dst, ops[i] );
        block *= ops[i].cols();
        msrc = dst.subView(0, msrc.size()/ops[i].rows()*ops[i].cols());
        }
}


template <typename OP>
void applyRight (index_t block, view<const real_t> src, view<real_t> dst, view<const OP> ops)
{
    index_t max_size=0;
    {
    index_t size=src.size();
    for (int i=ops.size()-1; i>=0;--i)
        {
        size/=ops[i].rows();
        size*=ops[i].cols();
        max_size=std::max(max_size, size);
        }
    }
    
    std::vector<real_t> tmp_vec(max_size);
    view<real_t>        tmp(tmp_vec);
    
    applyRight (block, src, dst, tmp, ops);
}


template <typename OP>
void applyRight (index_t block, const view<const real_t> src,view<real_t> dst, view<real_t> , OP op)
{ applyRight(block,src,dst,op); }

inline
void applyRight1 (index_t , view<const real_t> , view<real_t> , view<real_t>)
{}

template <typename OP, typename ...OPS>
void applyRight1 (index_t block, view<const real_t> src, view<real_t> dst, view<real_t> tmp, OP op, OPS...ops)
{
    index_t newblock=block* op.cols();
    index_t newsize=src.size()*op.cols()/op.rows();
    if (sizeof...(OPS)%2) {
        applyRight ( block,    src, tmp, op);
        applyRight1( newblock, tmp.subView(0,newsize) , dst, tmp, ops...);
    } else {
        applyRight ( block,    src, dst, op);
        applyRight1( newblock, dst.subView(0,newsize) , dst, tmp, ops...);
    }
}


template <typename OP,typename ...OPS>
index_t max_size(index_t src_size, OP op, OPS...ops)
{
    index_t nsize =src_size*op.cols()/op.rows();
    if constexpr (sizeof...(OPS)==0)
        return nsize;
    else
        return std::max(nsize,max_size(nsize, ops...) );
}


// apply many operators one at a time as in Kronecker composition, possibly of different types
// due to template mechanism the operators have to be passed in REVERSED order!!!
// applyRight(1, src,dst,  A , B ) computes  src* (B x A) and stores it in dst
// this is not avoidable without too much template glue code
// (see https://stackoverflow.com/questions/15904288/how-to-reverse-the-order-of-arguments-of-a-variadic-template-function )
template <typename ...OPS>
void applyRight (index_t block, view<const real_t> src, view<real_t> dst, OPS...ops)
{
    OwningView<real_t> tmp(max_size(src.size(), ops...));
    applyRight1( block, src, dst, tmp, ops...);
}





