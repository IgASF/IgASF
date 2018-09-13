/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <algebra/sparsity.h>

view<const index_t> Sparsity::rowPos  (index_t r) const      {return view<const index_t>(m_idxs.begin()+m_start[r],nnzs(r));}
view<index_t>       Sparsity::rowPos  (index_t r)            {return view<index_t>(m_idxs.begin()+m_start[r],nnzs(r));}

// infos
index_t Sparsity::cols() const {return m_cols;}
index_t Sparsity::rows() const {return m_start.size()-1;}
index_t Sparsity::nnzs() const {return m_idxs.size();}
// row info
index_t Sparsity::nnzs (index_t r) const {return m_start[r+1]-m_start[r];}
index_t Sparsity::start(index_t r) const {return m_start[r];}

// get full matrices
SMap     Sparsity::matrix (real_t* data)
{
    return SMap(m_start.size()-1,cols(),m_idxs.size(), m_start.begin(), m_idxs.begin(),data);
}
void     Sparsity::matrix (real_t* data, SMap&dst)
{
    new(&dst) SMap(m_start.size()-1,cols(),m_idxs.size(), m_start.begin(), m_idxs.begin(),data);
}


CSMap     Sparsity::matrix (const real_t *data) const
{
    return CSMap(m_start.size()-1,cols(),m_idxs.size(), m_start.begin(), m_idxs.begin(),data);
}

void   Sparsity::matrix (const real_t* data, CSMap&dst) const
{
    new(&dst) CSMap(m_start.size()-1,cols(),m_idxs.size(), m_start.begin(), m_idxs.begin(),data);
}

std::array<index_t, 2> Sparsity::posToRC (index_t pos) const
{
    assert(static_cast<size_t>(pos)<m_idxs.size());
    index_t r= std::upper_bound(m_start.begin(),m_start.end(), pos)-m_start.begin()-1;
    return {{r, m_idxs[pos]}};
}
index_t   Sparsity::RCToPos (index_t row, index_t col) const
{
    const auto pos=rowPos(row);
    auto cIt=std::find(pos.begin(),pos.end(),col);
    if (cIt!=pos.end())
        return cIt-m_idxs.begin();
    else return -1;
}
IMatrix Sparsity::posToRCMap () const
{
    IMatrix result(2,m_idxs.size());
    for (index_t i=0; i<result.cols();++i)
        {
        auto rc = posToRC (i);
        result(0,i)=rc[0];
        result(1,i)=rc[1];
        }
    return result;
}


Sparsity bilinearSparsity (const Sparsity &A, const Sparsity &B)
{
    Sparsity result(A.cols(),B.cols());
    result.m_start[0]=0;
    
    index_t AbegRow=0; // first row that intersect Acol in a non-zero
    index_t AendRow=0; // one past the last row with the same property

    std::vector<index_t> tmp;
    tmp.reserve(A.nnzs()+B.nnzs());
    
    for (index_t Acol=0; Acol< A.cols();++Acol) // rows of the result
    {
        while(AbegRow<A.rows() && A.rowPos(AbegRow).back()<Acol) ++AbegRow;
        AendRow=AbegRow;
        while(AendRow<A.rows() && A.rowPos(AendRow).front()<=Acol) ++AendRow;

        if (AbegRow==AendRow) continue;
    
        index_t BcolBeg = B.rowPos(AbegRow).front();
        index_t BcolEnd = B.rowPos(AendRow-1).back();
        for (;BcolBeg<=BcolEnd;++BcolBeg)
            tmp.emplace_back(BcolBeg);
        result.m_start[Acol+1]=static_cast<index_t>(tmp.size());
    }
    result.m_idxs.swap(OwningView<index_t>(tmp.begin(),tmp.end())); // result is acquiring ownership
    return result;
}


Sparsity kroneckerSparsity(const Sparsity &eOp, const Sparsity &iOp)
{
    Sparsity result(eOp.rows()*iOp.rows(),eOp.cols()*iOp.cols(), eOp.nnzs()*iOp.nnzs());
    
    result.m_start[0]=0;
    index_t curRow=1;
    index_t shift=0;
    for (index_t eR=0; eR<eOp.rows(); ++eR )     // row of the outer matrix
        {
        const view<const index_t> ePos=eOp.rowPos(eR);
        for (index_t iR=0; iR<iOp.rows(); ++iR ) // row of the inner matrix
            {
            const view<const index_t> iPos=iOp.rowPos(iR);
            for (index_t eC : ePos)              // non zero column of the outer matrix
                {
                result.m_idxs.subView(shift,iOp.nnzs(iR)).varray()=iPos.varray()+iOp.cols()*eC;
                shift += iOp.nnzs(iR);
                }
            result.m_start[curRow]=shift;
            ++curRow;
            }
        }
    return result;
}
