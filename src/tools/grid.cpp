/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <tools/grid.h>

CartesianGrid::CartesianGrid( view<view<const real_t>> points)
{
    index_t size=0;
    for (size_t i=0;i<points.size();++i)
        size+=points[i].size();
    m_data.resize(size);
    size=0;
    view<real_t> vdata(m_data);
    for(size_t i=0;i<points.size();++i)
    {
        const index_t this_size=points[i].size();
        m_comp.push_back(vdata.subView(size,this_size));
        m_comp.back().matrix(this_size,1)=points[i].matrix(this_size,1);
        size+=this_size;
    }
}


view<real_t>&       CartesianGrid::operator[](int i)       {return m_comp[i];}
view<const real_t>  CartesianGrid::operator[](int i) const {return m_comp[i];}
int     CartesianGrid::domainDim() const {return m_comp.size();}
    
    
index_t CartesianGrid::numPoints() const
{
    index_t r=1;
    for (auto &v : m_comp)
        r*=v.size();
    return r;
}
 

FMatrix CartesianGrid::toPoints()  const
{
    index_t nPts=numPoints();
    FMatrix res(domainDim(),nPts);
    FMatrix pt(domainDim(),1);
    int cmp=0;
    index_t rep=1;
    index_t nBlk=nPts/m_comp[0].size();
    for (auto &c: m_comp)
    {
        const index_t cs=c.size();
        for (index_t b=0; b< nBlk;++b)
        for (index_t p=0;p<cs;++p)
            res.block(cmp,cs*rep*b+rep*p,1,rep).array()=c[p];
        rep  *=cs;
        nBlk /=cs;
        ++cmp;
    }
    return res;
}

void to_json  (Json& j, const CartesianGrid& p)
{
    j["type"]="Grid";
    for (int i=0; i< p.domainDim(); ++i)
        j["component"][i]=p[i];
}

void from_json(const Json& j, CartesianGrid& p)
{
    std::vector<std::vector<real_t> > data;
    if ( j.at("type").get<std::string>()!="Grid" )
        throw std::logic_error ("Expected Grid.");
    data = j["component"].get<std::vector<std::vector<real_t> >>();
    std::vector<view<const real_t>> vs(data.size());
    for (size_t i=0;i<data.size();++i)
    {vs[i]=view<const real_t>(data[i]);}
    p= CartesianGrid(vs);
}
