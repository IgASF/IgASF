/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <algebra/eigen.h>

namespace Eigen {

void from_json(const Json& j, FMatrix& p)
{
    index_t r=j.at("rows").get<index_t>();
    index_t c=j.at("cols").get<index_t>();
    std::vector<real_t> data=j.at("coefs").get<std::vector<real_t>>();

    p.resize(r,c);
    p= CFMap(data.data(),r,c);
}

void from_json(const Json& j, IMatrix& p)
{
    index_t r=j.at("rows").get<index_t>();
    index_t c=j.at("cols").get<index_t>();
    std::vector<index_t> data=j.at("coefs").get<std::vector<index_t>>();

    p.resize(r,c);
    p= CIMap(data.data(),r,c);
}

void from_json(const Json& j, SMatrix& p)
{
    index_t r=j.at("rows").get<index_t>();
    index_t c=j.at("cols").get<index_t>();
    std::vector<real_t>  coefs=j.at("coefs").get<std::vector<real_t>>();
    std::vector<index_t> cols =j.at("col_pos").get<std::vector<index_t>>();
    std::vector<index_t> begs =j.at("row_beg").get<std::vector<index_t>>();

    p=SMatrix(CSMap(r,c, coefs.size(), begs.data(), cols.data(), coefs.data()));
}


void to_json  (Json& j, const FMatrix& p)
{
    j=Json({{"type","matrix"}, {"rows", p.rows()}, {"cols", p.cols()}, {"coefs", std::vector<real_t>(p.data(), p.data()+p.size()) }});
}
void to_json  (Json& j, const FMap& p)
{
    j=Json({{"type","matrix"}, {"rows", p.rows()}, {"cols", p.cols()}, {"coefs", std::vector<real_t>(p.data(), p.data()+p.size()) }});
}
void to_json  (Json& j, const CFMap& p)
{
    j=Json({{"type","matrix"}, {"rows", p.rows()}, {"cols", p.cols()}, {"coefs", std::vector<real_t>(p.data(), p.data()+p.size()) }});
}
void to_json  (Json& j, const IMatrix& p)
{
    j=Json({{"type","matrix"}, {"rows", p.rows()}, {"cols", p.cols()}, {"coefs", std::vector<index_t>(p.data(), p.data()+p.size()) }});
}
void to_json  (Json& j, const IMap& p)
{
    j=Json({{"type","matrix"}, {"rows", p.rows()}, {"cols", p.cols()}, {"coefs", std::vector<index_t>(p.data(), p.data()+p.size()) }});
}
void to_json  (Json& j, const CIMap& p)
{
    j=Json({{"type","matrix"}, {"rows", p.rows()}, {"cols", p.cols()}, {"coefs", std::vector<index_t>(p.data(), p.data()+p.size()) }});
}


void to_json(Json& j, const SMatrix& p)
{
    std::vector<real_t> coefs (p.valuePtr(), p.valuePtr()+p.nonZeros());
    std::vector<index_t> cols(p.innerIndexPtr(), p.innerIndexPtr()+p.nonZeros());
    std::vector<index_t> rows(p.outerIndexPtr(), p.outerIndexPtr()+p.rows()+1);
    j=Json({
        {"type","sparse matrix"},
        {"rows", p.rows()},
        {"cols", p.cols()},
        {"coefs",  coefs},
        {"col_pos", cols},
        {"row_beg", rows}
        });
}
}

