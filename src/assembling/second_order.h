/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/mmatrix.h>
#include <bases/basis.h>
#include <maps/geomap.h>
#include <assembling/model.h>


struct EqCoef
{
    EqCoef()  : dim(0), A(0,0), B(0,1), C(0) {}
    EqCoef(int d,FMatrix a, FMatrix b, real_t c)
        : dim(d), A(a), B(b), C(c) {}

    bool hasA()const;
    bool hasB()const;
    bool hasC()const;

    int     dim;
    FMatrix A;
    FMatrix B;
    real_t  C;
};
void to_json  (Json& j, const EqCoef& p);
void from_json(const Json& j, EqCoef& p);

struct SecondOrderModel;
void to_json  (Json& j, const SecondOrderModel& p);
void from_json(const Json& j, SecondOrderModel& p);

struct SecondOrderModel :public Model
{
    EqCoef coefs;
    GeoMap *geo;

    SecondOrderModel(EqCoef c, GeoMap *g) : coefs(c), geo(g) {}
    SecondOrderModel(EqCoef c) : coefs(c), geo(nullptr) {}
    SecondOrderModel(SecondOrderModel &&) = default;
    SecondOrderModel& operator=(SecondOrderModel &&o)
    { std::swap(coefs,o.coefs); geo=std::move(o.geo); return *this;}

    private:
        std::vector<Part> initWithGeo(const TensorQuadrature  &) const;
        std::vector<Part> initNoGeo  (const TensorQuadrature  &) const;
    virtual std::vector<Part> initParts(const TensorQuadrature  &) const;
    virtual void done     ( ) const;

    static thread_local std::vector<real_t> m_data;
    virtual Json toJson() const {Json r;  to_json(r,*this); return r;}
};
