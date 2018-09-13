/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <maps/geomap.h>
#include <maps/rational_map.h>
#include <maps/bc_map.h>

// default implementation gets points from grid
// and computes on all points
FMatrix GeoMap::evaluate(const CartesianGrid &grid) const
{
  auto pts=grid.toPoints();
  return this->evaluate(pts);
}

FMatrix GeoMap::jacobian(const CartesianGrid &grid) const
{
  auto pts=grid.toPoints();
  return this->evaluate(pts);
}

REGISTER_SUBTYPE(GeoMap,RationalMap);
REGISTER_SUBTYPE(GeoMap,BasisCoefficientMap);

