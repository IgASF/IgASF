/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <vector>
#include <algebra/eigen.h>
#include <tools/registry.h>

struct  QuadTemplate
{
    std::vector<real_t> nodes;
    std::vector<real_t> weights;
    std::string         name;

    friend void to_json  (Json& j, const QuadTemplate& p);
    friend void from_json(const Json& j, QuadTemplate& p);
};

void to_json  (Json& j, const QuadTemplate& p);
void from_json(const Json& j, QuadTemplate& p);

struct Rules
{
    static QuadTemplate Gauss  (int order);
};



