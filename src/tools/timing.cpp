/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <atomic>

// assembling
std::atomic<long> time_compute_structure=0;
std::atomic<long> time_eval_coef=0;
std::atomic<long> time_eval_bases=0;
std::atomic<long> time_assemble=0;
std::atomic<long> time_macro_split=0;
std::atomic<long> time_macro_setup=0;
std::atomic<long> time_add_macro=0;
std::atomic<long> time_geo_compute=0;
std::atomic<long> time_geo_transform=0;
// matrix free application
std::atomic<long> time_apply_trial;
std::atomic<long> time_apply_kronecker;
