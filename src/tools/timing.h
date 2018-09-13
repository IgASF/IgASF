/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include<atomic>
#include<chrono>

// all in microseconds
// i.e. 1e-6 s
extern std::atomic<long> time_compute_structure;
extern std::atomic<long> time_eval_coef;
extern std::atomic<long> time_geo_compute;
extern std::atomic<long> time_geo_transform;
extern std::atomic<long> time_eval_bases;
extern std::atomic<long> time_assemble;
extern std::atomic<long> time_macro_setup;
extern std::atomic<long> time_add_macro;

template <class T>
__attribute__((always_inline)) inline void DoNotOptimize( T &value) {
  asm volatile("" : "+m"(const_cast<T &>(value)));
}

enum TIMES
{
    compute_structure=0,
    eval_coef,
    geo_compute,
    geo_transform,
    eval_bases,
    assemble,
    macro_setup,
    add_macro,
    real,
    cpu,
    count
};
