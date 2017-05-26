/* This file is part of IGsim, a simple c++ simulation library.
 *
 * Copyright (C) 2016 Like Ma <milkpku@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public License
 * v. 2.0. If a copy of the MPL was not distributed with this file, You can
 * obtain one at http://mozilla.org/MPL/2.0/.
 * This should *NOT* be contained in a IGSIM_*_H ifdef, since it may be defined
 * differently based on when it is included
 */

#ifndef IGSIM_UNITTEST_DEFINES_H
#define IGSIM_UNITTEST_DEFINES_H

#define REPEAT_N 20
#define NUM_EQ(A, B) EXPECT_LT(abs(A-B), 1e-10)

#endif //IGSIM_UNITTEST_DEFINES_H
