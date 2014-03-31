// Copyright 2013 Harry Hill
//
// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

#ifndef OOK_H_
#define OOK_H_

#include "ook/norms.h"
#include "ook/message.h"
#include "ook/version.h"
#include "ook/options.h"
#include "ook/state.h"

#include "ook/stream_observer.h"

#include "ook/factorisations/cholesky.h"
#include "ook/factorisations/ldlt.h"
#include "ook/factorisations/gmw81.h"
#include "ook/factorisations/tools.h"

#include "ook/line_search_method.h"
#include "ook/newton.h"
#include "ook/steepest_descent.h"
#include "ook/bfgs.h"
#include "ook/fletcher_reeves.h"

#include "ook/line_search/conditions.h"
#include "ook/line_search/more_thuente.h"
#include "ook/line_search/backtracking.h"

#include "ook/finite_differences/central_difference.h"
#include "ook/finite_differences/finite_differences.h"
#include "ook/finite_differences/forward_difference.h"
#include "ook/finite_differences/backward_difference.h"

#endif
