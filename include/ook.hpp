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

#ifndef OOK_HPP_
#define OOK_HPP_

#include "ook/matrix.hpp"
#include "ook/vector.hpp"

#include "ook/message.hpp"
#include "ook/options.hpp"
#include "ook/version.hpp"

#include "ook/stream_observer.hpp"

#include "ook/bfgs.hpp"
#include "ook/lbfgs.hpp"
#include "ook/line_search_method.hpp"
#include "ook/newton.hpp"
#include "ook/nonlinear_cg.hpp"
#include "ook/steepest_descent.hpp"

#include "ook/line_search/backtracking.hpp"
#include "ook/line_search/conditions.hpp"
#include "ook/line_search/mcsrch.hpp"

#include "ook/finite_differences/backward_difference.hpp"
#include "ook/finite_differences/central_difference.hpp"
#include "ook/finite_differences/finite_differences.hpp"
#include "ook/finite_differences/forward_difference.hpp"

/// \mainpage Introduction
///
/// OOK is a library for developing optimisation algorithms.
///
/// Currently implemented optimisers are line search based methods:
///
/// - Steepest descent ook::steepest_descent.
/// - Fletcher-Reeves ook::fletcher_reeves.
/// - BFGS ook::bfgs.
/// - Newton ook::newton.
/// - LBFGS ook::lbfgs.

#endif
