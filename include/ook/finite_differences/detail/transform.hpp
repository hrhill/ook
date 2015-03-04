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

#ifndef OOK_FINITE_DIFFERENCES_DETAIL_TRANSFORM_HPP_
#define OOK_FINITE_DIFFERENCES_DETAIL_TRANSFORM_HPP_

namespace ook{
namespace finite_differences{
namespace detail{

template <typename In, typename Out, typename F>
void transform(const In& x, Out& y, F f)
{
	#pragma omp parallel for shared(x, y) firstprivate(f)
	for(size_t i=0; i < y.size(); ++i){
		y[i] = f(x[i]);
	}
}

}
} // ns finite differences
} // ns ook

#endif
