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

#ifndef OOK_STREAM_OBSERVER_H_
#define OOK_STREAM_OBSERVER_H_

#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "linalg/norms.hpp"

#include "ook/state.h"

namespace ook{

/// \brief Observer that outputs the state to the given stream.
template <typename Stream>
struct stream_observer
{
    stream_observer(Stream& stream)
    :
        stream_(stream)
    {}

    template <typename State>
    void operator()(const State& s)
    {
        if (s.tag == state_tag::init)
        {
            stream_ << std::setw(6) << "\nn"
                    << std::setw(6) << "nfev"
                    << std::scientific
                    << std::setw(14) << "a"
                    << std::setw(14) << "fx"
                    << std::setw(14) << "max ||dfx||"
                    << std::setw(14) << "max ||dx||\n";
        }

        if (s.tag == state_tag::iterate)
        {
            stream_ << std::setw(6) << s.iteration
                    << std::setw(6) << s.nfev
                    << std::scientific
                    << std::setw(14) << s.a
                    << std::setw(14) << s.fx
                    << std::setw(14) << linalg::norm_infinity(s.dfx)
                    << std::setw(14) << linalg::norm_infinity(s.dx);
        }

        if (s.tag == state_tag::final)
        {
            stream_ << "\nstatus : " << s.msg << "\n";
            stream_ << std::setw(8) << "iter"
                    << std::setw(8) << "nfev"
                    << std::setw(16) << "fx"
                    << std::setw(16) << "max ||dfx||"
                    << std::setw(16) << "max ||dx||\n";

            stream_ << std::setw(8) << s.iteration
                    << std::setw(8) << s.nfev
                    << std::scientific
                    << std::setw(16) << s.fx
                    << std::setw(16) << linalg::norm_infinity(s.dfx)
                    << std::setw(16) << linalg::norm_infinity(s.dx) << "\n";
        }
    }

private:
    Stream& stream_;
};

} // ns ook

#endif
