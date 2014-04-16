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

namespace ook{

/// \brief Observer that outputs the state to the given stream.
template <typename Stream>
struct stream_observer
{
    stream_observer(Stream& stream)
    :
        stream_(stream)
    {}

    Stream& stream_;

    template <typename State>
    void operator()(const State& state){
        stream_ << state << std::endl;
    }
};

} // ns ook

#endif
