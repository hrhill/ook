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

#include <string>

namespace ook{
struct version{
    static constexpr int major = 0;
    static constexpr int minor = 4;
    static constexpr int patch = 0;

    static std::string
    string(){
        return "v" + std::to_string(major) +
               "." + std::to_string(minor) +
               "." + std::to_string(patch);
    }
};

}

#endif
