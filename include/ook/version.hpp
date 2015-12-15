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

#ifndef OOK_VERSION_HPP_
#define OOK_VERSION_HPP_

#include <string>

#define OOK_MAJOR_VERSION 0
#define OOK_MINOR_VERSION 10
#define OOK_PATCH_LEVEL 0
#define OOK_VERSION ( OOK_MAJOR_VERSION * 100000 + OOK_MINOR_VERSION * 100 + OOK_PATCH_LEVEL )


namespace ook{

struct version{
    static constexpr int major = OOK_MAJOR_VERSION;
    static constexpr int minor = OOK_MINOR_VERSION;
    static constexpr int patch = OOK_PATCH_LEVEL;

    static std::string
    string(){
        return "v" + std::to_string(major) +
               "." + std::to_string(minor) +
               "." + std::to_string(patch);
    }
};

}

#endif
