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
