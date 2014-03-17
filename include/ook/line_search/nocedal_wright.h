#ifndef OOK_LINE_SEARCH_NOCEDAL_WRIGHT_H_
#define OOK_LINE_SEARCH_NOCEDAL_WRIGHT_H_

namespace ook{
namespace line_search{

struct nocedal_wright{
    template <typename F, typename T, typename Options>
    static
    std::tuple<message, T, T, T>
    search(F phi, T phi0, T dphi0, T a, const Options& opts){

    }
};

}// line_search
}// ook

#endif
