#ifndef OOK_OPTIONS_H_
#define OOK_OPTIONS_H_

#include <cassert>
#include <string>
#include <stdexcept>

namespace ook{

#define ASSERT_CHECK_AND_THROW(predicate) assert((predicate));\
    if (!(predicate))\
        throw std::invalid_argument(std::string(#predicate) + std::string(" must hold."));

struct options{
    options(const double ftolIn, const double gtolIn, const double xtolIn, const double stpminIn, const double stpmaxIn)
    :
        ftol(ftolIn), 
        gtol(gtolIn), 
        xtol(xtolIn), 
        stpmin(stpminIn), 
        stpmax(stpmaxIn),
        max_iteration(10000),
        max_line_search_attempts(20)
    {
        /* Check the input arguments for errors. */
        ASSERT_CHECK_AND_THROW(ftol >= 0.);
        ASSERT_CHECK_AND_THROW(gtol >= 0.); 
        ASSERT_CHECK_AND_THROW(xtol >= 0.);
        ASSERT_CHECK_AND_THROW(stpmin >= 0.); 
        ASSERT_CHECK_AND_THROW(stpmax > stpmin); 
    }
    const double ftol;
    const double gtol;
    const double xtol;
    const double stpmin;
    const double stpmax;
    const uint max_iteration;
    const uint max_line_search_attempts;
};

inline
void
validate_arguments(const double stp, const double g, const options& opts)
{
    ASSERT_CHECK_AND_THROW(!(stp < opts.stpmin));
    ASSERT_CHECK_AND_THROW(!(stp > opts.stpmax));
    ASSERT_CHECK_AND_THROW(!(g > 0.0));
}

#undef ASSERT_CHECK_AND_THROW

} // ns ook

#endif