#ifndef LINE_SEARCH_H_
#define LINE_SEARCH_H_

#include <string>
#include <stdexcept>
#include <cassert>

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
        stpmax(stpmaxIn)
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
};

#undef ASSERT_CHECK_AND_THROW

int dcsrch(const double finit, const double ginit, double& stp, double f, double g, std::string& task, const options& opts);
int dcstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, 
            double& stp, const double& fp, const double& dp, bool& brackt, const double& stpmin, const double& stpmax);

#endif
