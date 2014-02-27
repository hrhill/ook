#ifndef OOK_STREAM_OBSERVER_H_
#define OOK_STREAM_OBSERVER_H_

namespace ook{

template <typename Stream>
struct stream_observer
{
    stream_observer(Stream& stream)
    :
        stream_(stream)
    {}

    Stream& stream_;
/*
    stream << std::endl;
    stream << std::setw(6) << "n"
           << std::setw(6) << "nfev"
           << std::scientific
           << std::setw(14) << "a"
           << std::setw(14) << "fx"
           << std::setw(14) << "max ||dfx||"
           << std::setw(14) << "max ||dx||" << std::endl;

template <typename Stream, typename X>
void
report(Stream& stream, message msg, int iteration, int nfev_total, int nfev, double a, double fx, const X& dfx, const X& dx)
{
    stream << std::setw(6) << iteration
              << std::setw(6) << nfev_total
              << std::scientific
              << std::setw(14) << a
              << std::setw(14) << fx
              << std::setw(14) << ook::norm_infinity(dfx)
              << std::setw(14) << ook::norm_infinity(dx) << std::endl;
}

template <typename Stream, typename X>
void
final_report(Stream& stream, message msg, int iteration, int nfev_total,double fx, const X& dfx, const X& dx)
{
    stream << "status : " << msg << std::endl;
    stream << std::setw(8) << "iter"
           << std::setw(8) << "nfev"
           << std::setw(16) << "fx"
           << std::setw(16) << "max ||dfx||"
           << std::setw(16) << "max ||dx||" << std::endl;
    stream << std::setw(8) << iteration
           << std::setw(8) << nfev_total
           << std::scientific
           << std::setw(16) << fx
           << std::setw(16) << ook::norm_infinity(dfx)
           << std::setw(16) << ook::norm_infinity(dx) << std::endl;
}
*/
};

} // ns ook

#endif
