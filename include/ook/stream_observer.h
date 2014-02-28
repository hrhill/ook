#ifndef OOK_STREAM_OBSERVER_H_
#define OOK_STREAM_OBSERVER_H_

#include <iostream>
#include <fstream>

namespace ook{

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
