#ifndef OOK_OPTIMISER_STATUS_H_
#define OOK_OPTIMISER_STATUS_H_

namespace ook{

enum class optimiser_status{
    keep_going,
    problem_with_options,
    problem_with_objective_function,
    max_iterations_reached,
    no_step_possible
};

const char* status_names[] = {"keep_going", 
                              "problem_with_options", 
                              "problem_with_objective_function",
                              "max_iterations_reached",
                              "no_step_possible"};

template <typename Enumeration>
typename std::underlying_type<Enumeration>::type
as_integer(Enumeration const value)
{
    return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}

// Add a pretty printer here for optimiser_status

} //ook

#endif