#include <tuple>
#include <vector>

typedef std::vector<double> vector_t;

enum class task{
    evaluate,    
    converged,
    warn
};

struct options{
    double ftol;
    double gtol;
    double stpmin;
    double stpmax;
};

void
csrch(double phi, double dphi, const options& opts, task& msg)
{
    task msg = task::evaluate;

    return std::make_pair(step, msg);
}

/*
    // Call phi
    std::tuple<double, double> phi_x = phi(stp);
    options opts = {1e-03, 1e-01, 1e-10, 1e-16, 1e16};


*/


/*
task t;
do{
    // get descent direction
    state = descent_direction(objective_function, state); // potentially keep state within the optimiser class, e.g., B in BFGS.
    auto phi = [&state](const double alpha) -> std::pair<double, vector_t>{
        return objective_function(state.x + alpha * state.p);
    };
    do{
        std::tie(stp, msg) = line_search(state, opts);
        if (msg == task::evaluate)
            state = phi(stp)
    }while()

}while(msg == task::get_new_direction);
*/