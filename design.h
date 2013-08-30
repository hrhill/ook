#include <tuple>

struct real_type{};
struct vector_type{};
struct matrix_type{};

namespace ook{
    // concepts
    // 1. Objective function
    struct
    objective_function
    {
        typedef vector_type domain_type;
        typedef real_type range_type;

        range_type
        operator()(const domain_type& x) const {
            // some function of x
            return range_type();
        }
    };

    struct
    objective_gradient_function
    {
        typedef vector_type domain_type;
        typedef real_type range_type;
        typedef vector_type gradient_type;

        std::tuple<range_type, gradient_type>
        operator()(const domain_type& x) const {
            // some function of x
            return std::make_tuple(range_type(), gradient_type());
        }
    };

    struct
    objective_gradient_hessian_function
    {
        typedef vector_type domain_type;
        typedef real_type range_type;
        typedef vector_type gradient_type;
        typedef matrix_type hessian_type;

        std::tuple<range_type, gradient_type, hessian_type>
        operator()(const domain_type& x) const {
            // some function of x
            return std::make_tuple(range_type(), gradient_type(), hessian_type());
        }
    };

    // tags
    struct function_tag{};
    struct function_gradient_tag{};
    struct function_gradient_hessian_tag{};

    template <typename ObjectiveFunction>
    struct objective_function_traits{};

    template <>
    struct objective_function_traits<objective_function>{
        typedef function_tag gradient_information;
    };

    template <>
    struct objective_function_traits<objective_gradient_function>{
        typedef function_gradient_tag gradient_information;
    };

    template <>
    struct objective_function_traits<objective_gradient_hessian_function>{
        typedef function_gradient_hessian_tag gradient_information;
    };


    template <typename ObjectiveFunction>
    struct 
    optimiser_state{
        typedef typename ObjectiveFunction::range_type range_type;
        typedef typename ObjectiveFunction::domain_type domain_type;
        range_type fx;
        domain_type x;
    };

    namespace detail{
        template <typename ObjectiveFunction>
        void
        optimise_impl(ObjectiveFunction object, typename ObjectiveFunction::domain_type init, function_tag){
            optimiser_state<ObjectiveFunction> state;
        }

        template <typename ObjectiveFunction>
        void
        optimise_impl(ObjectiveFunction object, typename ObjectiveFunction::domain_type init, function_gradient_tag){
            optimiser_state<ObjectiveFunction> state;
        }

        template <typename ObjectiveFunction>
        void
        optimise_impl(ObjectiveFunction object, typename ObjectiveFunction::domain_type init, function_gradient_hessian_tag){
            optimiser_state<ObjectiveFunction> state;
        }
    } // ns detail

    // top level call for unconstrained optimisation
    template <typename ObjectiveFunction>
    void
    optimise(ObjectiveFunction object, typename ObjectiveFunction::domain_type init){
        // Based on the traits, we can decide which scheme to pick
        typedef typename objective_function_traits<ObjectiveFunction>::gradient_information gradient_information;
        return detail::optimise_impl(object, init, gradient_information());
    }
}

int main(){
    ook::optimise(ook::objective_function(), vector_type());
}
