#ifndef PARTICLE_SWARM_H_
#define PARTICLE_SWARM_H_

#include <vector>
#include <random>
#include <utility>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <exception>

namespace smartodds{
namespace optim{
namespace pso{

template <typename X>
struct particle{
    typedef X value_type;

    particle(){}
    particle(const X& x)
    :
        position(x)
    {}

    X position;
    X velocity;
    X best_position;
    typename X::value_type f_position;
    typename X::value_type f_best_position;

    friend
    bool operator<(const particle& x, const particle& y){
        return x.f_position < y.f_position;
    }

    friend
    std::ostream& operator<<(std::ostream& out, const particle& x){
        return out << x.f_position;
    }
};

struct Identity{
    template <typename X>
    const X& operator()(const X& x){
        // print
        return x;
    }
};

template <typename T>
T
bound_value(const T& x, const T& lower, const T& upper){
    return std::min(std::max(x, lower), upper);
}

template <typename T>
T 
bound_vector(const T& x, const T& lower, const T& upper){
    const int n = x.size();
    T bx(x.size());
    for (int i = 0; i < n; ++i){
        bx(i) = bound_value(x(i), lower(i), upper(i));
    }
    return bx;
}

template <typename X>
X
generate_vector(std::mt19937& rng, const X& lb, const X& ub)
{
    assert(ub.size() == lb.size());
    X v(ub.size());
    auto urand = std::bind(std::uniform_real_distribution<double>(0.0, 1.0), std::ref(rng));
    int idx = 0;
    std::generate(v.begin(), v.end(), [&](){
                                        const double u = lb(idx) + (ub(idx) - lb(idx)) * urand();
                                        ++idx;
                                        return u;
                                    });
    return v;
}

struct parameters{
    const int n_particles;
    const double w;
    const double phi_pos;
    const double phi_swm;

    friend
    std::ostream& operator<<(std::ostream& out, const parameters& p){
        return out << "{ n_particles = " << p.n_particles 
                   << ", w  = " << p.w
                   << ", phi_pos = " << p.phi_pos
                   << ", phi_swm = " << p.phi_swm
                   << " }";
    }
};

} // ns pso

template <typename F, typename X>
std::pair<int, X>
particle_swarm(F f, const X& lb, const X& ub, const int max_iterations, 
                const pso::parameters& params = {20, 0.75, 2.05, 2.05})
{
    using namespace pso;
    typedef typename X::value_type real_type;

    // Initialize the swarm using some random values
    const int n_particles = params.n_particles;// > 2 * lb.size() + 1 ? params.n_particles : 2 * lb.size() + 1;
    std::vector<particle<X>> particles(n_particles);
    std::mt19937 rng(std::time(0));

    real_type eps = sqrt(std::numeric_limits<real_type>::epsilon());

    std::for_each(particles.begin(), particles.end(), [&](particle<X>& p){
        p.position = generate_vector<X>(rng, lb, ub);
        p.velocity = generate_vector<X>(rng, lb - ub,  ub - lb);
        p.best_position = p.position;
        p.f_position = f(p.position);
        p.f_best_position = p.f_position;
    });

    particle<X> best = *std::min_element(particles.begin(), particles.end());

    X d = ub - lb;
    std::transform(d.begin(), d.end(), d.begin(), [](const real_type& xi){ 
                                                        return fabs(xi); 
                                                    });
    auto bound_velocity = [&d](const X& x){
        return bound_vector(x, X(-d), d);
    };
    auto bound_position = [&lb, &ub](const X& x){
        return bound_vector(x, lb, ub);
    };

    // Loop section, continue will some criteria not satisfied
    std::string line = "--------------------------------------------------------------------------------";
    std::cout << std::endl << std::endl;
    std::cout << std::right;
    std::cout << std::setw(6) << "i" <<
                 std::setw(16) << "f(x)" <<
                 std::setw(16) << "f_min - f_max" << std::endl << line << std::endl;

    int iteration = 0;
    
    real_type range_f = std::numeric_limits<real_type>::max();
    real_type delta_x = std::numeric_limits<real_type>::max();

    do {
        std::for_each(particles.begin(), particles.end(), 
            [&best, f, bound_velocity, bound_position, params, &rng](particle<X>& p){

                auto runif = std::bind(std::uniform_real_distribution<real_type>(0.0, 1.0), rng);
                const double r_pos = runif();
                const double r_swm = runif();

                p.velocity = params.w * p.velocity + params.phi_pos * r_pos * (p.best_position - p.position) 
                                                   + params.phi_swm * r_swm * (best.position - p.position);
                //p.velocity = bound_velocity(p.velocity);
                p.position += p.velocity;

                p.f_position = f(p.position);
                if (p.f_position < p.f_best_position){
                    p.best_position = p.position;
                    p.f_best_position = p.f_position;
                }
                if (p.f_position < best.f_best_position){
                    best = p;
                }
            });

        auto minmax = minmax_element(particles.begin(), particles.end());
        range_f = minmax.second->f_position - minmax.first->f_position;

        std::cout << std::setw(6) << iteration
                  << std::setw(16) << std::setprecision(8) << std::scientific << best.f_best_position
                  << std::setw(16) << std::setprecision(8) << std::scientific << range_f << " " << best.position << std::endl;
        /*
        std::cout << "all" << std::endl;
        std::for_each(particles.begin(), particles.end(), [](const particle<X>& p){
            std::cout << p.position << std::endl;
        });
        std::cout << "best " << best.position << std::endl;
        */
        ++iteration;
    }while (iteration <= max_iterations and range_f > eps);

    return std::make_pair(0, best.position);
}

} // ns optim
} // ns smartodds

#endif