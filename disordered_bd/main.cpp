#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <random>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

// Declare parameters in the Params namespace
namespace Params
{
    // Initial number of particles
    unsigned n_init = 0;

    // Birth rate
    double T_birth = 0;

    // Death rate
    double T_death = 0;
    //double T_death_per_particle = 0;

    // Maximum time
    double t_max = 0;

    // Output file
    std::string outfile_name;

    // Output interval
    double output_interval = 0.01;
}

void parse_command_line(int &argc, char ** &argv)
{
    using namespace Params;

    // Options description object
    po::options_description desc("Allowed options");

    // Add all the possible cmd line options to the desc cobject
    desc.add_options()
        ("help,h", "produce help message")
        ("n_init,n", po::value<unsigned>(&n_init)->required(), "initial number of particles")
        ("birth,b", po::value<double>(&T_birth)->required(), "birth rate")
        //("death,d", po::value<double>(&T_death_per_particle)->required(), "death rate per particle")
        ("death,d", po::value<double>(&T_death)->required(), "death rate")
        ("t_max,t", po::value<double>(&t_max)->required(), "maximum time")
        ("outfile,f", po::value<std::string>(&outfile_name)->required(), "output file name")
        ("interval,i", po::value<double>(&output_interval), "output interval")
        ("last_only,l","only output the last timestep");

    // Map for the the variables
    po::variables_map vm;

    // Parse the command line and store the args in the variables map
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // If the help arg is given, output the options and exit
    if(vm.count("help"))
    {
        std::cout << desc << '\n';
        // TODO change this to use exceptions
        exit(1);
    }

    // Check for required args and output any errors
    po::notify(vm);

    // Set the output interval to t_max if the flag is set
    if(vm.count("last_only"))
    {
        output_interval = t_max;
    }
}

int main(int argc, char **argv)
{
    using namespace Params;

    // Parse the command line args and store them in the Params namespace
    parse_command_line(argc, argv);

    // Make a file object with the given filename
    std::ofstream outfile(outfile_name.c_str());

    // Declare the storage for the urns and initialise to n_init
    int n = n_init;

    // "cryptographically" RNG used to seed the other RNGs.
    // Need this because other seed methods such as getpid() + time() etc can
    // produce duplicate seeds
    std::random_device rd;

    // RNG for uniform random distribution
    std::mt19937 rng_uniform(rd());

    // Uniform distribution object
    std::uniform_real_distribution<double> uniform_dist(0,1);

    // RNG for discrete distribution based on T
    // (dist constructed inside time loop because it is different for each
    // timestep)
    std::mt19937 rng_discrete(rd());

    // RNG for birth parameter
    //std::mt19937 rng_uniform_birth(rd());

    // Uniform distribution object for birth parameter
    //std::uniform_real_distribution<double> uniform_dist_birth(100,1000);

    // Get a random birth rate
    //T_birth = uniform_dist_birth(rng_uniform_birth);

    // Storage for uniformly random number used for timestep
    double r1 = 0;

    // Total of rates
    double T0 = 0;

    // Time
    double time = 0;

    // Time increment
    double dt = 0;

    // Output the initial state

    //outfile << time;
    //outfile << " " << n;
    //outfile << std::endl;

    // Index of the last output
    unsigned k = 0;

    // Timestepping loop
    while(time < t_max)
    {
        // Reset T0
        T0 = 0;

        T0 += T_birth;
        T0 += T_death;

        //T_death = T_death_per_particle*n;
        //T0 += T_death_per_particle*n;

        // Get a uniformly random number
        r1 = uniform_dist(rng_uniform);

        // Set next timestep from exponential distribution using inverse
        // transform sampling (inverse of cdf of exp dist.)
        dt = 1./T0*log(1./r1);

        // Set up weighted discrete distribution with the rates
        // discrete_distribution normalises the weights so we don't have to
        std::discrete_distribution<>
            discrete_dist({T_birth, T_death});

        // Randomly choose event
        unsigned event = discrete_dist(rng_discrete);

        // Output must happen here to avoid saving the previous state of the
        // urns

        // Do outputs at appropriate intervals until the time of the last
        // output exceeds the new actual time (t+dt)
        for(unsigned l = k+1; l*output_interval < time+dt && l*output_interval <= t_max; ++l)
        {
            outfile << l*output_interval << " " << n << std::endl;
            //outfile << n << std::endl;

            // Update k to l (by the end of the loop k should equal the
            // index of the last output performed). Only NEEDS to be done on
            // the last iteration but it would (probably) be more expensive
            // to do a test than to just set k each time.
            k = l;
        }

        // Perform stuff due to event

        if(event == 0)
        {
            ++n;
        }
        else if(event == 1)
        {
            --n;
        }
        else
        {
            std::cerr << "ERROR\n";
            exit(1);
        }

        // Increment time with the timestep
        time += dt;
    }

    outfile.close();

    return 0;
}

