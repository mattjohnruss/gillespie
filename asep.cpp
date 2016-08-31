#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <random>

#include <boost/program_options.hpp>
//#include <Eigen/Dense>

namespace po = boost::program_options;

// Declare parameters in the Params namespace
namespace Params
{
    // Urn capacity
    unsigned omega = 1;

    // Number of urns
    unsigned n_urns = 0;

    // Number of particles initially in the first urn
    unsigned n_init = 0;

    // p
    double p = 0;

    // Inflow rate
    double T_inflow = 0;

    // Outflow rate
    double T_outflow = 0;

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
        ("n_urns,n", po::value<unsigned>(&n_urns)->required(), "number of urns")
        ("omega,c", po::value<unsigned>(&omega), "urn capacity")
        ("n_init", po::value<unsigned>(&n_init), "initial number of particles")
        ("hop right,p", po::value<double>(&p)->required(), "hop right rate")
        ("inflow,a", po::value<double>(&T_inflow)->required(), "inflow rate")
        ("outflow,b", po::value<double>(&T_outflow)->required(), "outflow rate")
        ("t_max,t", po::value<double>(&t_max)->required(), "maximum time")
        ("outfile,f", po::value<std::string>(&outfile_name)->required(), "output file name")
        ("output_interval,i", po::value<double>(&output_interval), "output interval")
        ("last_only,l", "only output the last timestep");

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

    //std::ofstream testoutfile("testoutput.dat");

    // Declare the storage for the urns
    std::vector<unsigned> n(n_urns,0);

    // Add n_init particles to each urn
    for(unsigned i = 0; i < n_urns; ++i)
    {
        n[i] = n_init;
    }

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

    // Storage for uniformly random number used for timestep
    double r1 = 0;

    // Storage for uniformly random number used for selecting event
    double r2 = 0;

    // Index of event to be performed
    unsigned event = 0;

    // Temporary variable used for selecting event
    double sum_temp = 0;

    // First n_urns-1 entires are "hop left", next n_urns-1 are "hop right",
    // next n_urns are removal, then 1 inflow, 1 outflow
    const unsigned n_hop_right = n_urns-1;
    const unsigned n_inflow    = 1;
    const unsigned n_outflow   = 1;

    const unsigned n_events = n_hop_right + n_inflow + n_outflow;

    // Make the vector of "rates" of the different events
    std::vector<double> T(n_events);

    // **** TEST
    //std::ofstream sinksfile((outfile_name + "-sinks").c_str());
    //for(unsigned i = 0; i < n_urns; ++i)
    //{
    //    sinksfile << T_removal[i] << '\n';
    //}
    //sinksfile.close();

    // Total of rates
    double T0 = 0;

    // Time
    double time = 0;

    // Time increment
    double dt = 0;

    unsigned total_particles = n_init;

    // Output the initial state

    outfile << time;
    for(unsigned j = 0; j < n_urns; ++j)
    {
        outfile << " " << n[j];
    }
    outfile << std::endl;

    //testoutfile << time;
    //for(unsigned j = 0; j < n_urns; ++j)
    //{
    //    testoutfile << " " << n[j];
    //}
    //testoutfile << std::endl;

    // Index of the last output
    unsigned k = 0;

    // Index of urn, used when selecting events
    unsigned urn = 0;

    // Timestepping loop
    while(time < t_max)
    {
        // Reset T0
        T0 = 0;

        // Reset total_particles
        total_particles = 0;

        // Calculate the "reaction" rates for all events
        for(unsigned j = 0; j < n_urns; ++j)
        {
            // There are only n_urns-1 move events in each direction
            if(j < n_urns-1)
            {
                // Hop right
                T[j] = p*(n[j]*(omega-n[j+1]))/((double) omega*omega);

                // Add hopping contributions to T0
                T0 += T[j];
            }

            // Calculate total_particles
            total_particles += n[j];
        }

        // Inflow
        T[n_hop_right] = T_inflow*(omega-n[0])/((double) omega);
        T0 += T[n_hop_right];

        // Outflow
        T[n_hop_right + n_inflow] = T_outflow*n[n_urns-1]/((double) omega);
        T0 += T[n_hop_right + n_inflow];

        // Get a uniformly random number
        r1 = uniform_dist(rng_uniform);

        // Set next timestep from exponential distribution using inverse
        // transform sampling (inverse of cdf of exp dist.)
        dt = 1./T0*log(1./r1);

        // Efficient way of selecting event a la Anderson 2007

        // Draw uniform random number
        r2 = uniform_dist(rng_uniform);

        // Set temp sum variable to zero
        sum_temp = 0;

        // Set event varialbe to zero
        event = 0;

        for(; ; ++event)
        {
            sum_temp += T[event];
            if(sum_temp > T0*r2)
            {
                break;
            }
        }

        //std::cout << event << std::endl;

        // Output must happen here to avoid saving the previous state of the
        // urns

        // Do outputs at appropriate intervals until the time of the last
        // output exceeds the new actual time (t+dt)
        for(unsigned l = k+1; l*output_interval < time+dt && l*output_interval <= t_max; ++l)
        {
            outfile << l*output_interval;

            for(unsigned j = 0; j < n_urns; ++j)
            {
                outfile << " " << n[j];
            }

            outfile << std::endl;

            // Update k to l (by the end of the loop k should equal the
            // index of the last output performed). Only NEEDS to be done on
            // the last iteration but it would (probably) be more expensive
            // to do a test than to just set k each time.
            k = l;
        }

        // Perform stuff due to event

        // Hop right
        if(event < n_hop_right)
        {
            urn = event;

            --n[urn];
            ++n[urn+1];
        }
        else if(event < (n_hop_right + n_inflow))
        {
            urn = 0;

            ++n[urn];

            // Update running total of particles
            ++total_particles;
        }
        else if(event < (n_hop_right + n_inflow + n_outflow))
        {
            urn = n_urns-1;

            --n[urn];

            // Update running total of particles
            --total_particles;
        }
        else
        {
            std::cerr << "ERROR\n";
            exit(1);
        }

        // Increment time with the timestep
        time += dt;

        //// Full output for testing

        //testoutfile << time;
        //for(unsigned j = 0; j < n_urns; j++)
        //{
        //    testoutfile << " " << n[j];
        //}

        //testoutfile << std::endl;
    }

    outfile.close();
    //testoutfile.close();

    return 0;
}
