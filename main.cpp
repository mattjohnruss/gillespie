#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include <boost/random.hpp>
#include <boost/random/random_device.hpp> // boost/random.hpp doesn't include this header
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv)
{
    // Declare parameter variables so we can pass them to add_options() below

    // Number of urns
    unsigned n_urns = 0;

    // Number of particles initially in the first urn
    unsigned n_init = 0;

    // a
    double a = 0;

    // b
    double b = 0;

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

    // Desired mean and variance of the lognormal distribution
    double lognormal_mean = 1;
    double lognormal_variance = 1;

    // Options description object
    po::options_description desc("Allowed options");

    // Add all the possible cmd line options to the desc cobject
    desc.add_options()
        ("help,h", "produce help message")
        ("n_urns,n", po::value<unsigned>(&n_urns)->required(), "number of urns")
        ("n_init", po::value<unsigned>(&n_init), "initial number of particles")
        ("a,a", po::value<double>(&a)->required(), "diffusion rate")
        ("b,b", po::value<double>(&b)->required(), "advection rate")
        ("inflow,i", po::value<double>(&T_inflow)->required(), "inflow rate")
        ("outflow,o", po::value<double>(&T_outflow)->required(), "outflow rate")
        ("t_max,t", po::value<double>(&t_max)->required(), "maximum time")
        ("outfile,f", po::value<std::string>(&outfile_name)->required(), "output file name")
        ("interval,d", po::value<double>(&output_interval), "output interval")
        ("lognormal_mean,m", po::value<double>(&lognormal_mean), "mean of the lognormal distribution")
        ("lognormal_variance,v", po::value<double>(&lognormal_variance), "variance of the lognormal distribution")
        ("last_only,l","only output the last timestep");

    // Map for the the variables
    po::variables_map vm;

    // Parse the command line and store the args in the variables map
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // If the help arg is given, output the options and exit
    if(vm.count("help"))
    {
        std::cout << desc << '\n';
        return 1;
    }

    // Check for required args and output any errors
    po::notify(vm);

    // Flag for lognormal sinks
    bool lognormal_sinks = false;

    // Set the flag if we have mean and variance on the command line
    if(vm.count("lognormal_mean") && vm.count("lognormal_variance"))
    {
        lognormal_sinks = true;
    }
    else if(vm.count("lognormal_mean") || vm.count("lognormal_variance"))
    {
        std::cout << "lognormal mean or variance specified without the other!\n";
        return 1;
    }

    // Set the output interval to t_max if the flag is set
    if(vm.count("last_only"))
    {
        output_interval = t_max;
    }

    // Make a file object with the given filename
    std::ofstream outfile(outfile_name.c_str());

    //std::ofstream testoutfile("testoutput.dat");

    // Declare the storage for the urns
    std::vector<unsigned> n(n_urns,0);

    // Add n_init particles to the first urn
    //n[0] = n_init;

    // "cryptographically" RNG used to seed the other RNGs.
    // Need this because other seed methods such as getpid() + time() etc can
    // produce duplicate seeds
    boost::random::random_device rd;

    // RNG for uniform random distribution
    boost::mt19937 rng_uniform(rd());

    // Uniform distribution object
    boost::random::uniform_real_distribution<double> uniform_dist(0,1);

    // Variate generator for uniform distribution
    boost::variate_generator<
        boost::mt19937&,
        boost::random::uniform_real_distribution<double> >
            uniform_gen(rng_uniform, uniform_dist);

    // RNG for discrete distribution based on T
    // (dist constructed inside time loop because it is different for each
    // timestep)
    boost::mt19937 rng_discrete(rd());

    // RNG for lognormal distribution for the sink strengths
    boost::mt19937 rng_lognormal(rd());

    // Calculate the lognormal parameters from desired mean and variance

    // Calculate the m and s params for the distribution based on mean and var
    double lognormal_m =
        std::log(std::pow(lognormal_mean,2)/std::sqrt(lognormal_variance + std::pow(lognormal_mean,2)));
    double lognormal_s =
        std::sqrt(std::log(1+lognormal_variance/std::pow(lognormal_mean,2)));

    // Lognormal distribution object
    boost::random::lognormal_distribution<double>
        lognormal_dist(lognormal_m, lognormal_s);

    // Variate generator for lognormal sink distribution
    boost::random::variate_generator<
        boost::mt19937&,
        boost::random::lognormal_distribution<double> >
            lognormal_gen(rng_lognormal, lognormal_dist);

    // Storage for uniformly random number used for timestep
    double r1 = 0;

    // First n_urns-1 entires are "hop left", next n_urns-1 are "hop right",
    // next n_urns are removal, then 1 inflow, 1 outflow
    const unsigned n_hop_left  = n_urns-1;
    const unsigned n_hop_right = n_urns-1;
    const unsigned n_removal   = n_urns;
    const unsigned n_inflow    = 1;
    const unsigned n_outflow   = 1;

    const unsigned n_events = n_hop_left + n_hop_right + n_removal + n_inflow + n_outflow;

    // Make the vector of "rates" of the different events
    std::vector<double> T(n_events);

    // Storage for probabilities for the events - must be stored separately
    // since it depends on T0
    std::vector<double> probs(n_events);

    std::vector<double> T_removal(n_removal);

    if(lognormal_sinks)
    {
        for(unsigned i = 0; i < n_removal; i++)
        {
            T_removal[i] = lognormal_gen();
        }
    }
    else
    {
        for(unsigned i = 0; i < n_removal; i++)
        {
            T_removal[i] = 0.;
        }
    }

    // Total of rates
    double T0 = 0;

    // Time
    double time = 0;

    // Time increment
    double dt = 0;

    unsigned total_particles = n_init;

    // Output the initial state

    outfile << time;
    for(unsigned j = 0; j < n_urns; j++)
    {
        outfile << " " << n[j];
    }
    outfile << std::endl;

    //testoutfile << time;
    //for(unsigned j = 0; j < n_urns; j++)
    //{
    //    testoutfile << " " << n[j];
    //}
    //testoutfile << std::endl;

    // Index of the last output
    unsigned k = 0;

    // Timestepping loop
    while(time < t_max)
    {
        // Reset T0
        T0 = 0;

        // Reset total_particles
        total_particles = 0;

        for(unsigned j = 0; j < n_urns; j++)
        {
            // Calculate the "reaction" rates for all events

            // There are only n_urns-1 move events in each direction
            if(j < n_urns-1)
            {
                // Hop left
                T[j] = a*n[j+1];

                // Hop right
                T[n_hop_left + j] = (a+b)*n[j];

                // Add hopping contributions to T0
                T0 += T[j] + T[n_hop_left + j];
            }

            // Assign the rates of removal events (have to do this in the loop
            // since removal can only occur if the urn is non-empty!)
            //if(n[j] > 0)
            //{
                T[n_hop_left + n_hop_right + j] = T_removal[j]*n[j];
                T0 += T_removal[j];
            //}
            //else
            //{
            //    T[n_hop_left + n_hop_right + j] = 0;
            //}

            // Calculate total_particles
            total_particles += n[j];
        }

        // Inflow
        T[n_hop_left + n_hop_right + n_removal] = T_inflow;
        T0 += T_inflow;

        // Outflow (can only remove a particle if the last urn is non-empty)
        //if(n[n_urns-1] > 0)
        //{
            T[n_hop_left + n_hop_right + n_removal + n_inflow] =
                T_outflow*n[n_urns-1];
            T0 += T_outflow*n[n_urns-1];
        //}
        //else
        //{
        //    T[n_hop_left + n_hop_right + n_removal + n_inflow] = 0;
        //}

        // Calculate the probabilities of all events
        for(unsigned j = 0; j < n_events; j++)
        {
            probs[j] = T[j]/T0;
        }

        // Get two uniformly random numbers
        r1 = uniform_gen();

        // Set next timestep
        dt = 1./T0*log(1./r1);

        // Set up weighted discrete distribution with the probabilities
        boost::random::discrete_distribution<>
            discrete_dist(probs.begin(),probs.end());

        // we don't make a variate generator here because it makes no difference
        // to the values of anything (it's just encapsulation of a dist and a
        // rng) and because we are making a distribution every timestep

        // Randomly choose event
        unsigned event = discrete_dist(rng_discrete);

        // Output must happen here to avoid saving the previous state of the
        // urns

        // Do outputs at appropriate intervals until the time of the last
        // output exceeds the new actual time (t+dt)
        for(unsigned l = k+1; l*output_interval < time+dt && l*output_interval <= t_max; l++)
        {
            outfile << l*output_interval;

            for(unsigned j = 0; j < n_urns; j++)
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

        // Hop left
        if(event < n_hop_left)
        {
            unsigned urn = event;

            n[urn+1]--;
            n[urn]++;
        }
        // Hop right
        else if(event < (n_hop_left + n_hop_right))
        {
            unsigned urn = event - n_hop_left;

            n[urn]--;
            n[urn+1]++;
        }
        // Removal
        else if(event < (n_hop_left + n_hop_right + n_removal))
        {
            unsigned urn = event - (n_hop_left + n_hop_right);

            n[urn]--;

            // Update running total of particles
            total_particles--;
        }
        else if(event < (n_hop_left + n_hop_right + n_removal + n_inflow))
        {
            unsigned urn = 0;

            n[urn]++;

            // Update running total of particles
            total_particles++;
        }
        else if(event < (n_hop_left + n_hop_right + n_removal + n_inflow + n_outflow))
        {
            unsigned urn = n_urns-1;

            n[urn]--;

            // Update running total of particles
            total_particles--;
        }
        else
        {
            std::cerr << "ERROR\n";
            exit(1);
        }

        // Increment time with the timestep
        time += dt;

        //why should we do this???? pretty sure we shouldn't
        //if(total_particles == 0)
        //{
        //    break;
        //}

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

