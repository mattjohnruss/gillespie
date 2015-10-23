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
    // Number of urns
    unsigned n_urns = 0;

    // Number of particles initially in the first urn
    unsigned n_init = 0;

    // p
    double p = 0;

    // q
    double q = 0;

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

    // Flag for lognormal sinks
    bool lognormal_sinks = false;

    // Desired mean and variance of the lognormal distribution
    double lognormal_mean = 1;
    double lognormal_variance = 1;
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
        ("n_init", po::value<unsigned>(&n_init), "initial number of particles")
        ("hop right,p", po::value<double>(&p)->required(), "hop right rate")
        ("hop left,q", po::value<double>(&q)->required(), "hop left rate")
        ("inflow,a", po::value<double>(&T_inflow)->required(), "inflow rate")
        ("outflow,b", po::value<double>(&T_outflow)->required(), "outflow rate")
        ("t_max,t", po::value<double>(&t_max)->required(), "maximum time")
        ("outfile,f", po::value<std::string>(&outfile_name)->required(), "output file name")
        ("interval,i", po::value<double>(&output_interval), "output interval")
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
        // TODO change this to use exceptions
        exit(1);
    }

    // Check for required args and output any errors
    po::notify(vm);

    // Set the flag if we have mean and variance on the command line
    if(vm.count("lognormal_mean") && vm.count("lognormal_variance"))
    {
        lognormal_sinks = true;
    }
    else if(vm.count("lognormal_mean") || vm.count("lognormal_variance"))
    {
        std::cout << "lognormal mean or variance specified without the other!\n";
        // TODO change this to use exceptions
        exit(1);
    }

    // Set the output interval to t_max if the flag is set
    if(vm.count("last_only"))
    {
        output_interval = t_max;
    }
}

//void generate_correlated_ics(std::vector<unsigned> &n, double &mean, Eigen::Matrix3d &cholesky_cov)
//{
//    if(n.size() != 3)
//    {
//        // Lolcats
//        exit(1);
//    }
//
//    Eigen::Vector3d n_ics;
//    std::random_device rd_ics;
//    std::mt19937 rng_ics(rd_ics());
//    std::normal_distribution<double> dist_ics(0,1);
//    for(unsigned i = 0; i < 3; ++i)
//    {
//        n_ics(i) = dist_ics(rng_ics);
//    }
//
//    // apparently safe to do this in one line in Eigen
//    n_ics = cholesky_cov*n_ics;
//
//    for(unsigned i = 0; i < 3; ++i)
//    {
//        n[i] = static_cast<unsigned>(std::round(mean + n_ics(i)));
//    }
//}

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

    // ***********************
    // Generate correlated ICs
    // ***********************

    //double mean = 10;

    //// Hardcode the precomputed Cholesky decomposition of the covariance matrix
    //Eigen::Matrix3d cholesky_cov;
    //cholesky_cov << 1, 0, 0,
    //                1./2, std::sqrt(3.)/2., 0,
    //                1./2, 1./(2.*std::sqrt(3.)), std::sqrt(2./3);

    //// Generate the correlated ICs
    //generate_correlated_ics(n, mean, cholesky_cov);

    // Add n_init particles to the first urn
    //n[0] = n_init;

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

    // RNG for lognormal distribution for the sink strengths
    std::mt19937 rng_lognormal(rd());

    // Calculate the lognormal parameters from desired mean and variance

    // Calculate the m and s params for the distribution based on mean and var
    double lognormal_m =
        std::log(std::pow(lognormal_mean,2)/std::sqrt(lognormal_variance + std::pow(lognormal_mean,2)));
    double lognormal_s =
        std::sqrt(std::log(1+lognormal_variance/std::pow(lognormal_mean,2)));

    // Lognormal distribution object
    std::lognormal_distribution<double>
        lognormal_dist(lognormal_m, lognormal_s);

    std::mt19937 rng_uniform_sinks(rd());
    //double ep = 0.5;
    //std::uniform_real_distribution<double> uniform_sinks_dist(1. - ep,1. + ep);

    std::mt19937 rng_norm_pert_sinks(rd());
    std::normal_distribution<double> norm_pert_sinks_dist(0.5,0.25);

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
    const unsigned n_hop_left  = n_urns-1;
    const unsigned n_hop_right = n_urns-1;
    const unsigned n_removal   = n_urns;
    const unsigned n_inflow    = 1;
    const unsigned n_outflow   = 1;

    const unsigned n_events = n_hop_left + n_hop_right + n_removal + n_inflow + n_outflow;

    // Make the vector of "rates" of the different events
    std::vector<double> T(n_events);

    std::vector<double> T_removal(n_removal);

    if(lognormal_sinks)
    {
        for(unsigned i = 0; i < n_removal; ++i)
        {
            T_removal[i] = lognormal_dist(rng_lognormal);
        }
    }
    else
    {
        unsigned sink_interval = 10;
        ///double sink_strength = 1;

        //for(unsigned i = 0; i < n_removal; ++i)
        //{
        //    T_removal[i] = 0;
        //    //T_removal[i] = uniform_sinks_dist(rng_uniform_sinks);
        //    //T_removal[i] = norm_pert_sinks_dist(rng_norm_pert_sinks);
        //}

        T_removal[0] = 0.0;

        for(unsigned i = 1; i < n_removal; ++i)
        {
            if(!(i % sink_interval))
            {
                //T_removal[i] = sink_strength;
                //T_removal[i] = norm_pert_sinks_dist(rng_norm_pert_sinks);
                T_removal[i] = 0.5;
                //std::cout << "sink at urn " << i << std::endl;
            }
            else
            {
                T_removal[i] = 0.0;
            }

            //T_removal[i] = 0.;
            //T_removal[i] = uniform_sinks_dist(rng_uniform_sinks);
        }

        T_removal[0] = 0;
        T_removal[n_removal-1] = 0;
    }

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
                // Hop left
                T[j] = q*n[j+1];

                // Hop right
                T[n_hop_left + j] = p*n[j];

                // Add hopping contributions to T0
                T0 += T[j] + T[n_hop_left + j];
            }

            // Assign the rates of removal events (have to do this in the loop
            // since removal can only occur if the urn is non-empty!)
            //if(n[j] > 0)
            //{
                T[n_hop_left + n_hop_right + j] = T_removal[j]*n[j];
            //    T[n_hop_left + n_hop_right + j] = T_removal[j];
                T0 += T[n_hop_left + n_hop_right + j];
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
            T0 += T[n_hop_left + n_hop_right + n_removal + n_inflow];
        //}
        //else
        //{
        //    T[n_hop_left + n_hop_right + n_removal + n_inflow] = 0;
        //}

        // Get a uniformly random number
        r1 = uniform_dist(rng_uniform);

        // Set next timestep from exponential distribution using inverse
        // transform sampling (inverse of cdf of exp dist.)
        dt = 1./T0*log(1./r1);

        // Set up weighted discrete distribution with the rates
        // discrete_distribution normalises the weights so we don't have to
        //std::discrete_distribution<>
        //    discrete_dist(T.begin(),T.end());

        // Randomly choose event
        //event = discrete_dist(rng_discrete);

        // More efficient way of selecting event a la Anderson 2007

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

        // Hop left
        if(event < n_hop_left)
        {
            urn = event;

            --n[urn+1];
            ++n[urn];
        }
        // Hop right
        else if(event < (n_hop_left + n_hop_right))
        {
            urn = event - n_hop_left;

            --n[urn];
            ++n[urn+1];
        }
        // Removal
        else if(event < (n_hop_left + n_hop_right + n_removal))
        {
            urn = event - (n_hop_left + n_hop_right);

            --n[urn];

            // Update running total of particles
            --total_particles;
        }
        else if(event < (n_hop_left + n_hop_right + n_removal + n_inflow))
        {
            urn = 0;

            ++n[urn];

            // Update running total of particles
            ++total_particles;
        }
        else if(event < (n_hop_left + n_hop_right + n_removal + n_inflow + n_outflow))
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

