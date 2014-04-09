#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>

#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "config.h"

int main(int argc, char **argv)
{
    if(!(argc == 10))
    {
        std::cerr << "Usage: " << argv[0]
                  << " n_urns n_init a b inflow outflow t_max outfile interval"
                  << std::endl;
        exit(1);
    }

    // Set the initial number of "molecules" of each "species" (in this case,
    // the number of molecules in each urn)
    unsigned n_urns = 0;

    std::istringstream iss(argv[1]);
    iss >> n_urns;

    // Set the number of molecules initially in the first urn
    unsigned n_init = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[2]);
    iss >> n_init;

    // Set a
    double a = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[3]);
    iss >> a;

    // Set b
    double b = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[4]);
    iss >> b;

    // Set the inflow rate
    double T_inflow = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[5]);
    iss >> T_inflow;

    // Set the outflow rate
    double T_outflow = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[6]);
    iss >> T_outflow;

    // Set the maximum time
    double t_max = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[7]);
    iss >> t_max;

    // Set the output file
    std::string outfile_name;
    iss.str("");
    iss.clear();
    iss.str(argv[8]);
    iss >> outfile_name;

    // Set the output interval
    double output_interval;
    iss.str("");
    iss.clear();
    iss.str(argv[9]);
    iss >> output_interval;

    // Make a file object with the given filename
    std::ofstream outfile(outfile_name.c_str());

    //std::ofstream testoutfile("testoutput.dat");

    // Declare the storage for the urns
    std::vector<unsigned> n(n_urns,100);

    // Add n_init molecules to the first urn
    //n[0] = n_init;

    // Uniform random number generator
    boost::mt19937 rng(time(0) + getpid());
    boost::random::uniform_real_distribution<double> uniform_dist(0,1);

    boost::variate_generator<
        boost::mt19937&,
        boost::random::uniform_real_distribution<double> >
            uniform_gen(rng, uniform_dist);

    // Discrete distribution based on T
    boost::mt19937 rng2(time(0) + getpid()+1);

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

    for(unsigned i = 0; i < n_removal; i++)
    {
        T_removal[i] = 0.;
    }

    // Total of rates
    double T0 = 0;

    // Time
    double time = 0;

    // Time increment
    double dt = 0;

    unsigned total_particles = n_init;

    // Output the initial state

    outfile << 0 << " " << time;
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
    while(time < t_max && total_particles > 0)
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
            if(n[j] > 0)
            {
                T[n_hop_left + n_hop_right + j] = T_removal[j];
                T0 += T_removal[j];
            }
            else
            {
                T[n_hop_left + n_hop_right + j] = 0;
            }

            // Calculate total_particles
            total_particles += n[j];
        }

        // Inflow
        T[n_hop_left + n_hop_right + n_removal] = T_inflow;
        T0 += T_inflow;

        // Outflow (can only remove a particle if the last urn is non-empty)
        if(n[n_urns-1] > 0)
        {
            T[n_hop_left + n_hop_right + n_removal + n_inflow] =
                T_outflow; //*n[n_urns-1];
            T0 += T_outflow; //*n[n_urns-1];
        }
        else
        {
            T[n_hop_left + n_hop_right + n_removal + n_inflow] = 0;
        }

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

        // Randomly choose event
        unsigned event = discrete_dist(rng2);

        // Output must happen here to avoid saving the previous state of the
        // urns

        // Do outputs at appropriate intervals until the time of the last
        // output exceeds the new actual time (t+dt)
        for(unsigned l = k+1; l*output_interval < time+dt; l++)
        {
            outfile << l << " " << l*output_interval;

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

