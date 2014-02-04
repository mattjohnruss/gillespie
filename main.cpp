#include <iostream>
#include <vector>

#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

int main(int argc, char **argv)
{
    if(!(argc == 6))
    {
        std::cerr << "Usage: " << argv[0] << " n_urns n_init a b n_max\n";
        exit(1);
    }

    // Step 1: Set the initial number of "molecules" of each "species" (in this
    // case, the number of molecules in each urn)
    unsigned n_urns = 0;

    std::stringstream iss(argv[1]);
    iss >> n_urns;

    unsigned n_init = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[2]);

    iss >> n_init;

    // Set the problem parameters a, b
    double a = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[3]);
    iss >> a;

    double b = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[4]);
    iss >> b;

    // Set the maximum number of iterations
    unsigned n_max = 0;
    iss.str("");
    iss.clear();
    iss.str(argv[5]);
    iss >> n_max;

    // Declare the storage for the urns
    std::vector<unsigned> n(n_urns,0);

    // Add n_init molecules to the first urn
    n[0] = n_init;

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

    // First 2n_urns-1 entires are "hop left", next 2n_urns-1 are "hop right"
    std::vector<double> T(2*(n_urns-1));

    // Inefficient
    std::vector<double> probs(2*(n_urns-1));

    // Total of rates
    double T0 = 0;

    // Time
    double time = 0;

    // Time increment
    double dt = 0;

    // add later
    //double T_eject;

    bool event_happened = false;

    for(unsigned i = 0; i < n_max; i++)
    {
        // Reset T0
        T0 = 0;

        event_happened = false;

        for (unsigned j = 0; j < n_urns-1; j++)
        {
            // Calculate the "reaction" rates

            // Hop left
            T[j] = a*n[j+1];

            // Hop right
            T[n_urns-1 + j] = (a+b)*n[j];

            // Add to T_0 (inline for speed)
            T0 += T[j] + T[n_urns-1 + j];
        }

        std::cout << "At time " << time << ":\n";
        std::cout << "\nEvent\t\tProbability:\n";

        // Calculate the probabilities
        for (unsigned j = 0; j < 2*(n_urns-1); j++)
        {
            probs[j] = T[j]/T0;
            std::cout << j << "\t\t" << probs[j] << std::endl;
        }

        std::cout << std::endl;

        // Get two uniformly random numbers
        r1 = uniform_gen();

        // Set next timestep
        dt = 1./T0*log(1./r1);

        // Set up weighted discrete distribution with the probabilities
        boost::random::discrete_distribution<>
            discrete_dist(probs.begin(),probs.end());

        // Randomly choose event
        unsigned event = discrete_dist(rng2);

        std::cout << "T0: " << T0 << "\nEvent: " << event;

        // Perform stuff due to event
        if (event < (n_urns-1))
        {
            unsigned urn = event;
            std::cout << " (hop left from urn "
                      << urn+1 << " to urn " << urn << ")\n\n";
            if(n[urn+1] > 0)
            {
                n[urn+1]--;
                n[urn]++;

                event_happened = true;
            }
            else
            {
                std::cout << "No change since n[event+1] empty!\n";
            }
        }
        else
        {
            unsigned urn = event - (n_urns-1);
            std::cout << " (hop right from urn "
                      << urn << " to urn " << urn+1 << ")\n\n";
            if(n[urn] > 0)
            {
                n[urn]--;
                n[urn+1]++;

                event_happened = true;
            }
            else
            {
                std::cout << "No change since n[event] empty!\n";
            }
        }

        // Output

        std::cout << "Urn:\t\tn:\n";
        for(unsigned j = 0; j < n_urns; j++)
        {
            std::cout << j << "\t\t" << n[j] << std::endl;
        }

        // Increment time with the timestep
        if(event_happened)
        {
            time += dt;
        }
    }

    return 0;
}

