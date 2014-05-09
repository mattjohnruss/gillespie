#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc, char **argv)
{
    // Parse command line arguments
    // ----------------------------

    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " n_files file_prefix\n";
        exit(1);
    }

    unsigned n_files;
    std::string file_prefix;

    std::istringstream args_stream(argv[1]);
    args_stream >> n_files;

    args_stream.str("");
    args_stream.clear();
    args_stream.str(argv[2]);
    args_stream >> file_prefix;

    // Storage for data reading
    // ------------------------

    std::ifstream file;
    std::string line_data;

    double temp;

    std::vector<double> nodes;
    std::vector<std::vector<std::vector<double> > > data(n_files,std::vector<std::vector<double> > ());

    exit(0);
    bool first_done = false;
    unsigned n_nodes = 0;

    // Loop over all the files and read the data
    // -----------------------------------------

    for(unsigned i = 0; i < n_files; i++)
    {
        char file_suffix[30];
        std::sprintf(file_suffix, "%04i.dat", i+1);

        std::cout << file_prefix + file_suffix << std::endl;

        file.open((file_prefix + file_suffix).c_str());

        while(std::getline(file, line_data))
        {
            std::istringstream data_stream(line_data);

            data_stream >> temp;

            if(!first_done)
            {
                nodes.push_back(temp);
            }

            data_stream >> temp;
            data[i].push_back(temp);
        }

        if(!first_done)
        {
            n_nodes = nodes.size();
            for(unsigned j = 1; j < n_files; j++)
            {
                data[j].reserve(n_nodes);
            }
        }

        first_done = true;
        file.close();
    }

    // ------------------------------------------------------------------------

    for(unsigned i = 0; i < n_files; i++)
    {
        for(unsigned j = 0; j < 12; j++)
        {
            std::cout << data[i][j] << " ";
        }

        std::cout << std::endl;
    }
    exit(0);

    // Calculate statistics
    // --------------------

    double inv_n_files = 1./(double)n_files;
    double bias_correction = (double)n_files/(double)(n_files+1);

    std::vector<double> mean(n_nodes,0.);
    std::vector<double> mean_of_sqs(n_nodes,0.);
    std::vector<double> mean_of_sqs_of_res(n_nodes,0.);

    std::vector<double> variance(n_nodes,0.);
    std::vector<double> variance_of_res(n_nodes,0.);

    std::ofstream mean_file((file_prefix + "_mean.dat").c_str());
    std::ofstream mean_of_res_file((file_prefix + "_mean_of_res.dat").c_str());

    std::ofstream variance_file((file_prefix + "_variance.dat").c_str());
    std::ofstream variance_of_res_file((file_prefix + "_variance_of_res.dat").c_str());

    for(unsigned j = 0; j < n_nodes; j++)
    {
        for(unsigned i = 0; i < n_files; i++)
        {
            mean[j] += data[i][j];
            mean_of_sqs[j] += pow(data[i][j],2);
            mean_of_sqs_of_res[j] += pow(data[i][j] - (1. - 1./50.*nodes[j]),2);
        }

        mean[j] *= inv_n_files;
        mean_of_sqs[j] *= inv_n_files;
        mean_of_sqs_of_res[j] *= inv_n_files;

        mean_file << nodes[j] << " "
                  << mean[j] << std::endl;

        mean_of_res_file << nodes[j] << " "
                         << mean[j] - (1. - 1./50.*nodes[j]) << std::endl;

        variance_file << nodes[j] << " "
                      << bias_correction*(mean_of_sqs[j] - pow(mean[j],2)) << std::endl;

        variance_of_res_file << nodes[j] << " "
                             << bias_correction*(mean_of_sqs_of_res[j] -
                                     pow(mean[j] - (1. - 1./50.*nodes[j]),2))
                             << std::endl;
    }

    mean_file.close();
    mean_of_res_file.close();

    variance_file.close();
    variance_of_res_file.close();

    return 0;
}
