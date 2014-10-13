#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>

unsigned number_of_digits(unsigned n)
{
    unsigned digits = 0;

    do
    {
        n /= 10;
        digits++;
    }
    while (n != 0);

    return digits;
}

int main(int argc, char **argv)
{
    // Parse command line arguments
    // ----------------------------

    // If we don't have the required number of arguments print a message and exit
    if(argc != 3 && argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " n_files file_prefix [leading_zeros]\n";
        exit(1);
    }

    // The number of files (sample size)
    unsigned n_files;

    // Prefix that each file name begins with
    std::string file_prefix;

    // Stream for program arguments
    std::istringstream args_stream(argv[1]);

    // Put the first argument into n_files
    args_stream >> n_files;

    // Clear the stream and put the second argument into file_prefix
    args_stream.str("");
    args_stream.clear();
    args_stream.str(argv[2]);
    args_stream >> file_prefix;

    // Storage for data reading
    // ------------------------

    // Storage for data. Indices: files, lines, fields
    std::vector<std::vector<std::vector<double> > > data(n_files);

    // Storage for nodes
    std::vector<double> nodes;

    // File stream
    std::ifstream file;

    // String to store each line
    std::string line_data;

    // First iteration of file loop?
    bool first_it = true;

    // Number of nodes
    unsigned n_nodes = 0;

    // Temporary storage for data
    double temp_cell = 0;

    // Temporary storage for lines
    std::vector<double> temp_line;

    // Storage for the number of zeros to pad the file names with
    unsigned n_padding_zeros = 0;

    if(argc == 4)
    {
        args_stream.str("");
        args_stream.clear();
        args_stream.str(argv[3]);
        args_stream >> n_padding_zeros;
    }
    else
    {
        n_padding_zeros = number_of_digits(n_files);
    }

    // C format string for padding with the correct number of zeros
    char zeros_fmt_str[8];
    std::sprintf(zeros_fmt_str, "%%0%ii.dat", n_padding_zeros);

    // Loop over the files
    for(unsigned f = 0; f < n_files; f++)
    {
        // Temp c string for suffix of current file
        char file_suffix[30];

        // Print the correct suffix into file_suffix
        std::sprintf(file_suffix, zeros_fmt_str, f+1);

        // Open the file
        file.open((file_prefix + file_suffix).c_str());

        // Output error if the file wasn't opened successfully
        if(!file.is_open())
        {
            std::cerr << "File " << file_prefix + file_suffix << " not open!\n";
            exit(1);
        }

        // Loop over the lines of the current file and put the line in line_data
        while(std::getline(file, line_data))
        {
            // Clear the temporary line before reading the new line
            temp_line.clear();

            // Put the line into a stream
            std::istringstream line_stream(line_data);

            // Get the first item out of the stream ("node location")
            line_stream >> temp_cell;

            // If we're parsing the first file add the node location to the
            // nodes vector. Otherwise ignore it and move on to the data.
            if(first_it)
            {
                nodes.push_back(temp_cell);
            }

            // Read the other items (the actual fields) from the line stream
            while(line_stream >> temp_cell)
            {
                temp_line.push_back(temp_cell);
            }

            // Push back the temporary line
            data[f].push_back(temp_line);
        }

        // Allocate the memory for the rest of the data now that we know the
        // sizes of everything after the first iteration
        if(first_it)
        {
            // Set the number of nodes now that we know it
            n_nodes = nodes.size();

            // Print the number of nodes in the first file
            std::cout << "There are " << n_nodes << " nodes in the first file\n";

            // Loop over the files again (starting at the second file)
            for(unsigned f2 = 1; f2 < n_files; f2++)
            {
                // Reserve memory for each file
                data[f2].reserve(n_nodes);
            }
        }

        // We've now done the first iteration so set first_it to false
        first_it = false;

        // Close the file
        file.close();
    }

    std::cout << "Read " << n_files << " files\n";

    // Calculate statistics
    // --------------------

    // Storage for the number of fields
    // Get it from the first line of the first file (assume all lines of all
    // files have the same number of fields)
    unsigned n_fields = data[0][0].size();

    // Storage for mean and covariance
    std::vector<std::vector<double> > mean(n_nodes, std::vector<double>(n_fields,0.));
    std::vector<std::vector<std::vector<double> > > covariance(
            n_nodes, std::vector<std::vector<double> >(n_fields, std::vector<double>(n_fields,0.)));

    // Calculate 1/n_files
    double inv_n_files = 1./(double)n_files;
    double inv_n_files_m1 = 1./((double)n_files - 1.);

    // File streams for mean, covariance and correlation
    std::ofstream mean_file((file_prefix + "_mean.dat").c_str());
    std::ofstream covariance_file((file_prefix + "_covariance.dat").c_str());
    std::ofstream correlation_file((file_prefix + "_correlation.dat").c_str());

    mean_file << std::setprecision(10);
    covariance_file << std::setprecision(10);
    correlation_file << std::setprecision(10);

    // Loop over the nodes (timesteps)
    for(unsigned node = 0; node < n_nodes; node++)
    {
        // Output the node locations (once per node)
        mean_file << nodes[node];

        // Loop over the fields
        for(unsigned field = 0; field < n_fields; field++)
        {
            // Loop over the files
            for(unsigned f = 0; f < n_files; f++)
            {
                // Add the value from each file (i.e. sample) to the mean and mean_of_squares
                mean[node][field] += data[f][node][field]*inv_n_files;
            }

            // Loop over the fields again for the covariance
            for(unsigned field2 = 0; field2 < n_fields; field2++)
            {
                // Loop over the files again --- have to do this separately
                // from the above loop since we need the means in the covariance
                // calculation
                for(unsigned f = 0; f < n_files; f++)
                {
                   covariance[node][field][field2] +=
                       (data[f][node][field] - mean[node][field]) *
                       (data[f][node][field2] - mean[node][field2])*inv_n_files_m1;
                }

                // At this point in the loops, covariance[node][field][field2] has
                // all the data it needs and just needs to be output

                // Output the covariance
                covariance_file << covariance[node][field][field2] << " ";
            }

            // New line bewteen rows of the cov matrix
            covariance_file << std::endl;

            // Output the mean
            mean_file << " " << mean[node][field];
        }

        // New line
        mean_file << std::endl;

        // Double new line between cov matrices at different timesteps
        covariance_file << std::endl << std::endl;
    }

    // Close the mean and variance files
    mean_file.close();
    covariance_file.close();

    // Post-process the covariance into the correlation
    // ------------------------------------------------

    // Loop over the nodes (timesteps)
    for(unsigned node = 0; node < n_nodes; node++)
    {
        // Loop over the fields
        for(unsigned field = 0; field < n_fields; field++)
        {
            // Loop over the fields again
            for(unsigned field2 = 0; field2 < n_fields; field2++)
            {
                correlation_file
                    << covariance[node][field][field2]/(
                            sqrt(covariance[node][field][field])*sqrt(covariance[node][field2][field2]))
                    << " ";
            }

            // New line bewteen rows of the corr matrix
            correlation_file << std::endl;
        }
        // Double new line between corr matrices at different timesteps
        correlation_file << std::endl << std::endl;
    }

    correlation_file.close();

    // Cross-correlation
    // -----------------
    //
    // is it really the cross-corr? It is the correlation coefficient between n_i at time t
    // and n_j at time t+k, as a function of k

    std::ofstream cross_corr_1_file("cross_corr_1.dat");

    // The reference node corresponding to time t, calculated manually from t=30, dt=0.01
    // 1/0.01 * (30 + 0.01)
    unsigned ref_node = 3001;

    // The number of lag steps
    unsigned n_k = n_nodes - ref_node;

    // Storage for the cross correlation
    std::vector<std::vector<double> > cross_corr_1(n_k, std::vector<double>(n_fields,0.));

    // Loop over the lag
    for(unsigned k = 0; k < n_k; ++k)
    {
        cross_corr_1_file << k << " ";

        // Loop over the fields
        for(unsigned field = 0; field < n_fields; ++field)
        {
            // Loop over the files
            for(unsigned f = 0; f < n_files; ++f)
            {
                // Add the contribution from each file
                cross_corr_1[k][field] +=
                    (data[f][ref_node][0] - mean[ref_node][0]) *
                    (data[f][ref_node + k][field] - mean[ref_node + k][field])*inv_n_files_m1;
            }

            // Divide by the the variances of fields 1 and 2
            cross_corr_1[k][field] /=
                (sqrt(covariance[ref_node][0][0] *
                      covariance[ref_node][field][field]));

            // Output the cross correlation
            cross_corr_1_file << cross_corr_1[k][field] << " ";
        }

        // New line after each lag
        cross_corr_1_file << std::endl;
    }

    // Close the correlation file
    cross_corr_1_file.close();

    return 0;
}
