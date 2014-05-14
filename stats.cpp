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

    // Loop over the files
    for(unsigned f = 0; f < n_files; f++)
    {
        // Temp c string for suffix of current file
        char file_suffix[30];

        // Print the correct suffix into file_suffix
        std::sprintf(file_suffix, "%04i.dat", f+1);

        // TODO remove debug output
        std::cout << file_prefix + file_suffix << std::endl;

        // Open the file
        file.open((file_prefix + file_suffix).c_str());

        // Output error if the file wasn't opened successfully
        if(!file.is_open())
        {
            std::cout << "File " << file_prefix + file_suffix << " not open!\n";
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

            // TODO remove debug output
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

    // Calculate statistics
    // --------------------

    // Storage for the number of fields
    // Get it from the first line of the first file (assume all lines of all
    // files have the same number of fields)
    unsigned n_fields = data[0][0].size();

    // Storage for mean, mean_of_squares and variance
    std::vector<std::vector<double> > mean(n_nodes, std::vector<double>(n_fields,0.));
    std::vector<std::vector<double> > mean_of_squares(n_nodes, std::vector<double>(n_fields,0.));
    std::vector<std::vector<double> > variance(n_nodes, std::vector<double>(n_fields,0.));

    double inv_n_files = 1./(double)n_files;
    double bias_correction = (double)n_files/(double)(n_files+1);

    std::ofstream mean_file((file_prefix + "_mean.dat").c_str());
    std::ofstream variance_file((file_prefix + "_variance.dat").c_str());

    for(unsigned j = 0; j < n_nodes; j++)
    {
        mean_file << nodes[j] << " ";
        variance_file << nodes[j] << " ";

        for(unsigned field = 0; field < n_fields; field++)
        {
            for(unsigned f = 0; f < n_files; f++)
            {
                mean[j][field] += data[f][j][field];
                mean_of_squares[j][field] += pow(data[f][j][field],2);
            }

            mean[j][field] *= inv_n_files;
            mean_of_squares[j][field] *= inv_n_files;

            mean_file << mean[j][field] << " ";

            variance_file << bias_correction*(mean_of_squares[j][field] - pow(mean[j][field],2))
                          << " ";
        }

        mean_file << std::endl;
        variance_file << std::endl;
    }

    mean_file.close();

    variance_file.close();

    //unsigned n_lines;
    //unsigned n_fields;

    //for(unsigned f = 0; f < n_files; f++)
    //{
    //    n_lines = data[f].size();

    //    for(unsigned line = 0; line < n_lines; line++)
    //    {
    //        n_fields = data[f][line].size();

    //        std::cout << nodes[line] << " ";

    //        for(unsigned field = 0; field < n_fields; field++)
    //        {
    //            std::cout << data[f][line][field] << " ";
    //        }

    //        std::cout << std::endl;
    //    }

    //    std::cout << "\n-------------------------------------------\n\n";
    //}

    return 0;
}
