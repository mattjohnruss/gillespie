#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

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

    // Temp storage

    std::ifstream file;
    std::vector<std::vector<std::vector<double> > > file_data;
    std::string line;

    // Loop over files and read data

    for(unsigned f = 0; f < n_files; f++)
    {
        char file_suffix[30];
        std::sprintf(file_suffix, "%04i.dat", f+1);

        std::cout << file_prefix + file_suffix << std::endl;

        file.open((file_prefix + file_suffix).c_str());

        std::vector<std::vector<double> > data;
        while(std::getline(file,line))
        {
            std::stringstream  lineStream(line);
            std::string        cell;

            std::vector<double> temp_v;
            double temp_d;

            while(std::getline(lineStream,cell,' '))
            {
                // You have a cell
                std::stringstream  cellStream(cell);
                cellStream >> temp_d;
     
                temp_v.push_back(temp_d);
            }

            data.push_back(temp_v);
        }

        file_data.push_back(data);
        file.close();
    }

    //for(unsigned f = 0; f < n_files; f++)
    //{
    //    for(unsigned i = 0; i <= 10000; i++)
    //    {
    //        for(unsigned j = 0; j < 12; j++)
    //        {
    //            std::cout << file_data[f][i][j] << " ";
    //        }
    //        std::cout << std::endl;
    //    }
    //    std::cout << std::endl << std::endl;
    //}


    // Calculate statistics
    
    std::vector<std::vector<double> > mean(10001, std::vector<double>(12, 0.0));
    std::vector<std::vector<double> > mean_of_sqs(10001, std::vector<double>(12, 0.0));
    std::vector<std::vector<double> > variance(10001, std::vector<double>(12, 0.0));

    double inv_n_files = 1./(double)n_files;
    double bias_correction = (double)n_files/(double)(n_files+1);

    std::ofstream mean_file((file_prefix + "_mean.dat").c_str());
    std::ofstream variance_file((file_prefix + "_variance.dat").c_str());

    for(unsigned i = 0; i <= 10000; i++)
    {
        for(unsigned j = 0; j < 12; j++)
        {
            for(unsigned f = 0; f < n_files; f++)
            {
                mean[i][j] += file_data[f][i][j];
                mean_of_sqs[i][j] += pow(data[f][i][j],2);
            }
            mean[i][j] *= inv_n_files;
            mean_of_sqs[i][j] *= inv_n_files;

            mean_file << mean[j] << std::endl;

            variance_file << nodes[j] << " "
                << bias_correction*(mean_of_sqs[j] - pow(mean[j],2)) << std::endl;

        }
        std::cout << std::endl << std::endl;
    }
 }
