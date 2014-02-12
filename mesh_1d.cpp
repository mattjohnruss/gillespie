#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>

#include "mesh_1d.h"

template <typename type_, typename xtype_>
Mesh_1D<type_, xtype_>::Mesh_1D() :
    nodes_(),
    vars_(),
    num_nodes_(0),
    num_vars_(0)
{
}

template <typename type_, typename xtype_>
Mesh_1D<type_, xtype_>::~Mesh_1D()
{
}

template <typename type_, typename xtype_>
type_& Mesh_1D<type_, xtype_>::operator()(const std::size_t& node_index,
        const std::size_t& var_index)
{
    return vars_[num_vars_*node_index + var_index];
}

template <typename type_, typename xtype_>
const type_& Mesh_1D<type_, xtype_>::operator()(
        const std::size_t& node_index, const std::size_t& var_index) const
{
    return vars_[num_vars_*node_index + var_index];
}

template <typename type_, typename xtype_>
void Mesh_1D<type_, xtype_>::load_data(const std::string& path)
{
    load_data(path.c_str());
}

template <>
void Mesh_1D<double>::load_data(const char* path)
{
    std::ifstream data;
    std::string line_data;

    data.open(path);
    if (data.is_open())
    {
        /*if (!std::getline(data, line_data))
          {
          data.close();
          throw Exceptions::Generic("Data IO",
          "Unspecified error during reading of data file");
          }
          if (line_data != "# Mesh_1D")
          {
          data.close();
          throw Exceptions::Generic("Data IO",
          "Data must be in Mesh_1D output format with "
          "\"# Mesh_1D\" at head of file");
          }*/

        bool possible_block(false), num_vars_found(false);
        double temp(0.);

        num_nodes_ = 0;
        nodes_.clear();
        num_vars_ = 0;
        vars_.clear();

        while (std::getline(data, line_data))
        {
            if (line_data[0] == '#')
                continue;
            if (line_data == "")
            {
                if (possible_block)
                    break;
                else
                    possible_block = true;
                continue;
            }
            possible_block = false;
            std::istringstream iss(line_data);
            if (!(iss >> temp))
                std::cout << "Error\n";
                //throw Exceptions::Generic("Data IO",
                //        "Error during reading of data file, probably "
                //        "missing data");
            nodes_.push_back(temp);
            while (iss >> temp)
            {
                vars_.push_back(temp);
                if (!num_vars_found)
                    ++num_vars_;
            }
            num_vars_found = true;
        }
        data.close();

        num_nodes_ = nodes_.size();
    }
    else
        std::cout << "Error\n";
        //throw Exceptions::Generic("Data IO", "Unable to open data file");
}

template <typename type_, typename xtype_>
void Mesh_1D<type_, xtype_>::dump(std::ostream& out,
        const std::size_t& skip) const
{
    out << "# Mesh_1D\n";
    out << "# Node ";
    out << "# Data\n";
    for (std::size_t i = 0; i < num_nodes_; i += skip)
    {
        out << nodes_[i];
        for (std::size_t j = 0; j < num_vars_; ++j)
            out << ' ' << vars_[num_vars_*i + j];
        out << '\n';
    }
    out << "\n# Number of nodes: " << num_nodes_ << '\n';
    out << "# Skip size: " << skip << '\n';
    out << "# Number of variables: " << num_vars_ << "\n\n\n";
}

template class Mesh_1D<double>;
