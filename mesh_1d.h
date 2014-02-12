#ifndef MESH_1D_H_
#define MESH_1D_H_

#include <vector>
#include <string>
#include <iostream>

template<typename type_,typename xtype_=double>
class Mesh_1D
{
    protected:
        std::vector<xtype_> nodes_;
        std::vector<type_> vars_;
        std::size_t num_nodes_, num_vars_;

    public:
        Mesh_1D();
        ~Mesh_1D();

        type_& operator()(const std::size_t& node_index,
                const std::size_t& var_index);

        const type_& operator()(const std::size_t& node_index,
                const std::size_t& var_index) const;

        void load_data(const std::string& path);
        void load_data(const char* path);

        const type_ linear_interpolant(const xtype_& x,
                const std::size_t var) const;

        void dump(std::ostream &out = std::cout,
                const std::size_t& skip = 1) const;
};

#endif
