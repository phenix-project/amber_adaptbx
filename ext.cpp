#include <vector>
#include <string>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/shared_ptr.hpp>
#include "getmdgxfrc.h"
#include <scitbx/vec3.h>
#include <iostream>
#include <scitbx/array_family/ref_reductions.h>


//Function to call mdgx main routine.
void callMdgx (std::vector<double>& sites_cart, std::vector<double>& gradients,
               std::vector<double>& target, std::string prmtop, std::string crd )
{
        const char * p = prmtop.c_str();
        const char * c = crd.c_str();
        getmdgxfrc(p,c, sites_cart.data(), target.data(), gradients.data());
}

//Function to convert double flex array to vector of doubles
std::vector<double> ExtractVec (scitbx::af::const_ref<double> const& sites_cart){
        std::vector<double> xyz_flat;
        for (size_t i_seq = 0; i_seq < sites_cart.size(); i_seq++) {
                xyz_flat.push_back(sites_cart[i_seq]);
        }
        return xyz_flat;
}

//Function to print out entire vector of doubles
void printVec (std::vector<double> Vec) {
        std::cout << "[ ";
        for (unsigned int i=0; i < Vec.size(); ++i) {
                std::cout << Vec[i] << ", ";
        }
        std::cout << "]\n";
}


//boost::python registration
#include <boost/python.hpp>
BOOST_PYTHON_MODULE(amber_adaptbx_ext)
{
        using namespace boost::python;
    def("callMdgx", &callMdgx);
    def("ExtractVec", &ExtractVec);
    def("printVec", &printVec);
}
