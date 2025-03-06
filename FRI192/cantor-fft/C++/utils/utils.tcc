#include <utils/utils.hpp>

template <typename FieldT>
bool field_trace_binary(const FieldT &element)
{
    FieldT trace = element;
    FieldT element_squared = element;
    for (size_t i = 1; i < FieldT::extension_degree(); ++i)
    {
        element_squared = element_squared.squared();
        trace += element_squared;
    }
    return !trace.is_zero(); // === is one
}

template <typename FieldT>
bool field_trace_binary(const FieldT &element, size_t degree)
{
    FieldT trace = element;
    FieldT element_squared = element;
    for (size_t i = 1; i < degree; ++i)
    {
        element_squared = element_squared.squared();
        trace += element_squared;
    }
    return !trace.is_zero(); // === is one
}

template <typename FieldT>
std::vector<FieldT> cantor_basis(size_t m)
{  
    FieldT beta_m;
    std::vector<FieldT> basis;
    basis.reserve(m);
    basis.resize(m);
    size_t max_bm_degree = (FieldT::extension_degree() == 192) ? 64 : FieldT::extension_degree();

    if(FieldT::extension_degree() == 192){
        FieldT alpha(7); // 7 is the primitive element in GF(2^192)
        // gamma is a primitive element in 2 ^ 64
        FieldT gamma = alpha ^ (libff::bigint<4>("340282366920938463481821351505477763073")); // 340282366920938463481821351505477763073 = (2**192 - 1)//(2**64 - 1)
        while (!field_trace_binary(gamma, 64)){
            gamma *= gamma;
        }
        beta_m = gamma;
    }
    else {
            do
        {
            beta_m = FieldT::random_element();
        } while (!field_trace_binary(beta_m));
    }


    for (size_t i = max_bm_degree; i > m + 1; --i)
    {
        beta_m += beta_m.squared();
    }

    for (int i = m - 1; i >= 0; --i)
    {
        beta_m += beta_m.squared();
        basis[i] = beta_m;
    }

    return basis;
}

template <typename FieldT>
std::vector<FieldT> divide  (std::vector<FieldT> &g, 
                            const std::vector<size_t> &nz_S, 
                            const size_t input_size, 
                            const size_t offset){
    std::vector<FieldT> q;
    
}

template <typename T>
void my_print_vector(const std::vector<T> &v)
{
    std::cout <<"{ ";
    for (auto const& elem : v)
    {
        std::cout << elem << " ";
    }
    std::cout << "}" << std::endl;
}

