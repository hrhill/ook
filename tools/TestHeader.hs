module TestHeader where

templateText :: String
templateText = 
    "#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_@IFDEF@_H_\n\
    \#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_@IFDEF@_H_\n\
    \\n\
    \#include <tuple>\n\
    \#include <limits>\n\
    \#include <vector>\n\
    \\n\
    \namespace ook{\n\
    \namespace test_functions{\n\n\
    \template <typename Vector, typename Matrix>\n\
    \struct @NAME@\n\
    \{\n\
    \    typedef Vector vector_type;\n\
    \    typedef typename vector_type::value_type real_type;\n\
    \\n\
    \    std::tuple<real_type, vector_type>\n\
    \    operator()(const vector_type& x) const\n\
    \    {\n\
    \        real_type f(1.0);\n\
    \        vector_type df(@N@, 1.0);\n\
    \        return std::make_tuple(f, df);\n\
    \    }\n\
    \\n\
    \    static const int n = @N@;\n\
    \    static const int m = @M@;\n\
    \    static real_type f_min;\n\
    \    static real_type tolerance;\n\
    \    static std::vector<real_type> minima;\n\
    \    static std::vector<real_type> x0;\n\
    \};\n\
    \\n\
    \template <typename Vector>\n\
    \typename Vector::value_type\n\
    \@NAME@<Vector>::f_min = @FMIN@;\n\
    \\n\
    \template <typename Vector>\n\
    \typename Vector::value_type\n\
    \@NAME@<Vector>::tolerance = @TOL@;\n\
    \\n\
    \template <typename Vector>\n\
    \std::vector<typename Vector::value_type>\n\
    \@NAME@<Vector>::minima = @MINIMA@;\n\
    \\n\
    \template <typename Vector>\n\
    \std::vector<typename Vector::value_type>\n\
    \@NAME@<Vector>::x0 = @XO@;\n\
    \\n\n\
    \} // ns test_functions\n\    
    \} // ns ook\n\n\
    \#endif\n"