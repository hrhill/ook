module TestHeader where

templateText :: String
templateText = 
    "#ifndef @IFDEF@_H_\n\
    \#define @IFDEF@_H_\n\
    \\n\
    \#include <vector>\n\
    \\n\
    \template <typename Vector>\n\
    \struct @NAME@\n\
    \{\n\
    \    typedef Vector vector_type;\n\
    \    typedef typename vector_type::value_type real_type;\n\
    \\n\
    \    std::tuple<real_type, vector_type>\n\
    \    operator()(const vector_type& x) const\n\
    \    {\n\
    \        real_type f(0.0);\n\
    \        vector_type df(@N@, 0.0);
    \        @BODY@\n\
    \        return std::make_pair(f, df);\n
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
    \\n\
    \#endif\n"