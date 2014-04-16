module TestHeader where

templateText :: String
templateText =
    "// Copyright 2013 Harry Hill\n\
    //\n\
    // This file is part of ook.\n\
    //\n\
    // ook is free software: you can redistribute it and/or modify\n\
    // it under the terms of the GNU Lesser General Public License as published by\n\
    // the Free Software Foundation, either version 3 of the License, or\n\
    // (at your option) any later version.\n\
    //\n\
    // ook is distributed in the hope that it will be useful,\n\
    // but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
    // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
    // GNU Lesser General Public License for more details.\n\
    //\n\
    // You should have received a copy of the GNU Lesser General Public License\n\
    // along with ook.  If not, see <http://www.gnu.org/licenses/>.\n\
    \n\
    #ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_@IFDEF@_H_\n\
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
    \    typedef Matrix matrix_type;\n\
    \    typedef typename vector_type::value_type real_type;\n\
    \\n\
    \    std::tuple<real_type, vector_type, matrix_type>\n\
    \    operator()(const vector_type& x) const\n\
    \    {\n\
    \        real_type f(1.0);\n\
    \        vector_type df(@N@, 1.0);\n\
    \        matrix_type d2f(@N@, @N@, 0.0);\n\
    \        return std::make_tuple(f, df, d2f);\n\
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
    \template <typename Vector, typename Matrix>\n\
    \typename Vector::value_type\n\
    \@NAME@<Vector, Matrix>::f_min = @FMIN@;\n\
    \\n\
    \template <typename Vector, typename Matrix>\n\
    \typename Vector::value_type\n\
    \@NAME@<Vector, Matrix>::tolerance = @TOL@;\n\
    \\n\
    \template <typename Vector, typename Matrix>\n\
    \std::vector<typename Vector::value_type>\n\
    \@NAME@<Vector, Matrix>::minima = @MINIMA@;\n\
    \\n\
    \template <typename Vector, typename Matrix>\n\
    \std::vector<typename Vector::value_type>\n\
    \@NAME@<Vector, Matrix>::x0 = @XO@;\n\
    \\n\n\
    \} // ns test_functions\n\
    \} // ns ook\n\n\
    \#endif\n"
