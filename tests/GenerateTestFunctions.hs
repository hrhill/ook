import Data.Char (toUpper)
import Control.Monad (forM_)
import Data.String.Utils (replace)

import TestHeader (templateText)

data Function = Function{
    name :: String,
    n :: String,
    m :: String,
    x0 :: String,
    minima :: String,
    fmin :: String,
    tol :: String,
    body :: String
} deriving (Show)

epsilon = "std::numeric_limits<double>::epsilon()"

rosenbrockBody = "    const real_type f1 = 10.0 * (x(1) - x(0) * x(0));\n\
                 \    const real_type f2 = 1.0 - x(0);\n\
                 \    f = f1 * f1 + f2 * f2;\n\
                 \    df(0) = -4.0 * 10.0 * f1 * x(0) - 2.0 * f2;\n\
                 \    df(1) = -2.0 * 10.0 * f1;\n"

freudensteinrothBody = "    const real_type x1 = x(0);\n\
                        \    const real_type x2 = x(1);\n\
                        \    const real_type f1 = -13 + x1 + ((5 - x2) * x2 - 2) * x2;\n\
                        \    const real_type f2 = -29 + x1 + ((x2 + 1) * x2 - 14) * x2;\n\
                        \    f = f1 * f1 + f2 * f2;\n\
                        \    df(0) = 2 * f1 + 2 * f2;\n\
                        \    df(1) = 2 * f1 * (10 * x2 - 3 * x2 * x2 - 2)\n\
                        \              + 2 * f2 * ( 3 * x2 * x2 + 2 * x2 - 14);\n"

rosenbrock = Function "rosenbrock" "2" "2" "{-1.2, 1}" "{1.0, 1.0}" "0.0" epsilon rosenbrockBody
freudenstein_roth = Function "freudenstein_roth" "2" "2" "{0.5, 2}" "{5.0, 4.0}" "0.0" "1e-04" freudensteinrothBody

functionList = [rosenbrock]

upperCase = map toUpper

replaceIfDef = replace "@IFDEF@" . upperCase . name
replaceName = replace "@NAME@" . name
replaceN = replace "@N@" . n
replaceM = replace "@M@" . m
replaceX0 = replace "@XO@" . x0 
replaceMinima = replace "@MINIMA@" . minima
replaceFMin = replace "@FMIN@" . fmin
replaceTol = replace "@TOL@" . fmin
replaceBody = replace "@BODY@" . body

doReplacements f = (replaceIfDef f . replaceBody f . replaceFMin f 
                        . replaceMinima f . replaceX0 f . replaceTol f
                            . replaceM f . replaceN f . replaceName f)

generateHeader f = writeFile filename filecontents
                    where
                        filename = "./test_functions/" ++ name f ++ ".h" 
                        filecontents = (doReplacements f) templateText

main :: IO ()
main = forM_ functionList generateHeader 