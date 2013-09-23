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
    tol :: String
} deriving (Show)

epsilon = "std::numeric_limits<typename Vector::value_type>::epsilon()"
infinity = "std::numeric_limits<typename Vector::value_type>::infinity()"
neg_infinity = "-std::numeric_limits<typename Vector::value_type>::infinity()"
vector = "std::vector<typename Vector::value_type>"
-- Function "name" "0" "0" "{x0}" "{minima}" "f_min" epsilon,
functionList = [
    --Function "rosenbrock"          "2" "2"  "{-1.2, 1.0}" "{1.0, 1.0}"         "0.0"     epsilon,
    Function "freudenstein_roth"   "2" "2"  "{0.5,  2.0}" "{5.0, 4.0}"         "0.0"     epsilon,
    Function "powell_badly_scaled" "2" "2"  "{0.0,  1.0}" "{1.098e-05, 9.106}" "0.0"     epsilon,
    Function "brown_badly_scaled"  "2" "3"  "{0.0,  1.0}" "{1e6, 2e-06}"       "0.0"     epsilon,
    Function "beale"               "2" "3"  "{1.0,  1.0}" "{3, 0.5}"           "0.0"     epsilon,
    Function "jenrich_sampson"     "2" "10" "{0.3,  0.4}" "{0.2578, 0.2578}"   "124.362" epsilon,

    Function "helical_valley"      "3" "3"  "{-1.0, 0.0, 0.0}" "{1.0, 0.0, 0.0}"                                            "0.0"         epsilon,
    Function "bard"                "3" "15" "{1.0, 1.0, 1.0}"  ("{0.8406, " ++ neg_infinity ++ ", " ++ neg_infinity ++ "}") "17.4286"     epsilon,
    Function "gaussian"            "3" "15" "{0.4, 1.0, 0.0}"  "{}"                                                         "1.12793e-08" epsilon,
    Function "meyer"               "3" "16" "{0.2, 4000, 250}" "{}"                                                         "87.9458"     epsilon,
    Function "gulf_rnd"            "3" "50" "{5, 2.5, 0.15}"   "{50, 25, 1.5}"                                              "0.0"         epsilon,
    Function "box_3d"              "3" "10" "{0, 10, 20}"      "{1, 10, 1}"                                                 "0.0"         epsilon,

    Function "powell_singular"     "4"  "4"  "{3, -1, 0, 1}"                  "{0.0, 0.0, 0.0, 0.0}" "0.0"         epsilon,
    Function "wood"                "4"  "6"  "{-3, -1, -3, -1}"               "{1.0, 1.0, 1.0, 1.0}" "0.0"         epsilon,
    Function "kowalik_osborne"     "4"  "11" "{0.25, 0.39, 0.415, 0.39}"      "{}"                   "3.07505e-04" epsilon,
    Function "brown_dennis"        "4"  "20" "{25, 5, -5, -1}"                "{}"                   "85822.2"     epsilon,
    Function "osborne_1"           "5"  "33" "{0.5, 1.5, -1, 0.01, 0.02}"     "{}"                   "5.46489e-05" epsilon,
    Function "biggs_exp6"          "6"  "13" "{1.0, 2.0, 1.0, 1.0, 1.0, 1.0}" "{1, 10, 1, 5, 4, 3}"  "0.0"         epsilon,

    Function "osborne_2"           "11" "65" "{1.3, 0.65, 0.65, 0.7, 0.6, 3, 5, 7, 2, 4.5, 5.5}" "{}"                  "4.01377e-02" epsilon,
    Function "watson"              "12" "31" (vector ++ "(12, 0.0)")                             "{}"                  "4.72238e-10" epsilon,
    Function "extended_rosenbrock"            "10" "10" "{-1.2, 1.0, -1.2, 1.0, -1.2, 1.0, -1.2, 1.0, -1.2, 1.0, }" "{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, }" "0.0" epsilon,
    Function "extended_powell_singular"       "0" "0" "{}" "{}" "0.0" epsilon,
    Function "penalty_i"                      "10" "11" "{}" "{}" "0.0" epsilon,
    Function "penalty_ii"                     "0" "0" "{}" "{}" "0.0" epsilon,
    Function "variable_dim"                   "0" "0" "{}" "{}" "0.0" epsilon,
    Function "trigonometric"                  "0" "0" "{}" "{}" "0.0" epsilon,
    Function "brown_almost_linear"            "0" "0" "{}" "{}" "0.0" epsilon,
    Function "discrete_boundary_value"        "0" "0" "{}" "{}" "0.0" epsilon,
    Function "discrete_integral_equation"     "0" "0" "{}" "{}" "0.0" epsilon,
    Function "broyden_tridiagonal"            "0" "0" "{}" "{}" "0.0" epsilon,
    Function "broyden_banded"                 "0" "0" "{}" "{}" "0.0" epsilon,
    Function "linear_full_rank"               "0" "0" "{}" "{}" "0.0" epsilon,
    Function "linear_rank_1"                  "0" "0" "{}" "{}" "0.0" epsilon,
    Function "linear_rank_1_with_0_cols_rows" "0" "0" "{}" "{}" "0.0" epsilon,
    Function "chebyquad"                      "0" "0" "{}" "{}" "0.0" epsilon]

upperCase = map toUpper

--replaceFunctions :: Function -> String -> String
replaceIfDef = replace "@IFDEF@" . upperCase . name
replaceName = replace "@NAME@" . name
replaceN = replace "@N@" . n
replaceM = replace "@M@" . m
replaceX0 = replace "@XO@" . x0 
replaceMinima = replace "@MINIMA@" . minima
replaceFMin = replace "@FMIN@" . fmin
replaceTol = replace "@TOL@" . tol

doReplacements :: Function -> String -> String
doReplacements f = (replaceIfDef f . replaceFMin f 
                        . replaceMinima f . replaceX0 f . replaceTol f
                            . replaceM f . replaceN f . replaceName f)

generateHeader :: Function -> IO ()
generateHeader f = writeFile filename filecontents
                    where
                        filename = "../include/test_functions/more_garbow_hillstrom/" ++ name f ++ ".h" 
                        filecontents = (doReplacements f) templateText

main :: IO ()
main = forM_ functionList generateHeader 

