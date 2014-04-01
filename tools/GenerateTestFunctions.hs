import Data.Char (toUpper)
import Control.Monad (forM_)
import Data.String.Utils (replace)

import TestHeader (templateText)

data Function = Function
    { name :: String
    , n :: String
    , m :: String
    , x0 :: String
    , minima :: String
    , fmin :: String
    , tol :: String
} deriving (Show)

epsilon = "std::numeric_limits<typename Vector::value_type>::epsilon()"
infinity = "std::numeric_limits<typename Vector::value_type>::infinity()"
neg_infinity = "-std::numeric_limits<typename Vector::value_type>::infinity()"
vector = "std::vector<typename Vector::value_type>"

rosenbrock = Function
    { name = "rosenbrock"
    , n = "2"
    , m = "2"
    , x0 = "{-1.2, 1.0}"
    , minima = "{1.0, 1.0}"
    , fmin = "0.0"
    , tol = epsilon}

freudensteinRoth = Function
    {name = "freudenstein_roth"
    , n = "2"
    , m = "2"
    , x0 = "{0.5,  2.0}"
    , minima = "{5.0, 4.0}"
    , fmin = "0.0"
    , tol = epsilon}

powellSingular = Function
    { name = "powell_badly_scaled"
    , n =  "2"
    , m =  "2"
    , x0 = "{0.0,  1.0}"
    , minima = "{1.098e-05, 9.106}"
    , fmin = "0.0"
    , tol = epsilon}

brownBadlyScaled = Function
    {name = "brown_badly_scaled"
    , n = "2"
    , m = "3"
    , x0 = "{0.0,  1.0}"
    , minima = "{1e6, 2e-06}"
    , fmin = "0.0"
    , tol =  epsilon}

beale = Function
    { name = "beale"
    , n = "2"
    , m = "3"
    , x0 = "{1.0,  1.0}"
    , minima = "{3, 0.5}"
    , fmin = "0.0"
    , tol = epsilon}

jenrichSampson = Function
    { name = "jenrich_sampson"
    , n = "2"
    , m = "10"
    , x0 = "{0.3,  0.4}"
    , minima = "{0.2578, 0.2578}"
    , fmin = "124.362"
    , tol = epsilon}


functionList = [
    rosenbrock,
    freudensteinRoth,
    powellSingular,
    brownBadlyScaled,
    beale,
    jenrichSampson]

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
        path = "./test_headers/"
        filename = path ++ name f ++ ".h"
        filecontents = (doReplacements f) templateText

main :: IO ()
main = forM_ functionList generateHeader

