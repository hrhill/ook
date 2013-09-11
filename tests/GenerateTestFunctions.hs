import TestHeader 


data Function = Function{
    name :: String,
    n :: Int,
    m :: Int,
    x0 :: String,
    minima :: String,
    fmin :: Double,
    body :: String
} deriving (Show)


rosenbrockBody = "const real_type f1 = 10.0 * (x(1) - x(0) * x(0));\n\
                 \const real_type f2 = 1.0 - x(0);\n\
                 \const real_type f = f1 * f1 + f2 * f2;\n\
                 \vector_type df(2);\n\
                 \df(0) = -4.0 * 10.0 * f1 * x(0) - 2.0 * f2;\n\
                 \df(1) = -2.0 * 10.0 * f1;\n\
                 \return std::make_pair(f, df);\n"
rosenbrock = Function "rosenbrock" 2 2 "{-1.2, 1}" "{1.0, 1.0}" 0.0 rosenbrockBody

functionList = [rosenbrock]

replaceTokens :: String -> Function -> (String, String)
replaceTokens s f = (filename, contents) where
                    filename = name f
                    contents = s

main :: IO ()
main = do
    let results = map (replaceTokens templateText) functionList 
--    putStrLn show results
    putStrLn templateText