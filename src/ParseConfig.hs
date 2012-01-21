module ParseConfig where

import Data.Monoid
import Data.String.Utils
import Wagner
import Data.List

data ConfigData = ConfigData { dataFile :: FilePath
                                , methodName :: String
                                , convergence :: Bool
                                , iterations :: Int
                                , logIterations :: Bool
                                , logMotifIndices :: Bool
                                , logMotifs :: Bool
                                }
                  deriving (Show)
                           
defaultConfig = ConfigData { dataFile = "../data/lexA_e_coli_120.csv"
                           , methodName = "greedy"
                           , convergence = True
                           , iterations = 0
                           , logIterations = True
                           , logMotifIndices = True
                           , logMotifs = True
                           }
  
--parseConfig :: FilePath -> IO ConfigData
--parseConfig fp = do keyVals <- fmap (tableFromConfig) $ readFile fp
--                    return defaultConfig { dataFile = 

stripComments :: [String] -> [String]
stripComments = (filter (not . null)) . (map (takeWhile (/='#')))
        
--tableFromConfig :: String -> [(String,String)]
--tableFromConfig fileContents = keyValTable . stripComments . lines


keyValTable lines = zipWith (,) keys vals
  where [keys, vals] = transpose $ map ((map strip) . (split ":")) lines

interpretFunction :: String -> (Gestalt -> IO Gestalt)          
interpretFunction fname = case fname of
  "sa" -> sa
  "patrify" -> patrify
  "patrifySweep" -> patrifySweep
  "greedy" -> greedy