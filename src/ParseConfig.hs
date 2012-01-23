{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module ParseConfig (parseConfig) where
import Data.Monoid
import Data.String.Utils
import Wagner
import Data.List

type Key = String
type Value = String
type Table = [(Key, Value)]

readUpdateTable = [ ("greedy",greedy)
                  , ("sa",sa)
                  , ("patrify",patrify)
                  ]

data Config = Config { dataFile        :: String
                     , methodName      :: Update
                     , convergence     :: Bool
                     , iterations      :: Int
                     , logIterations   :: Bool
                     , logMotifIndices :: Bool
                     , logMotifs       :: Bool
                     }

                 
defaultTable =   [ ("dataFile"        , "../data/lexA_e_coli_120.csv")
                  , ("methodName"      , "greedy")
                  , ("convergence"     , "True")
                  , ("iterations"      , "0")
                  , ("logIterations"   , "True")
                  , ("logMotifIndices" , "True")
                  , ("logMotifs"       , "True")
                  ]

configFromTable :: Table -> Config
configFromTable table = 
  Config { dataFile        = extract "dataFile" table
         , methodName      = fromJust $ lookup fname readUpdateTable
         , convergence     = read $ extract "convergence" table
         , iterations      = read $ extract "iterations" table
         , logIterations   = read $ extract "logIterations" table
         , logMotifIndices = read $ extract "logMotifIndices" table
         , logMotifs       = read $ extract "logMotifs" table
         }
  where fname = extract "methodName" table
fromJust (Just x) = x

extract ::String -> Table -> String
extract keyString table = case lookup keyString table of 
  (Just val) -> val
  otherwise  -> extract keyString defaultTable
    
parseConfig :: FilePath -> IO Config
parseConfig fp = fmap (configFromTable . tableFromFile) $ readFile fp

stripComments :: [String] -> [String]
stripComments = filter (not . null) . map (takeWhile (/='#'))
        
--tableFromConfig :: String -> [(String,String)]
tableFromFile = keyValTable . stripComments . lines

-- addDefaults :: Table -> Table
-- addDefaults table = table ++ [(k,v) | (k,v) <- defaultConfig, 
--                               not $ k `elem` map fst table]

keyValTable :: [String] -> [(String, String)]
keyValTable lines = zip keys vals
  where [keys, vals] = transpose $ map (map strip . split ":") lines
