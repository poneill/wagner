{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module ParseConfig where
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
                  , ("greedySweep",greedySweep)                    
                  , ("patrifySweep",patrifySweep)                    
                  ]

data Config = Config { dataFile        :: FilePath
                     , configFile      :: FilePath
                     , methodName      :: String
                     , method          :: Update                       
                     , convergence     :: Bool
                     , numIterations   :: Int
                     , logIterations   :: Bool
                     , logMotifIndices :: Bool
                     , logMotifs       :: Bool
                     }

                 
defaultTable =   [ ("dataFile"         , "../data/lexA_e_coli_120.csv")
                 , ("configFile"       , "notApplicable.wg")
                  , ("methodName"      , "greedy")
                  , ("method"          , "greedy")                    
                  , ("convergence"     , "True")
                  , ("numIterations"   , "0")
                  , ("logIterations"   , "True")
                  , ("logMotifIndices" , "True")
                  , ("logMotifs"       , "True")
                  ]

configFromTable :: Table -> Config
configFromTable table = 
  Config { dataFile        = extract "dataFile" table
         , configFile      = extract "configFile" table
         , methodName      = fname
         , method          = fromJust $ lookup fname readUpdateTable
         , convergence     = read $ extract "convergence" table
         , numIterations   = read $ extract "iterations" table
         , logIterations   = read $ extract "logIterations" table
         , logMotifIndices = read $ extract "logMotifIndices" table
         , logMotifs       = read $ extract "logMotifs" table
         }
  where fname = extract "method" table

fromJust (Just x) = x

extract ::String -> Table -> String
extract keyString table = case lookup keyString table of 
  (Just val) -> val
  otherwise  -> extract keyString defaultTable
    
parseConfig :: FilePath -> IO Config
parseConfig fp = fmap (configFromTable . tableFromFile fp) $ readFile fp

stripComments :: [String] -> [String]
stripComments = filter (not . null) . map (takeWhile (/='#'))
        
tableFromFile :: FilePath -> String -> [(String,String)]
tableFromFile fp = ([("configFile", fp)] ++) . keyValTable . stripComments . lines

keyValTable :: [String] -> [(String, String)]
keyValTable lines = zip keys vals
  where [keys, vals] = transpose $ map (map strip . split ":") lines
