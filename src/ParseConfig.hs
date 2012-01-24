{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module ParseConfig where
import Data.Monoid
import Utils
import Wagner
import Data.List
import System.Environment

type Key = String
type Value = String
type Table = [(Key, Value)]

readUpdateTable = [ ("sa",sa)
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
                     , logFlag         :: Bool                       
                     , printFlag       :: Bool                       
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
                  , ("logFlag"         , "False")--Set true to /suppress/
                  , ("logPrint"        , "False")                    
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
         , logFlag         = read $ extract "logFlag" table
         , printFlag       = read $ extract "printFlag" table        
         }
  
  where fname = extract "method" table

fromJust (Just x) = x

extract ::String -> Table -> String
extract keyString table = case lookup keyString table of 
  (Just val) -> val
  otherwise  -> extract keyString defaultTable
  
-- pesco's really cheap and simple flags and options (tm)
clParts = getArgs >>= return . (\(a,b) -> (a,drop 1 b)) . break (=="--")
getArgs' = clParts >>= \(a,b)-> return ([h:t| h:t<-a, h/='-' || null t] ++ b)
getFlags = clParts >>= \(a,_)-> return (concat [t| '-':t <- a])
getFlag x = getFlags >>= return . elem x
--ends here

parseConfig :: IO Config
parseConfig = do 
  args <- getArgs'
  let configFile = head args
  flagless <- fmap (configFromTable . tableFromFile configFile) $ readFile configFile
  noLogging <- getFlag 'l'
  noPrinting <- getFlag 'p'  
  return flagless{logFlag = noLogging, printFlag = noPrinting}
    
stripComments :: [String] -> [String]
stripComments = filter (not . null) . map (takeWhile (/='#'))
        
tableFromFile :: FilePath -> String -> [(String,String)]
tableFromFile fp = ([("configFile", fp)] ++) . keyValTable . stripComments . lines

keyValTable :: [String] -> [(String, String)]
keyValTable lines = zip keys vals
  where [keys, vals] = transpose $ map (map trim . split ':') lines