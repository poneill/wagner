module WriteOutput where
import System.Time
import System.Environment
import Data.List
import ParseConfig
import Wagner 
import Data.String.Utils

writeOutput :: Gestalt -> Config -> (ClockTime,ClockTime) -> IO ()
writeOutput g config (tick,tock) = writeFile fp outputString 
  where cf = configFile config
        fp = makeFileName cf tick
        outputString = unlines [configString, timeString, motifString]
        configString = formatConfigFile $ cf
        timeString = formatTime tick tock
        motifString = formatMotifs g
        
makeFileName :: FilePath -> ClockTime -> String
makeFileName fp tick = fp' ++ tick' ++ ".log"
  where fp' = replace ".wg" "" fp
        tick' = replace " " "_" (show tick)
  
formatConfigFile :: FilePath -> String
formatConfigFile fp = "Simulation run with config file: " ++ fp

formatTime :: ClockTime -> ClockTime -> String
formatTime tick tock = unlines [timeBegan, timeFinished, timeTook]
  where timeBegan =    "Simulation began at: " ++ show tick
        timeFinished = "Simulation ended at: " ++ show tock
        timeTook = "Simulation took: " ++ timeDiffToString diff
        diff = diffClockTimes tock tick
  
formatMotifs :: Gestalt -> String
formatMotifs g@(Gestalt seqs mis) = unlines $ map unwords $ untuples
  where mis' = mmap (padLeft . show) mis
        motifs' = transpose $ recoverMotifs g
        tuples = zipWith zip mis' motifs'
        untuples  = mmap (\(x, y) -> x ++ " " ++ y) contents
        contents = zipWith zip mis' motifs'

mmap :: (a -> b) -> [[a]] -> [[b]]
mmap f = map (map f)

padLeft :: String -> String
padLeft str = (replicate (4 - length str) ' ') ++ str