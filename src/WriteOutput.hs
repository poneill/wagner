module WriteOutput where
import System.Time
import System.Environment
import Data.List
import ParseConfig
import Wagner 
import Utils
import Control.Monad

writeOutput :: Gestalt -> Config -> (ClockTime,ClockTime) -> IO ()
writeOutput g config (tick,tock) = do 
  let outputString = prepareOutput g config (tick,tock)
  let cf = configFile config
  let fp = makeFileName cf tick
  let suppressPrinting = printFlag config
  let suppressLogging = logFlag config
  unless (not suppressPrinting) $ putStrLn outputString
  unless (not suppressLogging) $ writeFile fp outputString 

prepareOutput :: Gestalt -> Config -> (ClockTime,ClockTime) -> String
prepareOutput g config (tick,tock) = outputMessage
  where cf = configFile config
        outputMessage = unlines [configLine, timeLine, motifLine, entropyLine, methodLine]
        configLine = formatConfigFile cf
        timeLine = formatTime tick tock
        motifLine = formatMotifs g
        methodLine = formatMethod config
        entropyLine = formatEntropy g
          
formatEntropy :: Gestalt -> String
formatEntropy g = unlines [total, perMotif]
  where total = "Total Entropy: " ++ show (gestaltEntropy g)
        perMotif = "Entropy by Motif: " ++ unwords entropyStrings
        entropyStrings = map show (gestaltEntropies g)
        
formatMethod :: Config -> String 
formatMethod config = unlines [fnameString, howLong]
  where fnameString = "Simulation ran with method: " ++ fname
        howLong = if converges then "until convergence" else iterString
        iterString = "for " ++ show iterations ++ " iterations" 
        iterations = numIterations config
        converges = convergence config
        fname = methodName config

makeFileName :: FilePath -> ClockTime -> String
makeFileName fp tick = process fp
  where process = addPath . addExt . addTimeStamp . removePath . removeExt
        removePath = replace "../config/" ""
        removeExt = replace ".wg" ""
        fpStripped = (removePath . removeExt) fp
        addPath = ("../log/" ++)
        addExt = (++ ".log")
        addTimeStamp = (++ tick')
        tick' = '_' : replace " " "_" (show tick)
  
formatConfigFile :: FilePath -> String
formatConfigFile fp = "Simulation run with config file: " ++ fp

formatTime :: ClockTime -> ClockTime -> String
formatTime tick tock = unlines [timeBegan, timeFinished, timeTook]
  where timeBegan =    "Simulation began at: " ++ show tick
        timeFinished = "Simulation ended at: " ++ show tock
        timeTook = "Simulation took: " ++ finalDiff
        finalDiff = if length stringDiff > 0 then stringDiff else "0 secs"
        stringDiff = timeDiffToString diff
        diff = diffClockTimes tock tick
  
formatMotifs :: Gestalt -> String
formatMotifs g@(Gestalt seqs mis) = unlines $ map unwords untuples
  where mis' = mmap (padLeft . show) mis
        motifs' = transpose $ recoverMotifs g
        tuples = zipWith zip mis' motifs'
        untuples  = mmap (\(x, y) -> x ++ " " ++ y) contents
        contents = zipWith zip mis' motifs'

mmap :: (a -> b) -> [[a]] -> [[b]]
mmap f = map (map f)

padLeft :: String -> String
padLeft str = replicate (4 - length str) ' ' ++ str