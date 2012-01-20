module Main where
import Wagner                         
import Utils
import System (getArgs)

main = do args <- getArgs
          let iterations = read $ head args
          let functionName = args !! 1
          let f = interpretFunction functionName
          seqs <- readSequences "../data/lexA_e_coli_120.csv"
          mis <- seedMotifs seqs
          let g = Gestalt seqs mis
          g' <- if iterations > 0
                then iterateN' iterations (>>= f) (return g)
                else converge g f

          print $ gestaltEntropy g'
          print $ motifIndices g'
          print $ recoverMotifs g'


interpretFunction :: String -> (Gestalt -> IO Gestalt)          
interpretFunction fname = case fname of
  "sa" -> sa
  "patrify" -> patrify
  "patrifySweep" -> patrifySweep
  "greedy" -> greedy