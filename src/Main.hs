module Main where
import Wagner                         
import Utils
import System (getArgs)

main = do args <- getArgs
          let iterations = read $ head args
          let methodName = args !! 1
          let f = case methodName of
                "sa" -> sa
                "patrify" -> patrify
                "patrifySweep" -> patrifySweep
          seqs <- readSequences "data/lexA_e_coli_120.csv"
          mis <- seedMotifs seqs
          let g = Gestalt seqs mis
          g' <- iterateN' iterations (>>= f) (return g)
          print $ gestaltEntropy g'
          print $ motifIndices g'
