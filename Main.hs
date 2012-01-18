module Main where
import Wagner                         
import Utils
import System (getArgs)

main = do args <- getArgs
          let iterations = read $ head args
          seqs <- readSequences "data/lexA_e_coli_120.csv"
          mis <- seedMotifs seqs
          let g = Gestalt seqs mis
          g' <- iterateN iterations (>>= sa) (return g)
          print $ gestaltEntropy g'
          print $ motifIndices g'
