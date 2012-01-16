module Main where
import Wagner                         
import Utils

main = do seqs <- readSequences "data/lexA_e_coli_120.csv"
          mis <- seedMotifs seqs
          let g = Gestalt seqs mis
          g' <- iterateN 1000 (>>= sa) (return g)
          print $ gestaltEntropy g'
