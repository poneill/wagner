module Main where
import Wagner                         
import Utils
import ParseConfig 
import System (getArgs)

main = do args <- getArgs
          let configFile = head args
          config <- parseConfig configFile
          let f = (method config)
          let fname = (methodName config)              
          seqs <- readSequences (dataFile config)
          let iterations = (numIterations config) 
          let converges = (convergence config)
          let iterLog = logIterations config
          let misLog = logMotifIndices config              
          let motifsLog = logMotifs config                            
          mis <- seedMotifs seqs
          let g = Gestalt seqs mis
          g' <- if converges
                then converge g f
                else iterateN' iterations (>>= f) (return g)
          print $ fname
          print iterations
          print converges
          print $ gestaltEntropy g'
          print $ motifIndices g'
          print $ recoverMotifs g'