module Main where
import Wagner                         
import Utils
import ParseConfig 
import WriteOutput
import System (getArgs)
import System.Environment
import System.Time

main = do tick <- getClockTime
          config <- parseConfig
          let f = method config
          let iterations = numIterations config
          let converges = convergence config
          seqs <- readSequences $ dataFile config              
          mis <- seedMotifs seqs
          --Input ends here
          let g = Gestalt seqs mis
          g' <- if converges
                then convergeCyclic g f iterations
                else iterateN' iterations (>>= f) (return g)
          --Output begins
          tock <- getClockTime
          let time = (tick,tock)
          writeOutput g' config time
