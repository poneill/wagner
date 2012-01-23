module Main where
import Wagner                         
import Utils
import ParseConfig 
import WriteOutput
import System (getArgs)
import System.Environment
import System.Time

-- pesco's really cheap and simple flags and options (tm)
clParts = getArgs >>= return . (\(a,b) -> (a,drop 1 b)) . break (=="--")
getArgs' = clParts >>= \(a,b)-> return ([h:t| h:t<-a, h/='-' || null t] ++ b)
getFlags = clParts >>= \(a,_)-> return (concat [t| '-':t <- a])
getFlag x = getFlags >>= return . elem x

main = do args <- getArgs'
          tick <- getClockTime
          let configFile = head args
          config <- parseConfig configFile
          let f = method config
          let iterations = numIterations config
          let converges = convergence config
          seqs <- readSequences $ dataFile config              
          mis <- seedMotifs seqs
          --Input ends here
          let g = Gestalt seqs mis
          g' <- if converges
                then converge g f
                else iterateN' iterations (>>= f) (return g)
          --Output begins
          tock <- getClockTime
          let time = (tick,tock)
          quiet <- getFlag 'q'
          writeOutput g' config time quiet
