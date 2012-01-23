module Main where
import Wagner                         
import Utils
import ParseConfig 
import WriteOutput
import System (getArgs)
import System.Environment
import System.Time

-- pesco's really cheap and simple flags and options (tm)
clparts = getArgs >>= return . (\(a,b) -> (a,drop 1 b)) . break (=="--")
getargs = clparts >>= \(a,b)-> return ([h:t| h:t<-a, h/='-' || null t] ++ b)
getflags = clparts >>= \(a,_)-> return (concat [t| '-':t <- a])
getflag x = getflags >>= return . elem x

main = do args <- getArgs
          tick <- getClockTime
          let configFile = head args
          config <- parseConfig configFile
          printing <- getflag 'p'
          let f = method config
          let fname = methodName config
          let iterations = numIterations config
          let converges = convergence config
          let iterLog = logIterations config
          let misLog = logMotifIndices config              
          let motifsLog = logMotifs config
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
          writeOutput g' config time 
          print fname
          print iterations
          print converges
          print $ gestaltEntropy g'
          print $ motifIndices g'
