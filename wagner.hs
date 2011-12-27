import Data.List
import Data.Char (isSpace)
import System.IO
import System.Random hiding (split)
import Control.Monad (replicateM)

type DistanceMatrix = [[Int]]
type Sequence = String
type PSSM = [[Float]]
type Motif = [Sequence]
type Sequences = [Sequence] --the idea here is that Sequences are just
                            --raw input, whereas a Motif is
                            --semantically significant alignment
type MotifIndex = [Index] --Records left endpoints of occurrence of
                        --motif in each sequence.  Can be used to
                        --recover Motif from Sequences
type MotifIndices = [MotifIndex]
type Index = Int
data Gestalt = Gestalt { sequences :: Sequences 
                       , motifIndices :: MotifIndices
                       }
             deriving Show
  
delta = "ACGT"
epsilon = 1e-10
numMotifs = 1
motifLength = 6
uniformProbs = replicate 4 0.25
indexOf :: Char -> Index
indexOf base = unpack $ lookup base (zip delta [0..3])
  where unpack (Just x) = x

log2 :: (Floating a) => a -> a
log2 = logBase 2

argMax :: (Ord b) => (a -> b) -> [a] -> a
argMax f = foldl1 (\x x' -> if f x' > f x then x' else x) 
  
trim :: String -> String --stole this from wikipedia for portability
trim = f . f
  where f = reverse . dropWhile isSpace

columnProbs :: Sequence -> [Float]
columnProbs column = [epsilon + fromIntegral (numBases base column) / n 
                          | base <- delta]
    where n = fromIntegral $ length column
          numBases base = length . filter (== base)

makePSSM :: Motif -> [Float] -> PSSM
makePSSM motif bgProbs = map (f . columnProbs) columns 
    where columns = transpose motif
          f column = zipWith (\c bg -> log2 (c / bg)) column bgProbs

seedMotif :: Sequences -> IO MotifIndex
seedMotif seqs = sequence [randomRIO (0,length seq) | seq <- seqs]

seedMotifs :: Sequences -> IO MotifIndices
seedMotifs  = fmap transpose . replicateM numMotifs . seedMotif

distanceMatrix :: Int -> MotifIndices -> DistanceMatrix
--Return a distance matrix for the nth sequence
distanceMatrix n mis = [[abs (i - j) | i <- nthIndices] | j <- nthIndices]
    where nthIndices = transpose mis !! n

rescoreSequence :: Sequence -> Sequences -> MotifIndices -> MotifIndex
--Accepts a sequence and its LOO MotifIndices, returns a MotifIndex for sequence
rescoreSequence seq seqs mis = [maxResponseOverSeq pssm seq]
  where pssm = maxPSSMoverSeq pssms seq
        pssms = map (`recoverPSSM` seqs) mis

score :: PSSM -> Sequence -> Float
score pssm seq = sum $ zipWith (\p s -> p !! indexOf s) pssm seq

maxResponseOverSeq :: PSSM -> Sequence -> Index
maxResponseOverSeq pssm seq = argMax (\n -> maxOverSequence pssm (drop n seq)) [0..]

maxPSSMoverSeq :: [PSSM] -> Sequence -> PSSM
maxPSSMoverSeq pssms seq = argMax (`maxOverSequence` seq) pssms

scoreSequence :: PSSM -> Sequence -> [Float] --scan PSSM over sequence
scoreSequence pssm seq = map (score pssm) longEnoughs
  where longEnoughs = takeWhile (\tail -> length tail >= m) (tails seq)
        m = length pssm

updateAlignment :: Gestalt -> IO Gestalt
updateAlignment gestalt = do { let seqs = sequences gestalt
                             ; let mis = motifIndices gestalt 
                             ; i <- randomRIO (0, length mis)
                             ; return (updateIthSequence gestalt i)
                              }
                           

updateIthSequence :: Gestalt -> Int -> Gestalt
updateIthSequence gestalt i = Gestalt seqs mis'
    where 
      seqs = sequences gestalt
      mis = motifIndices gestalt
      seq = seqs !! i
      seqsRest = removeNth seqs i
      mi = mis !! i
      misRest = removeNth mis i  
      mi' = rescoreSequence seq seqs misRest
      mis' = take (i - 1) mis ++ [mi'] ++ drop i mis
            
maxOverSequence :: PSSM -> Sequence -> Float --scan PSSM over sequence, take max
maxOverSequence pssm seq = maximum  $ scoreSequence pssm seq

recoverPSSM :: MotifIndex -> Sequences -> PSSM
recoverPSSM mi seqs = makePSSM (recoverMotif mi seqs) uniformProbs

recoverPSSMs :: Gestalt -> [PSSM]
recoverPSSMs gestalt = map (`recoverPSSM` seqs) mis
  where mis = motifIndices gestalt
        seqs = sequences gestalt

recoverMotif :: MotifIndex -> Sequences -> Motif
recoverMotif = zipWith (\m s -> (take motifLength . drop m) s)

selectSequence :: Sequences -> IO (Sequence, Sequences)
selectSequence seqs = do {
  i <- randomRIO (0, length seqs);
  return (seqs!!i, removeNth seqs i)
  }
  
readSequences :: FilePath -> IO Sequences
readSequences filePath = do
  content <- readFile filePath
  return (sanitizeFASTA content)
  
sanitizeFASTA :: String -> Sequences
sanitizeFASTA content = map (filter (/= ',')) relevantLines
  where ls = map trim (lines content)
        relevantLines = filter ((/= '>') . head) ls
        
removeNth :: [a] -> Int -> [a]
removeNth xs n = ys ++ tail zs
  where (ys,zs) = splitAt n xs   

main = do { seqs <- readSequences "data/lexA_e_coli_120.csv"
          ; mis <- seedMotifs seqs
          ; let gestalt = Gestalt seqs mis 
          }
