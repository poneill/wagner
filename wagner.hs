import Data.List
import Data.Char (isSpace)
import System.IO
import System.Random

type DistanceMatrix = [[Int]]
type Sequence = String
type PSSM = [[Float]]
type Motif = [Sequence]
type Sequences = [Sequence] --the idea here is that Sequences are just
                            --raw input, whereas a Motif is
                            --semantically significant alignment
type MotifIndex = [Int] --Records left endpoints of occurrence of
                        --motif in each sequence.  Can be used to
                        --recover Motif from Sequences


delta = "ACGT"
epsilon = 1e-10
numMotifs = 3
motifLength = 6

log2 :: (Floating a) => a -> a
log2 = logBase 2

trim      :: String -> String --stole this from wikipedia for portability
trim      = f . f
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
seedMotif seqs = sequence $ [randomRIO (0,length seq) | seq <- seqs]

seedMotifs :: Sequences -> IO [MotifIndex]
seedMotifs seqs = sequence $ replicate numMotifs (seedMotif seqs)

distanceMatrix :: Int -> [MotifIndex] -> DistanceMatrix
--Return a distance matrix for the nth sequence
distanceMatrix n mis = [[i - j | i <- nthIndices] | j <- nthIndices]
    where nthIndices = (transpose mis) !! n

--rescoreSequence :: Sequence -> [MotifIndex] -> [MotifIndex]

  
recoverMotif :: MotifIndex -> Sequences -> Motif
recoverMotif mi seqs = zipWith (\m s -> (take motifLength . drop m) s) mi seqs

readSequences :: FilePath -> IO Sequences
readSequences filePath = do
  content <- readFile filePath
  let ls = map trim (lines content)
  return (filter ((/= '>') . head) ls)