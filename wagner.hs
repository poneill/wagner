import Data.List
import Data.Char (isSpace)
import System.IO

type DistanceMatrix = [[Int]]
type Sequence = String
type PSSM = [[Float]]
type Motif = [Sequence]

delta = "ACGT"
epsilon = 1e-10
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

readMotif :: FilePath -> IO Motif
readMotif filePath = do
  content <- readFile filePath
  let ls = map trim (lines content)
  return (filter ((/= '>') . head) ls)