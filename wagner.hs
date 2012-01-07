--wagner.hs
import Data.List
import Data.Char (isSpace)
import System.IO
import System.Random hiding (split)
import Control.Monad (replicateM)
import Debug.Trace
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
type MotifIndexCol = [Index]
type MotifIndices = [MotifIndex]
type Index = Int
data Gestalt = Gestalt { sequences :: Sequences 
                       , motifIndices :: MotifIndices
                       }
             deriving Show
  
delta = "ACGT"
epsilon = 1e-10
numMotifs = 3
motifLength = 16
uniformProbs = replicate 4 0.25
indexOf :: Char -> Index
indexOf base = unpack $ lookup base (zip delta [0..3])
  where unpack (Just x) = x

range :: (Integral a) => a -> [a]
range n = [0..(n-1)]

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
    where nthIndices = mis !! n

rescoreSequence :: Sequence -> Sequences -> MotifIndices -> MotifIndex
--Accepts a sequence and its LOO MotifIndices, returns a MotifIndex for sequence
--by greedily assigning tfs in sequential order.
rescoreSequence seq seqs mis = [maxResponseOverSeq pssm seq | pssm <- pssms]
  where pssms = map (`recoverPSSM` seqs) $ transpose mis

score :: PSSM -> Sequence -> Float
score pssm seq = sum $ zipWith (\p s -> p !! indexOf s) pssm seq

maxResponseOverSeq :: PSSM -> Sequence -> Index
maxResponseOverSeq pssm seq = head $ elemIndices (maximum scores) scores
  where scores = scoreSequence pssm seq

maxPSSMoverSeq :: [PSSM] -> Sequence -> PSSM
maxPSSMoverSeq pssms seq = argMax (`maxOverSequence` seq) pssms

scoreSequence :: PSSM -> Sequence -> [Float] --scan PSSM over sequence
scoreSequence pssm seq = map (score pssm) longEnoughs
  where longEnoughs = takeWhile (\tail -> length tail >= m) (tails seq)
        m = length pssm

updateAlignment :: Gestalt -> IO Gestalt
updateAlignment gestalt = do { let seqs = sequences gestalt
                             ; let mis = motifIndices gestalt 
                             ; i <- randomRIO (0, length mis - 1)
                             ; return (updateIthSequence gestalt i)
                              }
                           

ivanizeIthSequence :: Gestalt -> Int -> IO Gestalt
ivanizeIthSequence g i = do { motifOrder <- orderMotifs pssms' seq seqs'
                            ; placements <- foldl folder (return []) motifOrder
                            ; let mi' = collocateMotifIndex placements
                            ; return (Gestalt seqs (insertAt i mi' mis'))
                            }
                         where mis = motifIndices g
                               seqs = sequences g
                               (_,mis') = separate i mis
                               (seq,seqs') = separate i seqs
                               pssms' = recoverPSSMs (Gestalt seqs' mis')
                               folder = \ma b -> ma >>= \x -> addToMIs seq x b
                               
orderMotifs :: [PSSM] -> Sequence -> Sequences -> IO [(Int, PSSM)]
-- establish an order in which the PSSMs are to be indexed.  For now,
-- they are just sorted by their max response over sequence
orderMotifs pssms seq seqs = return sorteds 
  where sorteds = sortBy f indexedPSSMs
        indexedPSSMs = zip [0..] pssms
        f p q 
          | maxOverSequence (snd p) seq < maxOverSequence (snd q) seq = LT
          | otherwise = GT


addToMIs :: Sequence -> [(Index,Index)] -> (Int,PSSM) -> IO [(Index,Index)]
-- [(i,j)] denotes the placement index j of the ith motif
addToMIs seq ijs (i,pssm) = return (ijs ++ [(i,j)])
  where j = maxResponseOverSeq pssm seq

collocateMotifIndex :: [(Index,Index)] -> MotifIndex
collocateMotifIndex = map snd . sort
  
separate :: Index -> [a] -> (a,[a])
separate i seqs = (seqs !! i, removeNth seqs i)

insertAt :: Index -> a -> [a] -> [a] 
insertAt i a as = take i as ++ [a] ++ drop i as 

replaceAt :: Index -> a -> [a] -> [a] 
replaceAt i a as = take i as ++ [a] ++ drop (i + 1) as 

updateIthSequence :: Gestalt -> Index -> Gestalt
updateIthSequence gestalt i = Gestalt seqs mis'
    where 
      seqs = sequences gestalt
      mis = motifIndices gestalt
      seq = seqs !! i
      seqsRest = removeNth seqs i
      mi = mis !! i
      misRest = removeNth mis i  
      mi' = rescoreSequence seq seqs misRest
      mis' = replaceAt i mi' mis 
            
maxOverSequence :: PSSM -> Sequence -> Float --scan PSSM over sequence, take max
maxOverSequence pssm seq = maximum  $ scoreSequence pssm seq

recoverPSSM :: MotifIndexCol -> Sequences -> PSSM
recoverPSSM mic seqs = makePSSM (recoverMotif mic seqs) uniformProbs

recoverPSSMs :: Gestalt -> [PSSM]
recoverPSSMs gestalt = map (`recoverPSSM` seqs) mics
  where mics = transpose $ motifIndices gestalt
        seqs = sequences gestalt

recoverMotif :: MotifIndexCol -> Sequences -> Motif
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

iterateN :: Int -> (a -> a) -> a -> a
iterateN n f x = iterate f x !! n

patSweep :: Gestalt -> Gestalt
patSweep g | trace ("patSweep"++ " " ++ show (motifIndices g)) False = undefined
patSweep g = foldl updateIthSequence g is
  where is = (range . length . motifIndices) g

ivanSweep :: Gestalt -> IO Gestalt
--ivanSweep g | trace ("ivanSweep"++ " " ++ show (motifIndices g)) False = undefined
ivanSweep g = foldl (\mg i -> mg >>= \g -> ivanizeIthSequence g i) (return g) is
  where is = (range . length . motifIndices) g
        
--springConstant :: MotifIndices -> Index -> Index -> 
springConstant mis i j k = 1 / variance $ map fromIntegral $ zipWith (-) is js
  where is = selectColumn mis' i
        js = selectColumn mis' j
        mis' = removeNth mis k 

selectColumn :: [[a]] -> Index -> [a]
selectColumn xss i = [xs !! i | xs <- xss]

variance :: (Floating a) => [a] -> a
variance xs = mean (map (**2) xs) - (mean xs) ** 2

mean :: (Fractional a) => [a] -> a
mean xs = sum xs / fromIntegral (length xs)

converge :: Gestalt -> IO Gestalt
converge g = converge' g (updateAlignment g)
  where converge' g mg = do { g' <- mg
                            ; if motifIndices g == motifIndices g'
                                 then return g
                              else converge' g' (updateAlignment g')
                            }
          
fixpoint :: (Eq b) => (a -> a) -> a -> (a -> b) -> a
fixpoint f a p = fst $ head $ dropWhile (\(x,y) -> p x /= p y) $ zip its (tail its)
  where its = iterate f a
  
main = do { seqs <- readSequences "data/lexA_e_coli_120.csv"
          ; mis <- seedMotifs seqs
          ; return (Gestalt seqs mis)
          }
