--wagner.hs
import Data.List
import Data.Char (isSpace)
import System.IO
import System.Random hiding (split)
import Control.Monad (replicateM)
import Debug.Trace
type DistanceMatrix = [[Int]]
type Sequence = String
type NamedSequence = (Int,String)
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
type Index = Int --An index is a position in a sequence
type Name = Int --A name is a tag denoting a motif element.
type NamedPSSM = (Name, PSSM)
type VarMatrix = [[Float]]
data Gestalt = Gestalt { sequences :: Sequences 
                       , motifIndices :: MotifIndices
                       }
             deriving Show
  
delta = "ACGT"
epsilon = 1/100
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

--argMax :: (Ord b) => (a -> b) -> [a] -> a

argMax f = foldl1 (\x x' -> if f x' > f x then x' else x)

--argMin :: (Ord b) => (a -> b) -> [a] -> a

argMin f = foldl1 (\x x' -> if f x' <= f x then x' else x)
  
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

varianceMatrix :: (Floating a) => MotifIndices -> [[a]]
varianceMatrix mis = map (map (variance . map fromIntegral) . transpose) (transpose dms)
  where dms = [distanceMatrix i mis | i <- [0..length mis - 1]]
  
rescoreSequence :: Sequence -> Sequences -> MotifIndices -> MotifIndex
--Accepts a sequence and its LOO MotifIndices, returns a MotifIndex for sequence
--by greedily assigning tfs in sequential order.
rescoreSequence seq seqs mis = [maxResponseOverSeq pssm seq | pssm <- pssms]
  where pssms = map (`recoverPSSM` seqs) $ transpose mis

score :: PSSM -> Sequence -> Float -- assumes sequence is as long as pssm
score pssm seq = sum $ zipWith (\p s -> p !! indexOf s) pssm seq

scoreAt :: PSSM -> Sequence -> Index -> Float
scoreAt pssm seq i = score pssm (drop i seq)

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
                           

{- in Ivanization, we update the placements of the motif indices for a
given sequence.  We do this by iteratively adding motif indices,
sampled randomly according to their maxResponseOverPSSM and their
z-score.-}

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
                               folder ma b = ma >>= \x -> addToMIs seq x b

potential :: Sequence -> NamedPSSM -> Index -> MotifIndex -> VarMatrix -> Float
--potential can't be larger than 700, or exp (-potential) will underflow
potential seq (i,pssm) pos mi varMatrix = (bindingEnergy + a * stringEnergy)/700
  where bindingEnergy =printBE $ - (scoreAt pssm seq pos) --bigger is worse
        stringEnergy =printSE $ sum [log $ epsilon + energyFromString j jpos
                           | (j, jpos) <- zip [0..] mi, j /= i]
        energyFromString j jpos =printEFS $ 1/ (epsilon + var j) * fromIntegral (pos - jpos) ** 2
        var j = varMatrix !! i !! j
        a = 1/10

-- patrifyIthSequence :: Gestalt -> Int -> IO Gestalt
-- patrifyIthSequence g i = return seqs mis''
--   where mis = motifIndices g
--         seqs = sequences g
--         (mi,mis') = separate i mis
--         (seq,seqs') = separate i seqs
--         pssms' = recoverPSSMs (Gestalt seqs' mis')
--         dm = distanceMatrix i mis'
--         mis'' = patrifySequence seq 
--         varMatrix = varianceMatrix mis
        
patrify :: Gestalt -> IO Gestalt
patrify (Gestalt seqs mis) = do
  seqNum <- randomRIO (0, length seqs - 1)
  let seq = seqs !! seqNum      
  let mi = mis !! seqNum      
  motifNum <- randomRIO (0, numMotifs - 1)
  let looPSSM = recoverNthPSSM (delete seq seqs) (delete mi mis) motifNum -- revise
  i' <- assignIthIndex (seqNum,seq) (motifNum, looPSSM) mis
  let  mi' = replaceAt motifNum i' mi
  let mis' = replaceAt seqNum mi' mis
  return (Gestalt seqs mis')


assignIthIndex :: NamedSequence -> NamedPSSM -> MotifIndices -> IO Index

assignIthIndex (seqNum,seq) (i,pssm) mis =
  sample positions (\pos ->printEnergy $ exp (- energy pos))
  where end = length seq - length pssm --check this
        positions = [0..end]
        mi = mis !! seqNum
        energy pos =printPotential $ potential seq (i,pssm) pos mi varMatrix
        varMatrix = varianceMatrix (delete mi mis)
    
toMotifIndex :: [NamedPSSM] -> MotifIndex
toMotifIndex = map fst . sortWith snd
  
chooseRandomly :: [a] -> IO a
chooseRandomly xs = do
  r <- randomRIO (0, length xs - 1)
  return (xs !! r)
  
sortWith :: (Ord b) => (a -> b) -> [a] -> [a]
sortWith f xs = map fst $ sortBy g $ map (\x -> (x, f x)) xs
  where g (x, fx) (y, fy) = compare fx fy
          
orderMotifs :: [PSSM] -> Sequence -> Sequences -> IO [NamedPSSM]
-- establish an order in which the PSSMs are to be indexed.  For now,
-- they are just sorted by their max response over sequence
orderMotifs pssms seq seqs = return sorteds 
  where sorteds = sortBy f indexedPSSMs
        indexedPSSMs = zip [0..] pssms
        f p q 
          | maxOverSequence (snd p) seq < maxOverSequence (snd q) seq = LT
          | otherwise                                                 = GT

orderMotifs' :: [PSSM] -> Sequence -> Sequences -> IO [(Int, PSSM)]
orderMotifs' pssms seq seqs = orderBySampling indexedPSSMs f 
  where f p = maxOverSequence (snd p) seq
        indexedPSSMs = zip [0..] pssms
                                                 
--orderBySampling :: (Random b, Ord b, Floating b) => [a] -> (a -> b) -> IO [a]
--orderBySampling [] f = return []
orderBySampling [a] f = return [a]
orderBySampling as f = do { a <- sample as f
                          ; let aless = delete a as
                          ; aless' <- orderBySampling aless f
                          ; return (a : aless')
                          }
                            

--sample :: (Random b, Ord b, Floating b) => [a] -> (a -> b) -> IO a 
sample as f = do { r <- randomRIO (0.0,1.0)
                 ; return (sample' as f r)
                 }

--sample' :: (Ord b, Floating b) => [a] -> (a -> b) -> b -> a
-- Pick an a according to a likelihood function (and an implicit
-- constant k)
sample' as f r = fst $ argMin snd $ filter ((>= r) . snd)  tups
              where k = 1
                    faks =printFaks $  map (\a -> f a ** k) as
                    tups =printTubs $ zip as (scanl1 (+) (map (/z) faks))
                    z =printZ $ sum faks
                    
addToMIs :: Sequence -> [(Index,Index)] -> NamedPSSM -> IO [(Index,Index)]
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

recoverNthPSSM :: Sequences -> MotifIndices -> Int -> PSSM
recoverNthPSSM seqs mis n = recoverPSSM mic seqs 
  where mic = (transpose mis) !! n

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

updateSweep :: Gestalt -> Gestalt
updateSweep g = foldl updateIthSequence g is
  where is = (range . length . motifIndices) g

ivanSweep :: Gestalt -> IO Gestalt
--ivanSweep g | trace ("ivanSweep"++ " " ++ show (motifIndices g)) False = undefined
ivanSweep g = foldl (\mg i -> mg >>= \g -> ivanizeIthSequence g i) (return g) is
  where is = (range . length . motifIndices) g
        
--springConstant :: MotifIndices -> Index -> Index -> 
springConstant mis i j k = 1 / (epsilon + variance (map fromIntegral $ zipWith (-) is js))
  where is = selectColumn mis' i
        js = selectColumn mis' j
        mis' = removeNth mis k 

selectColumn :: [[a]] -> Index -> [a]
selectColumn xss i = [xs !! i | xs <- xss]

variance :: (Floating a) => [a] -> a
variance xs = mean (map (**2) xs) - mean xs ** 2

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

-- debugging

-- printPotential x | trace ("printPotential"++ " " ++ show x) False = undefined
-- printTubs xs | trace ("printTubs"++ " " ++ show xs) False = undefined
-- printSE x | trace ("printSE"++ " " ++ show x) False = undefined
-- printFaks xs | trace ("printFaks"++ " " ++ show xs) False = undefined
-- printZ xs | trace ("printZ"++ " " ++ show xs) False = undefined
-- printEnergy x | trace ("printEnergy"++ " " ++ show x) False = undefined
-- printEFS x | trace ("printEFS"++ " " ++ show x) False = undefined
-- printBE x | trace ("printBE"++ " " ++ show x) False = undefined

printBE x = x
printEFS x = x
printSE x = x
printEnergy x = x
printPotential x = x
printTubs xs = xs
printFaks xs = xs
printZ xs = xs
