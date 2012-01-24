--wagner.hs
module Wagner where
import Prelude
import Data.List
import Data.Char (isSpace)
import System.IO
import System.Random hiding (split)
import Control.Monad (replicateM)
import Debug.Trace
import Utils
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
type Indices = [Index]
type Name = Int --A name is a tag denoting a motif element.
type NamedPSSM = (Name, PSSM)
type Assigner = NamedSequence -> NamedPSSM -> MotifIndices -> StatMatrices -> IO Index
type VarMatrix = [[Float]]
type MeanMatrix = [[Float]]
type StatMatrices = (MeanMatrix, VarMatrix)
type Update = (Gestalt -> IO Gestalt)
data Gestalt = Gestalt { sequences :: Sequences 
                       , motifIndices :: MotifIndices
                       }
             deriving Show
  

delta = "ACGT"
epsilon = 1/100
numMotifs = 5
motifLength = 6
uniformProbs = replicate 4 0.25

indexOf :: Char -> Index
indexOf 'A' = 0
indexOf 'C' = 1
indexOf 'G' = 2
indexOf 'T' = 3

spacer = motifLength
spacerPenalty = 100

assignSpacerPenalty :: Index -> Index -> Float
assignSpacerPenalty i j
  | abs (i - j) * 2 < spacer = spacerPenalty
  | otherwise = 1

gestaltEntropies :: Gestalt -> [Float]
gestaltEntropies g = map motifEntropy motifs
  where motifs = recoverMotifs g
gestaltEntropy :: Gestalt -> Float
gestaltEntropy g = sum $ gestaltEntropies g
        
motifEntropy :: Motif -> Float
motifEntropy motif = sum $ map colEntropy motif'
  where motif' = transpose motif

colEntropy :: (Ord a) => [a] -> Float
colEntropy col = (-1) * sum (map (\x -> x * log2 (x + epsilon)) freqs)
  where freqs =  [fromIntegral count / fromIntegral len | count <- counts]
        counts = getCounts col
        len = length col

entropyEpsilon = 10*10e-10
entropy xs = (- 1) * sum (map (\x -> x * log2 (x + entropyEpsilon)) xs)

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
seedMotif seqs = sequence [randomRIO (0,length seq - motifLength) | seq <- seqs]

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

bindingEnergyAt :: PSSM -> Sequence -> Index -> Float --lower is better
bindingEnergyAt pssm seq i = scoreToEnergy score
  where score = scoreAt pssm seq i

--reverseSigmoid x = 1 / (1 + exp (x))
scoreToEnergy x = - x

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
sampled randomly according to their maxResponseOverSeq and their
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

potential :: Sequence -> NamedPSSM -> Index -> MotifIndex -> MotifIndices -> StatMatrices -> Float
--potential can't be larger than 700, or exp (-potential) will underflow
--higher potential means lower probability state
potential seq (i,pssm) pos mi mis statMatrices = bE + a * sE
  where bE = bindingEnergyAt pssm seq pos --bigger is worse
        sE = springEnergy seq (i,pssm) pos mi mis statMatrices
        a = 0.1
        
potentials :: Sequence -> NamedPSSM -> Indices -> MotifIndex -> MotifIndices -> StatMatrices -> [Float]
--potential can't be larger than 700, or exp (-potential) will underflow
--higher potential means lower probability state
potentials seq (i,pssm) ps mi mis statMatrices = map (\p -> bE p + a * sE p) ps
  where bE = bindingEnergyAt pssm seq --bigger is worse
        sE p= springEnergy seq (i,pssm) p mi mis statMatrices
        a = 0.1

springEnergy :: Sequence -> NamedPSSM -> Index -> MotifIndex -> MotifIndices -> StatMatrices -> Float
springEnergy seq (i,pssm) pos mi mis (muMatrix,varMatrix) = sum [log (epsilon + 
                                                       energyFromSpring j jpos) + 
                                                      assignSpacerPenalty pos jpos 
                                                     | (j, jpos) <- zip [0..] mi, j /= i]
  where energyFromSpring j jpos = displacement (i, pos) (j,jpos) **2 / (epsilon + var j) 
        var j = varMatrix !! i !! j
        mu i j = muMatrix !! i !! j
        displacement (i,ipos) (j,jpos) = fromIntegral (ipos - jpos) - mu i j

makeStatMatrices :: MotifIndices -> StatMatrices 
makeStatMatrices mis = (meanMatrix mis, varianceMatrix mis)

meanMatrix :: MotifIndices -> [[Float]] --compute resting lengths
--matrix is symmetric, upper triangular; could just compute half of it
meanMatrix mis = [[mean [(mi!!i) - (mi!!j) | mi <- mis]
                  |i <- motifRange] 
                 | j <- motifRange]
  where motifRange = [0..numMotifs - 1]
        
        
assignIthWrapper :: Assigner -> Gestalt -> Int -> IO Gestalt
assignIthWrapper assigner (Gestalt seqs mis) seqNum = do
  let (seq, seqs') = separate seqNum seqs      
  let (mi, mis') = separate seqNum mis      
  motifNum <- randomRIO (0, numMotifs - 1)
  let looPSSM = recoverNthPSSM seqs' mis' motifNum
  let statMatrices = makeStatMatrices mis'
  i' <- assigner (seqNum,seq) (motifNum, looPSSM) mis statMatrices
  let  mi' = replaceAt motifNum i' mi
  let mis' = replaceAt seqNum mi' mis
  return (Gestalt seqs mis')
  
patrifyIthSeq :: Gestalt -> Int -> IO Gestalt
patrifyIthSeq = assignIthWrapper assignIthIndex 

greedyIthSeq :: Gestalt -> Int -> IO Gestalt
greedyIthSeq = assignIthWrapper $ compose4 return greedyAssignIthIndex'

sweepWrapper :: (Gestalt -> Int -> IO Gestalt) -> Gestalt -> IO Gestalt
sweepWrapper updater g = foldl' f (return g) is
  where numSeqs = length $ motifIndices g
        is = [0..numSeqs - 1]
        f mg i = mg >>= \g -> updater g i

patrifySweep :: Gestalt -> IO Gestalt
patrifySweep = sweepWrapper patrifyIthSeq
        
greedySweep :: Gestalt -> IO Gestalt
greedySweep = sweepWrapper greedyIthSeq

-- sweepify :: (Gestalt -> IO Gestalt) -> (Gestalt -> IO Gestalt)
-- sweepify method = \g -> foldl' f (return g) [0..(length (motifIndices g)) - 1]
--   where f mg i = mg >>= \g -> method g i
        
-- saSweep :: Gestalt -> IO Gestalt
-- saSweep = sweepify sa
        
patrify :: Gestalt -> IO Gestalt
patrify g = do
  seqNum <- randomRIO (0, length (sequences g) - 1)
  patrifyIthSeq g seqNum

sa :: Gestalt -> IO Gestalt
sa g = do
  seqNum <- randomRIO (0, length (sequences g) - 1)
  saIth g seqNum

saIth :: Gestalt -> Int -> IO Gestalt
saIth = assignIthWrapper saCore

saCore :: Assigner 
saCore (seqNum,seq) (motifNum,looPSSM) mis statMatrices = do
  proPos <- assignIthIndex (seqNum,seq) (motifNum, looPSSM) mis statMatrices
  let mi = mis !! seqNum
  let curPos = mi !! motifNum
  let curPot = potential seq (motifNum,looPSSM) curPos mi mis statMatrices
  let proPot = potential seq (motifNum,looPSSM) proPos mi mis statMatrices
  acceptProposed <- accept' curPot proPot
  let nextPos = if acceptProposed then proPos else curPos
  return nextPos


accept' :: Float -> Float -> IO Bool --old 
accept' current proposed = do
  r <- randomRIO (0.0,1.0)
  -- we have already mapped energies into probabilities at this point,
  -- so take the ratio rather than exp (current - proposed)
  let acceptProposed = (proposed < current) || (r < current / proposed) 
  return acceptProposed
      
assignIthIndex :: NamedSequence -> NamedPSSM -> MotifIndices -> StatMatrices -> IO Index
assignIthIndex (seqNum,seq) (i,pssm) mis statMatrices = sample positions likelihood
  where end = length seq - length pssm --check this
        positions = [0..end]
        (mi, mis') = separate seqNum mis
        energy pos = potential seq (i,pssm) pos mi mis' statMatrices
        likelihood pos = exp (- energy pos) --via Boltzmann distribution

greedyAssignIthIndex' :: NamedSequence -> NamedPSSM -> MotifIndices -> StatMatrices -> Index
greedyAssignIthIndex' (seqNum,seq) (i,pssm) mis statMatrices = fst $ argMax snd (zip positions ps)
  where end = length seq - length pssm --check this
        (mi, mis') = separate seqNum mis
        positions = [0..end]
        ps = potentials seq (i,pssm) positions mi mis' statMatrices
        likelihoods = map (\p -> exp (- p)) ps --via Boltzmann distribution


toMotifIndex :: [NamedPSSM] -> MotifIndex
toMotifIndex = map fst . sortWith snd
    
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
                                                 
                    
addToMIs :: Sequence -> [(Index,Index)] -> NamedPSSM -> IO [(Index,Index)]
-- [(i,j)] denotes the placement index j of the ith motif
addToMIs seq ijs (i,pssm) = return (ijs ++ [(i,j)])
  where j = maxResponseOverSeq pssm seq

collocateMotifIndex :: [(Index,Index)] -> MotifIndex
collocateMotifIndex = map snd . sort
  
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
  where mic = transpose mis !! n

recoverPSSMs :: Gestalt -> [PSSM]
recoverPSSMs gestalt = map (`recoverPSSM` seqs) mics
  where mics = transpose $ motifIndices gestalt
        seqs = sequences gestalt

recoverMotif :: MotifIndexCol -> Sequences -> Motif
recoverMotif = zipWith (\m s -> (take motifLength . drop m) s)

recoverMotifs :: Gestalt -> [Motif]
recoverMotifs g = map (`recoverMotif` seqs) mics
  where mics = transpose $ motifIndices g
        seqs = sequences g
        

selectSequence :: Sequences -> IO (Sequence, Sequences)
selectSequence seqs = do {
  i <- randomRIO (0, length seqs);
  return (seqs!!i, removeNth seqs i)
  }
  

updateSweep :: Gestalt -> Gestalt
updateSweep g = foldl updateIthSequence g is
  where is = (range . length . motifIndices) g

ivanSweep :: Gestalt -> IO Gestalt
--ivanSweep g | trace ("ivanSweep"++ " " ++ show (motifIndices g)) False = undefined
ivanSweep g = foldl (\mg i -> mg >>= \g -> ivanizeIthSequence g i) (return g) is
  where is = (range . length . motifIndices) g
        
springConstant :: MotifIndices -> Index -> Index -> Int -> Float
springConstant mis i j k = 1 / (epsilon + variance (map fromIntegral $ zipWith (-) is js))
  where is = selectColumn mis' i
        js = selectColumn mis' j
        mis' = removeNth mis k 

converge :: Gestalt -> (Gestalt -> IO Gestalt) -> IO Gestalt
converge g f = converge' g (f g)
  where converge' g mg = do { g' <- mg
                            ; if motifIndices g == motifIndices g'
                                 then return g
                              else converge' g' (f g')
                            }
                         
--convergeCyclic :: Gestalt -> (Gestalt -> IO Gestalt) -> Int -> IO Gestalt
convergeCyclic g f burnIn = converge' burnedIn []
  where converge' mg table = do { g' <- mg
                                ; let mis = motifIndices g'
                                ; if mis `elem` trace (show $ length table) table
                                  then return g'
                                  else converge' (f g') (mis:table)
                                }
        burnedIn = iterateN burnIn (>>= f) (return g)
