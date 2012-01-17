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
  where freqs =  [(fromIntegral count)/ (fromIntegral len) | count <- counts]
        counts = getCounts col
        len = length col

entropyEpsilon = 10*10e-10
entropy xs = (- 1) * (sum $ (map (\x -> x * log2 (x + entropyEpsilon)) xs))

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

bindingEnergyAt :: PSSM -> Sequence -> Index -> Float
bindingEnergyAt pssm seq i = exp (- score)
  where score = scoreAt pssm seq i

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

potential :: Sequence -> NamedPSSM -> Index -> MotifIndex -> MotifIndices -> VarMatrix -> Float
--potential can't be larger than 700, or exp (-potential) will underflow
potential seq (i,pssm) pos mi mis varMatrix = (bindingEnergy + a * stringEnergy)/700
  where bindingEnergy = printBE $ bindingEnergyAt pssm seq i --bigger is worse
        stringEnergy = printSE $ sum [log $ epsilon + energyFromString j jpos
                           | (j, jpos) <- zip [0..] mi, j /= i]
        energyFromString j jpos =printEFS $ 1/ (epsilon + var j) * fromIntegral ((pos - jpos) - 1) ** 2
        var j = varMatrix !! i !! j
        a = 0
        muMatrix = meanMatrix mis
        mu i j = muMatrix !! i !! j
        
meanMatrix :: MotifIndices -> [[Float]] --compute resting lengths
--matrix is symmetric, upper triangular; could just compute half of it
meanMatrix mis = [[mean [((mi!!i) - (mi!!j)) | mi <- mis]
                  |i <- motifRange] 
                 | j <- motifRange]
  where motifRange = [0..numMotifs - 1]
        
        
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
  let varMatrix = varianceMatrix (delete mi mis)
  i' <- assignIthIndex (seqNum,seq) (motifNum, looPSSM) mis varMatrix
  let  mi' = replaceAt motifNum i' mi
  let mis' = replaceAt seqNum mi' mis
  return (Gestalt seqs mis')

sa :: Gestalt -> IO Gestalt
sa (Gestalt seqs mis) = do
  seqNum <- randomRIO (0, length seqs - 1)
  let seq = seqs !! seqNum      
  let mi = mis !! seqNum      
  motifNum <- randomRIO (0, numMotifs - 1)
  let varMatrix = varianceMatrix (delete mi mis)
  let looPSSM = recoverNthPSSM (delete seq seqs) (delete mi mis) motifNum -- revise
  proPos <- assignIthIndex (seqNum,seq) (motifNum, looPSSM) mis varMatrix
  let curPos = mi !! motifNum
  let curPot = potential seq (motifNum,looPSSM) curPos mi mis varMatrix
  let proPot = potential seq (motifNum,looPSSM) proPos mi mis varMatrix
  r <- randomRIO (0.0,1.0)
  let accept = (proPot < curPot) || (r < curPot / proPot)
  let nextPos = if accept then proPos else curPos
  let  mi' = replaceAt motifNum nextPos mi
  let mis' = replaceAt seqNum mi' mis
  return (Gestalt seqs mis')


assignIthIndex :: NamedSequence -> NamedPSSM -> MotifIndices -> VarMatrix -> IO Index

assignIthIndex (seqNum,seq) (i,pssm) mis varMatrix =
  sample positions (\pos ->printLikelihood $ exp (- energy pos))--via Boltzmann distribution
  where end = length seq - length pssm --check this
        positions = [0..end]
        (mi, mis') = separate seqNum mis
        energy pos = printPotential $ potential seq (i,pssm) pos mi mis' varMatrix
    
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
  where mic = (transpose mis) !! n

recoverPSSMs :: Gestalt -> [PSSM]
recoverPSSMs gestalt = map (`recoverPSSM` seqs) mics
  where mics = transpose $ motifIndices gestalt
        seqs = sequences gestalt

recoverMotif :: MotifIndexCol -> Sequences -> Motif
recoverMotif = zipWith (\m s -> (take motifLength . drop m) s)

recoverMotifs :: Gestalt -> [Motif]
recoverMotifs g = map (\mic -> recoverMotif mic seqs) mics
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

converge :: Gestalt -> IO Gestalt
converge g = converge' g (updateAlignment g)
  where converge' g mg = do { g' <- mg
                            ; if motifIndices g == motifIndices g'
                                 then return g
                              else converge' g' (updateAlignment g')
                            }
-- debugging

printPotential x | trace ("printPotential"++ " " ++ show x) False = undefined
printPotential x = x
printTubs xs | trace ("printTubs"++ " " ++ show xs) False = undefined
printTubs xs = xs
printSE x | trace ("printSE"++ " " ++ show x) False = undefined
printSE x = x
printFaks xs | trace ("maxFax"++ " " ++ show (maximum xs / sum xs)) False = undefined
printFaks xs = xs
printZ xs | trace ("printZ"++ " " ++ show xs) False = undefined
printZ xs = xs
printLikelihood x | trace ("printLikelihood"++ " " ++ show x) False = undefined
printLikelihood x = x
printEFS x | trace ("printEFS"++ " " ++ show x) False = undefined
printEFS x = x
printBE x | trace ("printBE"++ " " ++ show x) False = undefined
printBE x = x








