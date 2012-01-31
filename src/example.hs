import Utils
import Wagner
main = do 
  seqs <- readSequences "../data/lexA_e_coli_120.csv"
  let mis = replicate 30 [52] ::[[Int]]
  let motifNum = 0 :: Int          
  let seqNum     = 0 :: Int          
  let i = seqNum
  let (mi,mis') = separate motifNum mis
  let (seq,seqs') = separate seqNum seqs
  let looPSSM = recoverNthPSSM seqs' mis' motifNum -- revise    
  let pssm = looPSSM
  let muMatrix = meanMatrix mis' :: MeanMatrix
  let varMatrix = varianceMatrix mis' :: VarMatrix
  let statMatrices = (muMatrix, varMatrix)
  let positions = [0..104] :: [Int]
  let g = Gestalt seqs mis
  return g