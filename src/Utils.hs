module Utils where 

import Data.List
import Data.Char
import System.Random hiding (split)
import Debug.Trace

range :: (Integral a) => a -> [a]
range n = [0..(n-1)]

log2 :: (Floating a) => a -> a
log2 = logBase 2

getCounts :: Ord a => [a] -> [Int]
getCounts xs = map length ((group . sort) xs)

argMax :: (Ord b) => (a -> b) -> [a] -> a
argMax f = foldl1 (\x x' -> if f x' > f x then x' else x)

argMin :: (Ord b) => (a -> b) -> [a] -> a
argMin f = foldl1 (\x x' -> if f x' <= f x then x' else x)
  
trim :: String -> String --stole this from wikipedia for portability
trim = f . f
  where f = reverse . dropWhile isSpace

chooseRandomly :: [a] -> IO a
chooseRandomly xs = do
  r <- randomRIO (0, length xs - 1)
  return (xs !! r)

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
-- greediness constant k)
--sample' as f r | trace ("printsample'"++ " " ++ show as++ " ") False = undefined
sample' as f r = fst $ argMin snd $ filter ((>= r) . snd)  tups
              where k = 1
                    faks =  printFaks $ map (\a -> f a ** k) as
                    z = sum faks
                    normFaks = printNormFaks $ map (/z) faks
                    tups = printTups $ zip as (scanl1 (+) normFaks)

--printFaks x | trace ("printFaks"++ " " ++ show x) False = undefined
printFaks x = x
--printNormFaks x | trace ("printNormFaks"++ " " ++ show x) False = undefined
printNormFaks x = x
--printTups x | trace ("printTups"++ " " ++ show x) False = undefined
printTups x = x 


separate ::  Int -> [a] -> (a,[a])
separate i seqs = (seqs !! i, removeNth seqs i)

insertAt ::  Int -> a -> [a] -> [a] 
insertAt i a as = take i as ++ [a] ++ drop i as 

replaceAt ::  Int -> a -> [a] -> [a] 
replaceAt i a as = take i as ++ [a] ++ drop (i + 1) as 

readSequences :: FilePath -> IO [String]
readSequences filePath = do
  content <- readFile filePath
  return (sanitizeFASTA content)
  
sanitizeFASTA :: String -> [String]
sanitizeFASTA content = map (filter (/= ',')) relevantLines
  where ls = map trim (lines content)
        relevantLines = filter ((/= '>') . head) ls
        
removeNth ::  [a] -> Int -> [a]
removeNth xs n = ys ++ tail zs
  where (ys,zs) = splitAt n xs   

split :: Char -> String -> [String]
split ch xs = reverse $ map reverse $ split' ch xs [[]]

split' :: Char -> String -> [String] -> [String]
split' ch xs acc
  | null xs = acc
  | ch == head xs = split' ch (tail xs) ([]:acc)
  | otherwise = split' ch (tail xs) (((head xs):(head acc)):(tail acc))
  
replace :: (Eq a) => [a] -> [a] -> [a] -> [a]
replace old new str
  | null str = []
  | old == take l str = new ++ replace old new (drop l str)
  | otherwise = head str : replace old new (tail str)
    where l = length old
  
iterateN' ::  Int -> (a -> a) -> a -> a
iterateN' 0 f x = x
iterateN' n f x = x `seq` iterateN' (n - 1) f (f x)

iterateNverbose ::  Int -> (a -> a) -> a -> (a,Int)
iterateNverbose n f x = (iterateN' n f x, n)

iterateN :: Int -> (a -> a) -> a -> a 
iterateN 0 f x = x
iterateN n f x = iterateN (n - 1) f (f x)

matrixMap :: (a -> b) -> [[a]] -> [[b]]
matrixMap f = map (map f) 

selectColumn ::  [[a]] -> Int -> [a]
selectColumn xss i = [xs !! i' | xs <- xss]
  where i' = fromIntegral i

variance :: (Real b, Floating b, Floating a) => [b] -> a
variance xs = mean (map (**2) xs) - mean xs ** 2

mean :: (Real a, Fractional b) => [a] -> b
mean xs = realToFrac (sum xs) / genericLength xs --Thanks, Don Stewart

fixpoint :: (Eq b) => (a -> a) -> a -> (a -> b) -> a
fixpoint f a p = fst $ head $ dropWhile (\(x,y) -> p x /= p y) $ zip its (tail its)
  where its = iterate f a

-- Extend composition operator to functions of more arguments
compose4 :: (b -> c) -> (a -> a1 -> a2 -> a3 -> b) -> a -> a1 -> a2 -> a3 -> c
compose4 = (.).(.).(.).(.)

normalize xs = map (/sum xs) xs

mmap :: (a -> b) -> [[a]] -> [[b]]
mmap f = map (map f)