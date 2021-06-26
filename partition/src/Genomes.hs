{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

-- |
-- Module      : Genomes
-- Description : Representation of a genome
-- Copyright   : (c) Gabriel Siqueira, 2021
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
--
-- A genome comprises of a gene list and a integer list (correspondent to the sizes of intergenic regions).
module Genomes where

-- ( Duo,
--   Gene,
--   Genome (genomeIsSigned),
--   IR (..),
--   Idx (..),
--   Sign (..),
--   allSubGenomes,
--   balanced,
--   breakGenome,
--   combineGenomes,
--   combineGenomesL,
--   duoByIdx,
--   irByIdx,
--   duosList,
--   duoIdx,
--   duoIr,
--   duoToBS,
--   estimateITD,
--   fromLists,
--   genePair,
--   genomeSize,
--   intToGene,
--   interleaveListRepresentation,
--   interleaveListToGenome,
--   maybeEstimateITD,
--   occurenceMax,
--   randomGenome,
--   rearrangeGenome,
--   readGenome,
--   writeGenome,
--   sliceGenome,
--   reversal,
--   transposition,
--   subGenCount,
--   subGenome,
--   toLists,
--   validBegin,
--   validEnd,
--   -- weigth,
--   subGenomeFind,
--   GenomeMap,
--   gmEmpty,
--   gmEmptyRev,
--   gmLookup,
--   gmInsert,
--   gmLookupInsert,
-- )

import Control.Exception (assert)
import Control.Monad.Random (MonadRandom, getRandomRs, getRandoms)
import Data.ByteString.Builder (intDec, toLazyByteString)
import Data.ByteString.Lazy (ByteString)
import qualified Data.ByteString.Lazy.Char8 as BS
import Data.Coerce (coerce)
import Data.Foldable (foldl', toList)
import Data.Hashable (Hashable)
import Data.IntMap (IntMap)
import qualified Data.IntMap as IntMap
import qualified Data.List as List
import Data.Maybe (fromJust)
import Data.Text (Text)
import qualified Data.Text as Text
import Data.Vector (Vector, (!))
import qualified Data.Vector as Vec
import qualified Data.Vector.Mutable as MVec
import LocalBase
import System.Random (Random)
import System.Random.Shuffle (shuffleM)

newtype IR = IR Int deriving newtype (Eq, Show, Read, Num, Ord, Random)

newtype Idx = Idx Int deriving newtype (Eq, Show, Read, Num, Ord, Enum, Random, Integral, Real)

newtype Gene = Gene Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Random)

type Gstring = Vector Gene

type IRList = Vector IR

data Sign = Signed | Unsigned deriving (Eq, Show)

-- | Representation of a genome the gstring must be a non empty sequence of genes and
--  size of irList must be size of gstring minus one
data Genome = Genome
  { gstring :: Gstring,
    irList :: IRList,
    genomeIsSigned :: Sign
  }
  deriving (Eq)

data Duo = Duo
  { genePair :: (Gene, Gene),
    duoIr :: IR,
    duoIdx :: Idx,
    parGenomeSize :: Size,
    duoIsSigned :: Sign
  }
  deriving (Show)

----------------------------------------
--     Instantiations of Type Classes --
----------------------------------------

instance Orientable Gene where
  getOri a = if a >= 0 then LR else RL
  invOri a = - a

instance Orientable Duo where
  getOri duo =
    case duoIsSigned duo of
      Signed
        | getOri a1 == LR && getOri a2 == LR -> LR
        | getOri a1 == RL && getOri a2 == RL -> RL
        | otherwise -> if canonicOri a1 < canonicOri a2 then LR else RL
      Unsigned -> if canonicOri a1 < canonicOri a2 then LR else RL
    where
      (a1, a2) = genePair duo

  invOri duo@Duo {..} = case duoIsSigned of
    Signed -> duo {genePair = (invOri a2, invOri a1), duoIdx = idx'}
    Unsigned -> duo {genePair = (a2, a1), duoIdx = idx'}
    where
      (a1, a2) = genePair
      idx' = coerce parGenomeSize - duoIdx

instance Semigroup Sign where
  Signed <> Signed = Signed
  Signed <> Unsigned = Signed
  Unsigned <> Signed = Signed
  Unsigned <> Unsigned = Unsigned

instance Show Genome where
  show g =
    unwords . (("(" ++ head str_s ++ ")") :) $
      zipWith (\ir a -> "- " ++ ir ++ " - (" ++ a ++ ")") str_i (tail str_s)
    where
      str_s = Vec.toList $ show <$> gstring g
      str_i = Vec.toList $ show <$> irList g

instance Orientable Genome where
  getOri g =
    if
        | gstring g < rg -> LR
        | gstring g > rg -> RL
        | otherwise -> if irList g <= (Vec.reverse . irList $ g) then LR else RL
    where
      rg = Vec.reverse . gstring $ g
  invOri g = Genome (Vec.reverse . gstring $ g) (Vec.reverse . irList $ g) (genomeIsSigned g)

--------------------------------
--      Random Generation     --
--------------------------------

randomGenome :: MonadRandom mon => Size -> Int -> mon Genome
randomGenome size lim = do
  -- coins <- getRandoms
  -- ls <- zipWith swaps coins . take n . map fromIntegral <$> getRandomRs (1, lim)
  ls <- take n <$> getRandomRs (1, coerce lim)
  li <- take (n + 1) <$> getRandomRs (1, 100)
  return $ fromLists True Unsigned ls li
  where
    n = fromIntegral size

-- swaps b v = if b then v else invOri v

rearrangeGenome :: MonadRandom mon => Genome -> mon Genome
rearrangeGenome g = do
  let (sign, ls, li) = toLists True g
      s = sum li
  ls' <- shuffleM ls
  x <-
    List.sort . map IR . take (coerce $ genomeSize g) <$> getRandomRs (0, coerce s)
  let li' = zipWith (-) (x ++ [s]) (0 : x)
  return $ fromLists True sign ls' li'

--------------------------
--      Conversions     --
--------------------------

readGenome :: Bool -> Sign -> ByteString -> ByteString -> Genome
readGenome extend sign s i =
  fromLists
    extend
    sign
    (map (Gene . readInt) . BS.splitWith (\x -> x == ',' || x == ' ') $ s)
    (map (IR . readInt) . BS.splitWith (\x -> x == ',' || x == ' ') $ i)
  where
    readInt = coerce . fst . fromJust . BS.readInt

writeGenome :: Bool -> Genome -> (ByteString, ByteString)
writeGenome rext g =
  ( BS.unwords . fmap (toLazyByteString . intDec . coerce) $ ls,
    BS.unwords . fmap (toLazyByteString . intDec . coerce) $ li
  )
  where
    (_, ls, li) = toLists rext g

duoToBS :: Duo -> [ByteString]
duoToBS Duo {..} =
  let (al, ar) = genePair
   in [geneToBS al, irToBS duoIr, geneToBS ar]

intToGene :: Int -> Gene
intToGene = coerce

geneToBS :: Gene -> ByteString
geneToBS = toLazyByteString . (<>) "g" . intDec . coerce

irToBS :: IR -> ByteString
irToBS = toLazyByteString . (<>) "i" . intDec . coerce

fromLists :: Bool -> Sign -> [Gene] -> [IR] -> Genome
fromLists extend sign ls_ li = Genome (Vec.fromList ls) (Vec.fromList li) sign
  where
    ls =
      if extend
        then 0 : (ls_ ++ [if null ls_ then 1 else maximum ls_ + 1])
        else ls_

toLists :: Bool -> Genome -> (Sign, [Gene], [IR])
toLists rext g = (genomeIsSigned g, Vec.toList . (if rext then Vec.slice 1 (coerce $ genomeSize g - 2) else id) $ gstring g, Vec.toList $ irList g)

duosList :: Genome -> [Duo]
duosList g = foldr toDuo [] . zip [1 ..] . lPairs . Vec.toList $ gstring g
  where
    irs = irList g
    toDuo (i, (al, ar)) l = Duo (al, ar) (irs ! (i -1)) (Idx i) (genomeSize g) (genomeIsSigned g) : l

interleaveListRepresentation :: Genome -> ([ByteString], Sign)
interleaveListRepresentation g = (interleavelists ls li, genomeIsSigned g)
  where
    ls = toList . fmap geneToBS $ gstring g
    li = toList . fmap irToBS $ irList g

interleaveListToGenome :: [ByteString] -> Sign -> Genome
interleaveListToGenome l = Genome (Vec.fromList g_g) (Vec.fromList g_ir)
  where
    (g_g, g_ir) = foldr go ([], []) l
    go l (g, ir)
      | BS.head l == BS.head "g" = ((readG . BS.tail $ l) : g, ir)
      | BS.head l == BS.head "i" = (g, (readG . BS.tail $ l) : ir)
    readG :: (Num a) => ByteString -> a
    readG = fromIntegral . fst . fromJust . BS.readInt

------------------------------------------
--           Operations                 --
------------------------------------------

transposition :: Idx -> Idx -> Idx -> IR -> IR -> IR -> Genome -> Genome
transposition i j k x y z g =
  assert (2 <= i)
    . assert (i < j)
    . assert (j < k)
    . assert (k <= coerce (genomeSize g))
    . assert (0 <= x && x <= (vi ! (coerce i - 2)))
    . assert (0 <= y && y <= (vi ! (coerce j - 2)))
    . assert (0 <= z && z <= (vi ! (coerce k - 2)))
    $ Genome vs' vi' (genomeIsSigned g)
  where
    vs' = Vec.modify updateG vs
    vi' = Vec.modify updateIR vi

    updateG v =
      do
        aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
        aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
        MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
        MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1

    updateIR v = do
      do
        aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
        aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
        MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
        MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1
      MVec.write v (coerce i - 1 - 1) (x + y_rest)
      MVec.write v (coerce $ i + k - j - 2) (z + x_rest)
      MVec.write v (coerce k - 1 - 1) (y + z_rest)

    vs = gstring g
    vi = irList g
    x_rest = (vi ! (coerce i - 2)) - x
    y_rest = (vi ! (coerce j - 2)) - y
    z_rest = (vi ! (coerce k - 2)) - z

reversal :: Idx -> Idx -> IR -> IR -> Genome -> Genome
reversal i j x y g =
  assert (2 <= i)
    . assert (i < j)
    . assert (j <= coerce (genomeSize g) - 1)
    . assert (0 <= x && x <= (vi ! (coerce i - 2)))
    . assert (0 <= y && y <= (vi ! (coerce j - 1)))
    $ Genome vs' vi' (genomeIsSigned g)
  where
    vs' = Vec.modify updateG vs
    vi' = Vec.modify updateIR vi

    updateG v = do
      mapM_ (\k -> MVec.swap v (coerce $ i + k - 1) (coerce $ j - k - 1)) [0 .. (j - i + 1) `div` 2 - 1]
      case genomeIsSigned g of
        Unsigned -> pure ()
        Signed -> mapM_ (MVec.modify v invOri . coerce) [i .. j]

    updateIR v = do
      mapM_ (\k -> MVec.swap v (coerce i + coerce k - 1) (coerce j - coerce k - 2)) [0 .. (j - i + 1) `div` 2 - 1]
      MVec.write v (coerce i - 2) (x + y)
      MVec.write v (coerce j - 1) (x_rest + y_rest)

    vs = gstring g
    vi = irList g
    x_rest :: IR
    x_rest = (vi ! (coerce i - 2)) - x
    y_rest :: IR
    y_rest = (vi ! (coerce j - 1)) - y

irByIdx :: Genome -> Idx -> IR
irByIdx g idx = irList g ! (coerce idx - 1)

duoByIdx :: Genome -> Idx -> Duo
duoByIdx g idx = Duo (a1, a2) ir idx (genomeSize g) (genomeIsSigned g)
  where
    i = coerce idx
    vs = gstring g
    vi = irList g
    a1 = vs ! (i -1)
    a2 = vs ! i
    ir = vi ! (i -1)

combineGenomes :: IR -> Genome -> Genome -> Genome
combineGenomes ir g h = Genome (gstring g <> gstring h) (irList g <> Vec.cons ir (irList h)) (genomeIsSigned g <> genomeIsSigned h)

combineGenomesL :: [IR] -> [Genome] -> Genome
combineGenomesL irs gs =
  foldl' (\x (ir, y) -> combineGenomes ir x y) (head gs) (zip irs (tail gs))

breakGenome :: Genome -> Idx -> (IR, Genome, Genome)
breakGenome g idx = (irList g ! (coerce idx - 1), sliceGenome g 1 idx, sliceGenome g (idx + 1) (coerce $ genomeSize g))

sliceGenome :: Genome -> Idx -> Idx -> Genome
sliceGenome g i j = Genome vs' vi' (genomeIsSigned g)
  where
    vs' = Vec.slice (coerce $ i - 1) (coerce $ j - i + 1) $ gstring g
    vi' = Vec.slice (coerce $ i - 1) (coerce $ j - i) $ irList g

occurence :: Genome -> Gene -> Int
occurence g a = sum . fmap (\x -> if x == a then 1 else 0) . gstring $ g

occurenceMax :: Genome -> Int
occurenceMax g = maximum . fmap (g `occurence`) . gstring $ g

------------------------------------------
--           Predicates                 --
------------------------------------------

-- valid paths of a suffix tree must start with a gene
-- must end in a gene
validBegin :: ByteString -> Bool
validBegin bs = BS.head bs == BS.head "g"

-- valid ByteStrings correspondent to nodes of a suffix tree
-- must end in a gene
validEnd :: ByteString -> Bool
validEnd bs = BS.head bs == BS.head "g"

balanced :: Genome -> Genome -> Bool
balanced g h = balancedGenes && balancedIR
  where
    balancedGenes =
      (List.sort . map canonicOri . Vec.toList $ gstring g)
        == (List.sort . map canonicOri . Vec.toList $ gstring h)
    balancedIR = sum (irList g) == sum (irList h)

genomeSize :: Genome -> Size
genomeSize = Size . Vec.length . gstring

------------------------------------------
--           Sub-Genome                 --
------------------------------------------

subGenome :: Genome -> Genome -> Bool
subGenome x g = or $ subGenomeFind x g

-- | For each gene of a genome g, return if an occurence of
--  a nonempty genome x starts in that gene
subGenomeFind :: Genome -> Genome -> Vector Bool
subGenomeFind x g =
  fmap (\i -> (xs `isPrefixOfV` Vec.drop i gs) && (xi `isPrefixOfV` Vec.drop i gi)) inds
  where
    inds = Vec.elemIndices (Vec.head $ gstring x) gs
    gs = gstring g
    gi = irList g
    xs = gstring x
    xi = irList x
    isPrefixOfV v1 v2 = length v1 <= length v2 && v1 == Vec.unsafeTake (length v1) v2

allSubGenomes :: Genome -> [Genome]
allSubGenomes g = go 0 1 []
  where
    vs = gstring g
    vi = irList g
    n = coerce $ genomeSize g
    go i k acc
      | i > n = acc
      | i + k > n = go (i + 1) 1 acc
      | otherwise =
        let g' = Genome (Vec.slice i k vs) (Vec.slice i (k -1) vi) (genomeIsSigned g)
         in go i (k + 1) (g' : acc)

maybeEstimateITD :: Genome -> Genome -> Maybe Dist
maybeEstimateITD g h =
  if balanced g h
    then Just $ estimateITD g h
    else Nothing

-- | Count the number of occurrences of a nonempty genome x
--  in a genome g
subGenCount :: Genome -> Genome -> Int
subGenCount x g = sum . fmap (\b -> if b then 1 else 0) $ subGenomeFind x g

-- weigth :: Genome -> Genome -> Genome -> Int
-- weigth g h x = subGenCount x g - subGenCount x h

------------------------------------------
--              Distances               --
------------------------------------------

estimateITD :: Genome -> Genome -> Dist
estimateITD g h = undefined

----------------------------------
--        Genome Map            --
----------------------------------

data MapType = Direct | Reverse

data GenomeMap v = GM MapType (IntMap [(Genome, v)])

instance Show v => Show (GenomeMap v) where
  show (GM _ m) = show . concat . IntMap.elems $ m

gmEmpty :: GenomeMap v
gmEmpty = GM Direct IntMap.empty

gmEmptyRev :: GenomeMap v
gmEmptyRev = GM Reverse IntMap.empty

gmLookup :: (Orientable v) => Genome -> GenomeMap v -> Maybe v
gmLookup g m@(GM Direct _) = gmLookup_ g k m
  where
    k = coerce . Vec.head . gstring $ g
gmLookup g m@(GM Reverse _) =
  case getOri g of
    LR -> gmLookup_ g' k m
    RL -> invOri <$> gmLookup_ g' k m
  where
    g' = canonicOri g
    k = coerce . Vec.head . gstring $ g'

gmLookup_ :: Genome -> Int -> GenomeMap v -> Maybe v
gmLookup_ g k (GM _ m) = fmap snd $ IntMap.lookup k m >>= List.find ((==) g . fst)

gmInsert :: Genome -> v -> GenomeMap v -> GenomeMap v
gmInsert g v m@(GM mtype _) = gmInsert_ g' k v m
  where
    g' = case mtype of Direct -> g; Reverse -> canonicOri g
    k = coerce . Vec.head . gstring $ g'

gmInsert_ :: Genome -> Int -> v -> GenomeMap v -> GenomeMap v
gmInsert_ g k val (GM t m) = GM t m'
  where
    m' = case IntMap.lookup k m of
      Nothing -> IntMap.insert k [(g, val)] m
      Just _ -> IntMap.adjust ([(g, val)] ++) k m

-- if not found (return False and insert) else (return True)
gmLookupInsert :: (Orientable v) => Genome -> v -> GenomeMap v -> (GenomeMap v, v, Bool)
gmLookupInsert g val m@(GM Direct _) =
  case gmLookup_ g k m of
    Nothing -> (gmInsert_ g k val m, val, False)
    Just v -> (m, v, True)
  where
    k = coerce . Vec.head . gstring $ g
gmLookupInsert g val m@(GM Reverse _) =
  case gmLookup_ g' k m of
    Nothing -> (gmInsert_ g' k val m, val, False)
    Just v -> (m, case getOri g of LR -> v; RL -> invOri v, True)
  where
    g' = canonicOri g
    k = coerce . Vec.head . gstring . canonicOri $ g'
