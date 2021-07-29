{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}

-- |
-- Module      : Partition
-- Description : Representation and construction of genome partitions
-- Copyright   : (c) Gabriel Siqueira, 2021
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
--
-- A partition of two genomes G and H is essentially composed of two genome sequences S and P, such that:
--     - the sequence S when combined gives us G
--     - the sequence P when combined gives us H
--     - it is possible to rearrange S to obtain P
module Partition
  ( Partition,
    PartitionType (..),
    getPartition,
    validPartition,
    breakpoints,
    blocks,
    reduced,
    cost,
    -- Partition.weigth,
    sizeTmin,
  )
where

import Control.Arrow (second)
import Control.Exception (assert)
import Data.ByteString (ByteString)
import qualified Data.ByteString.Char8 as BS
import Data.Coerce (coerce)
import Data.Foldable (toList)
import Data.HashSet (HashSet)
import qualified Data.HashSet as HashSet
import qualified Data.List as List
import qualified Data.Map as Map
import Data.Maybe
  ( catMaybes,
    fromJust,
    fromMaybe,
    listToMaybe,
    mapMaybe,
  )
import Data.Sequence (Seq ((:<|), (:|>)))
import qualified Data.Sequence as Seq
import qualified Data.Set as Set
import Data.Vector (Vector, (!))
import qualified Data.Vector as Vec
import Debug.Trace
import GTree
import Genomes as G
import LocalBase

-- | Representation of a partition
data Partition = Partition {partType :: PartitionType, gseq :: Seq Genome, hseq :: Seq Genome, gbps :: Seq IR, hbps :: Seq IR} deriving (Show)

-- | A partial partition may be invalid
makePartialPartition :: PartitionType -> Genome -> Genome -> Partition
makePartialPartition ptype g h = Partition ptype (pure g) (pure h) mempty mempty

-- | Recives: a partition, a duo correspondent to a new breakpoint, and
--     an indication of the original genome and index of an occurrence
--     of a genome containing the breakpoint
--     Returns: a partition with the new breakpoint
addBreakpoint :: Partition -> Duo -> GenomePosition -> Partition
addBreakpoint part duo gp =
  case gp of
    (G idx _ _) ->
      let (gseq', gbps') = aux (gseq part) (gbps part) idx
       in part {gseq = gseq', gbps = gbps'}
    (H idx _ _) ->
      let (hseq', hbps') = aux (hseq part) (hbps part) idx
       in part {hseq = hseq', hbps = hbps'}
  where
    aux :: Seq Genome -> Seq IR -> Idx -> (Seq Genome, Seq IR)
    aux kseq kbps idx = (kseq', kbps')
      where
        (kseq_front, x, kseq_back, kbps_front, kbps_back, idx') = go Seq.Empty kseq Seq.Empty kbps idx

        (ir, xl, xr) = breakGenome x (duoIdx duo + idx' - 1)
        kseq' = kseq_front <> (xl :<| xr :<| kseq_back)
        kbps' = kbps_front <> (assert (ir == duoIr duo) ir :<| kbps_back)

        go ::
          Seq Genome ->
          Seq Genome ->
          Seq IR ->
          Seq IR ->
          Idx ->
          (Seq Genome, Genome, Seq Genome, Seq IR, Seq IR, Idx)
        go ys (x :<| xs) ybs xbss@(xb :<| xbs) idx' =
          let size = coerce (genomeSize x)
           in if idx' > size
                then go (ys :|> x) xs (ybs :|> xb) xbs (idx' - size)
                else (ys, x, xs, ybs, xbss, idx')
        go ys (x :<| xs) ybs xbss idx' = (ys, x, xs, ybs, xbss, idx')
        go _ Seq.Empty _ _ _ = error patternError

-- | 2k-approximation for the intergenic partition problem
getPartition :: PartitionType -> Genome -> Genome -> Partition
getPartition ptype g h = go (makePartialPartition ptype g h) tmin0 breaks0
  where
    tmin0 = makeTmin ptype g h
    breaks0 = Breaks ptype HashSet.empty
    go :: Partition -> Tmin -> Breaks -> Partition
    go part tmin breaks =
      case getGenome tmin of
        Nothing -> part
        Just (x, gp) -> case ptype of
          MCISP -> assert (tmin' /= tmin) $ go part' tmin' breaks'
          RMCISP -> assert (tmin'' /= tmin) $ go part' tmin'' breaks'
          where
            (break, breaks') = getBreak breaks x
            tmin' = updateTmin g h gp tmin break
            tmin'' = updateTmin (invOri g) (invOri h) (invOri gp) tmin' (invOri break)
            part' = addBreakpoint part break gp

-- | A partition (s,p) is valid if
--  it is possible to rearrange s to obtain p
validPartition :: Partition -> Bool
validPartition = uncurry balanced . reduced

breakpoints :: Partition -> (Seq IR, Seq IR)
breakpoints part = (gbps part, hbps part)

blocks :: Partition -> (Seq Genome, Seq Genome)
blocks part = (gseq part, hseq part)

reduced :: Partition -> (Genome, Genome)
reduced part = (gr, hr)
  where
    gr = fromLists False sign gr_ls (toList $ gbps part)
    hr = fromLists False sign hr_ls (toList $ hbps part)
    (gr_ls, hr_ls) = (genomesToGenes ggs, genomesToGenes hhs)
    ggs = gseq part
    hhs = hseq part
    (sign, gmEmptyX) = case partType part of
      MCISP -> (Unsigned, gmEmpty)
      RMCISP -> (Signed, gmEmptyRev)

    genomesToGenes = map genomeToGene . toList
    genomeToGene g = fromMaybe (error "Error on reduced.") $ gmLookup g m
    m = fst $ foldl addGenome (gmEmptyX, 1 :: Gene) hhs
    addGenome (m, count) g = (m', count')
      where
        count' = if keepOld then count else count + 1
        (m', _, keepOld) = gmLookupInsert g count m

-- | The cost of a partition (s,p) is the number of
--  breakpoints from s
cost :: Partition -> Int
cost = length . fst . breakpoints

-- weigth :: Partition -> Genome -> Int
-- weigth part x = sumCounts bls1 - sumCounts bls2
--   where
--     sumCounts = sum . fmap (subGenCount' x)
--     subGenCount' x g =
--         case partType part of
--           MCISP -> subGenCount x g
--           RMCISP -> subGenCount x g + subGenCount (invOri x) g
--     (bls1, bls2) = blocks part

listGenomes :: PartitionType -> Genome -> Genome -> [Genome]
listGenomes ptype g h = go [] tmin0 breaks0
  where
    tmin0 = makeTmin ptype g h
    breaks0 = Breaks ptype HashSet.empty
    go :: [Genome] -> Tmin -> Breaks -> [Genome]
    go acc tmin breaks =
      case getGenome tmin of
        Nothing -> acc
        Just (x, gp) -> go acc' tmin'' breaks'
          where
            (break, breaks') = getBreak breaks x
            tmin' = updateTmin g h gp tmin break
            tmin'' = updateTmin (invOri g) (invOri h) (invOri gp) tmin' (invOri break)
            acc' = x : acc

sizeTmin :: PartitionType -> Genome -> Genome -> Size
sizeTmin ptype g h = Size . length $ listGenomes ptype g h

-- | Breaks represents a set of breakpoints from a partial partition
data Breaks = Breaks PartitionType (HashSet (Gene, Gene)) deriving (Show)

isBreak :: Duo -> Breaks -> Bool
isBreak duo (Breaks partType brks) = (genePair . canonicOri $ duo) `HashSet.member` brks

getBreak :: Breaks -> Genome -> (Duo, Breaks)
getBreak brks@(Breaks partType set) g =
  case dropWhile (not . (`isBreak` brks)) ls of
    [] -> (head ls, Breaks partType (HashSet.insert (genePair . canonicOri $ head ls) set))
    (l : _) -> (l, brks)
  where
    ls = duosList g
