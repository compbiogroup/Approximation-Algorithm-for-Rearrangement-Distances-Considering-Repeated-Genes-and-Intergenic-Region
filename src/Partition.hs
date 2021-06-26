{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
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
  -- ( Partition,
  --   PartitionType (..),
  --   getPartition,
  --   validPartition,
  --   breakpoints,
  --   blocks,
  --   reduced,
  --   cost,
  --   -- Partition.weigth,
  --   makeTmin,
  --   sizeTmin,
  -- )
where

import Control.Arrow (second)
import Control.Exception (assert)
import Data.ByteString.Lazy (ByteString)
import Data.Coerce (coerce)
import Data.Foldable (toList)
import Data.HashSet (HashSet)
import qualified Data.HashSet as HashSet
import Data.IntMap (IntMap)
import qualified Data.IntMap as IntMap
import qualified Data.List as List
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
import qualified Data.SuffixTree as ST
import Data.Tree (Forest, Tree)
import qualified Data.Tree as Tree
import Debug.Trace
import Genomes as G
import LocalBase

-- | Indication of diferent partition types
data PartitionType = MCISP | RMCISP deriving (Show, Read, Enum, Bounded, Eq)

-- | Representation of a partition
data Partition = Partition {partType :: PartitionType, gseq :: Seq Genome, hseq :: Seq Genome, gbps :: Seq IR, hbps :: Seq IR} deriving (Show)

-- | Position of a subgenome X in the indicated genome Y (G or H). The first value (Idx) is the index of X in Y. The second value (Size) is the size of Y minus the size of X. The third value (Ori) indicate whether X and Y are inverted.
data GenomePosition = G Idx Size Ori | H Idx Size Ori deriving (Show, Eq)

instance Orientable GenomePosition where
  getOri (G _ _ LR) = LR
  getOri (G _ _ RL) = RL
  getOri (H _ _ LR) = LR
  getOri (H _ _ RL) = RL

  invOri (G idx remSize LR) = G (coerce remSize - idx + 2) remSize RL
  invOri (G idx remSize RL) = G (coerce remSize - idx + 2) remSize LR
  invOri (H idx remSize LR) = H (coerce remSize - idx + 2) remSize RL
  invOri (H idx remSize RL) = H (coerce remSize - idx + 2) remSize LR

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
        Just (x, gp) -> assert (tmin'' /= tmin) go part' tmin'' breaks'
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

-----------------------------------------------------
-- Auxiliar structures to produce the partition    --
-----------------------------------------------------

-- | Tmin is the minimal set of set T, where
--  T is the set of genomes with different
--  weights in g or h
data Tmin = Tmin PartitionType Sign (HSSTree ByteString) deriving (Show, Eq)

-- | SuffixTree from the Hitting Set Algorithm, we use it to construct the set Tmin.
--  Leafs correspondent to suffixes started in G are marked with G and Leafs correspondent
--  to suffixes started in H are marked with H.
--  Each node has SumG with the number of leafs of the subtree marked as G and
--  SumH with the number of leafs of the subtree marked as H.
--  Size with the genome size
--  does not contain the characters marking the end of strings.
data HSSTree a = HSRoot
  { rSize :: Size,
    rChildren :: [Tree (HSLabel a)]
  }
  deriving (Eq)

instance (Show a) => Show (HSSTree a) where
  show HSRoot {..} =
    "Genome Size:"
      ++ show rSize
      ++ "\n"
      ++ Tree.drawForest (map (fmap show) rChildren)

type HSSubTree = Tree (HSLabel ByteString)

type SumG = Int

type SumH = Int

data HSLabel a = HSLabel
  { hsInfo :: Info,
    hsPref :: [a],
    hsSize :: Size
  }
  deriving (Eq)

data IsInG = InG | InH deriving (Eq, Show)

-- A leaf is available if it is proper (the correspondent infix is not empty)
-- and its first intergenic region is not a breakpoint or after a breakpoint
-- the lastBreakDistance is the number of characters after the last breakpoint
-- of the correspondent suffix
data Info
  = LeafInfo
      { lastBreakDistance :: Int,
        isInG :: IsInG,
        isRev :: Ori,
        available :: Bool
      }
  | NodeInfo
      { sumG :: SumG,
        sumH :: SumH
      }
  deriving (Eq)

instance (Show a) => Show (HSLabel a) where
  show HSLabel {..} =
    "Prefix:" ++ show hsPref
      ++ ", Size:"
      ++ show hsSize
      ++ ( case hsInfo of
             LeafInfo {..} ->
               (if isInG == InG then " - G" else " - H")
                 ++ ", up:"
                 ++ show lastBreakDistance
                 ++ (if isRev == LR then "" else " R")
                 ++ (if available then "" else " X")
             NodeInfo {..} -> " - sumG: " ++ show sumG ++ ", sumH:" ++ show sumH
         )

makeTmin :: PartitionType -> Genome -> Genome -> Tmin
makeTmin ptype g h = Tmin ptype (genomeIsSigned g) $ makeHSSTree (ST.construct str)
  where
    n = genomeSize g
    (l1, _) = interleaveListRepresentation g
    (l2, _) = interleaveListRepresentation h
    (rl1, _) = interleaveListRepresentation (invOri g)
    (rl2, _) = interleaveListRepresentation (invOri h)
    str = case ptype of
      MCISP -> l1 ++ ("$" : l2) ++ ["#"]
      RMCISP -> l1 ++ ("$" : l2) ++ ("#" : rl1) ++ ("%" : rl2) ++ ["&"]
    -- Markers for end of strings
    markers = case ptype of
      MCISP -> ["$", "#"]
      RMCISP -> ["$", "#", "%", "&"]
    getG t = case hsInfo . Tree.rootLabel $ t of
      LeafInfo {..} -> if isInG == InG then 1 else 0
      NodeInfo {..} -> sumG
    getH t = case hsInfo . Tree.rootLabel $ t of
      LeafInfo {..} -> if isInG == InH then 1 else 0
      NodeInfo {..} -> sumH

    -- Construct the HSSTree:
    makeHSSTree ST.Leaf = error patternError
    makeHSSTree (ST.Node edges) =
      HSRoot n ((\ts -> assert (fmap getG ts == fmap getH ts) ts) subTrees)
      where
        subTrees = concat $ mapMaybe makeSubHSSTree edges
        makeSubHSSTree e@(p_, t) =
          let p = ST.prefix p_
           in if
                  | not $ validBegin (head p) -> Nothing
                  | length p == 1 -> Just $ go Nothing e
                  | otherwise ->
                    Just $
                      go Nothing (ST.mkPrefix [head p], ST.Node [(ST.mkPrefix (tail p), t)])

    -- inG indicate whether the nodes suffix starts in G
    -- hasMarker indicate whether a marker was found
    go maybeHead (p_, st) =
      case st of
        ST.Leaf ->
          pure $
            flip Tree.Node [] $
              if
                  | "$" `elem` p__ ->
                    let p = takeWhile (/= "$") p__
                     in HSLabel
                          (LeafInfo 0 InG LR (head p__ /= "$"))
                          p
                          (Size . length $ p)
                  | "#" `elem` p__ ->
                    let p = takeWhile (/= "#") p__
                     in HSLabel
                          (LeafInfo 0 InH LR (head p__ /= "#"))
                          p
                          (Size . length $ p)
                  | "%" `elem` p__ ->
                    let p = takeWhile (/= "%") p__
                     in HSLabel
                          (LeafInfo 0 InG RL (head p__ /= "%"))
                          p
                          (Size . length $ p)
                  | "&" `elem` p__ ->
                    let p = takeWhile (/= "&") p__
                     in HSLabel
                          (LeafInfo 0 InH RL (head p__ /= "&"))
                          p
                          (Size . length $ p)
                  | otherwise -> error patternError
        (ST.Node edges) ->
          if size == 0
            then subTrees
            else [Tree.Node (HSLabel (NodeInfo sG sH) p (Size . length $ p)) subTrees]
          where
            size = Size . length $ p
            sG = sum . map getG $ subTrees
            sH = sum . map getH $ subTrees
            subTrees = concatMap makeSubHSSTree edges
            makeSubHSSTree = go maybeHead'
            (maybeHead', p) =
              (\h -> if validEnd h then (Nothing, p__) else (Just h, init p__)) $ last p__
      where
        p__ = case maybeHead of
          Just h -> h : ST.prefix p_
          Nothing -> ST.prefix p_

class WalkDownInfo w where
  goDown :: w -> HSSubTree -> [(w, HSSubTree)]

-- information stored while walking on the tree
data SearchDownInfo = SearchDown
  { subStrSize :: Size, -- node's substring size
    subStr :: [ByteString], -- node's substring
    accPrefSize :: Size, -- accummulated prefix size
    accPref :: [ByteString] -- accummulated prefix  (reversed)
  }

makeSearchDownInfo :: HSSubTree -> SearchDownInfo
makeSearchDownInfo t@(Tree.Node HSLabel {..} children) =
  if hsSize >= 2
    then SearchDown 2 [head hsPref, head . tail $ hsPref] hsSize (reverse hsPref)
    else SearchDown 1 [head hsPref] hsSize (reverse hsPref)

instance WalkDownInfo SearchDownInfo where
  goDown SearchDown {..} currentNode =
    let (Tree.Node HSLabel {..} children) = currentNode
     in case hsInfo of
          LeafInfo {..} -> []
          NodeInfo {..} -> fmap aux children
    where
      aux t@(Tree.Node HSLabel {..} children) =
        (,t) $
          SearchDown
            (min 2 hsSize + accPrefSize)
            (reverse $ take 2 hsPref ++ accPref)
            (Size (length hsPref) + accPrefSize)
            (reverse hsPref ++ accPref)

data UpdateInfo = Update
  { currentStr :: [ByteString], -- remaining of selected suffix
    breakDistance :: Idx, -- number of characters until the breakpoint is reached
    updateDistance :: Idx, -- number of characters until update starts
    leafInG :: IsInG -- whether the suffix's leaf is in G
  }
  deriving (Show)

makeUpdateInfo :: [ByteString] -> Idx -> UpdateInfo
makeUpdateInfo str bDist = Update str bDist 0 InG

instance WalkDownInfo (Maybe UpdateInfo) where
  goDown Nothing _ = []
  goDown (Just up@Update {..}) currentNode =
    case hsInfo of
      LeafInfo {..} -> []
      NodeInfo {..} -> map aux children
    where
      (Tree.Node HSLabel {..} children) = currentNode
      aux :: HSSubTree -> (Maybe UpdateInfo, HSSubTree)
      aux t@(Tree.Node HSLabel {..} children) =
        if hsPref `List.isPrefixOf` currentStr
          then (Just $ up {currentStr = drop (coerce hsSize) currentStr}, t)
          else (Nothing, t)

breakNode :: HSSubTree -> Idx -> IsInG -> HSSubTree
breakNode (Tree.Node HSLabel {..} children) idx inG =
  case hsInfo of
    LeafInfo {..} ->
      Tree.Node
        ( HSLabel
            ( if
                  | not available -> NodeInfo 0 0
                  | isInG == InG -> NodeInfo 1 0
                  | otherwise -> NodeInfo 0 1
            )
            pl
            (coerce idx - 1)
        )
        [Tree.Node (HSLabel hsInfo {available = False} pr (hsSize - coerce idx + 1)) children]
    NodeInfo {..} ->
      Tree.Node
        (HSLabel (NodeInfo sumG sumH) pl (coerce idx - 1))
        [ Tree.Node
            ( if inG == InG
                then HSLabel (NodeInfo (sumG - 1) sumH) pr (hsSize - coerce idx + 1)
                else HSLabel (NodeInfo sumG (sumH - 1)) pr (hsSize - coerce idx + 1)
            )
            children
        ]
  where
    (pl, pr) = splitAt (coerce idx - 1) hsPref

updateNode :: UpdateInfo -> HSSubTree -> (UpdateInfo, HSSubTree)
updateNode up@Update {..} t@(Tree.Node hsl@HSLabel {..} children) =
  if breakDistance <= 0
    then (up, t)
    else
      (info',) $
        if updateDistance <= 0 && breakDistance' <= -2
          then breakNode node' (- breakDistance') leafInG
          else node'
  where
    breakDistance' = breakDistance - coerce hsSize
    (info', node') =
      case hsInfo of
        LeafInfo {..} ->
          let lastBreakDistance' = max (coerce breakDistance) lastBreakDistance
           in ( up
                  { updateDistance = coerce lastBreakDistance - coerce hsSize,
                    breakDistance = breakDistance',
                    leafInG = isInG
                  },
                Tree.Node
                  ( hsl
                      { hsInfo =
                          hsInfo
                            { lastBreakDistance = lastBreakDistance',
                              available = available && breakDistance' <= -2
                            }
                      }
                  )
                  children
              )
        NodeInfo {..} ->
          if updateDistance <= 0 && breakDistance' > -2
            then
              if leafInG == InG
                then
                  ( up {breakDistance = breakDistance'},
                    Tree.Node (hsl {hsInfo = hsInfo {sumG = sumG - 1}}) children
                  )
                else
                  ( up {breakDistance = breakDistance'},
                    Tree.Node (hsl {hsInfo = hsInfo {sumH = sumH - 1}}) children
                  )
            else
              ( up {updateDistance = updateDistance - coerce hsSize, breakDistance = breakDistance'},
                t
              )

-- | Recover a sub-genome of some block also returns
--  a genome position of one of its occurrences
getGenome :: Tmin -> Maybe (Genome, GenomePosition)
getGenome (Tmin ptype sign HSRoot {..}) =
  fmap toGenome . safeMinimum
    . mapMaybe (\t -> walkDown (makeSearchDownInfo t) t)
    $ rChildren
  where
    toGenome :: ([ByteString], Size, Size, IsInG, Ori) -> (Genome, GenomePosition)
    -- Size is the size of the suffix after the sub-genome occurrence
    -- and the boolean indicate whether the occurrence is in A
    toGenome (bs, _, suf_size, inG, ori) =
      (g,) $
        if inG == InG
          then G idx (rSize - genomeSize g) LR
          else H idx (rSize - genomeSize g) LR
      where
        g_ = interleaveListToGenome bs sign
        (idx, g) = case ori of
          LR -> (coerce $ rSize - suf_size `div` 2 - genomeSize g + 1, g_)
          RL -> (coerce $ suf_size `div` 2 + 1, invOri g_)

    -- Find minimum element of T'.
    -- Note: a first node is never chosen,
    -- because it can not have a breakpoint (so is never unbalanced)
    walkDown :: SearchDownInfo -> HSSubTree -> Maybe ([ByteString], Size, Size, IsInG, Ori)
    walkDown gd@SearchDown {..} currentNode =
      let (Tree.Node HSLabel {..} children) = currentNode
       in case hsInfo of
            LeafInfo {..} ->
              if available
                then Just (subStr, subStrSize, hsSize - 2, isInG, isRev)
                else Nothing
            NodeInfo {..} ->
              if sumG /= sumH
                then
                  let inG = if sumG > sumH then InG else InH
                   in (\(ori, sp) -> (subStr, subStrSize, sp, inG, ori))
                        . second (subtract 2)
                        <$> sum_pfx inG 0 currentNode
                else safeMinimum . mapMaybe (uncurry walkDown) $ goDown gd currentNode
    sum_pfx :: IsInG -> Size -> HSSubTree -> Maybe (Ori, Size)
    sum_pfx inG acc (Tree.Node HSLabel {..} children) =
      case hsInfo of
        LeafInfo {..} ->
          let acc' = acc + hsSize
           in if inG == isInG && lastBreakDistance < coerce acc' - 1
                then Just (isRev, acc')
                else Nothing
        NodeInfo {..} -> listToMaybe $ mapMaybe (sum_pfx inG (hsSize + acc)) children
    safeMinimum [] = Nothing
    safeMinimum l = Just $ List.minimumBy (\(_, a, _, _, _) (_, b, _, _, _) -> a `compare` b) l

listGenomes :: Genome -> Genome -> Tmin -> [Genome]
listGenomes g h tmin0@(Tmin ptype _ _) = go [] tmin0 breaks0
  where
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

updateTmin :: Genome -> Genome -> GenomePosition -> Tmin -> Duo -> Tmin
updateTmin g h gp tmin duo = foldl aux tmin suffixesAndIndexes
  where
    (suffixesAndIndexes, inG, ori) = case gp of
      (G gidx _ ori) -> (select gidx g, InG, ori)
      (H gidx _ ori) -> (select gidx h, InH, ori)
    n = genomeSize g
    select gidx =
      let prefSize = gidx + duoIdx duo - 1
          i = 2 * (coerce n - prefSize) - 1
       in zip (repeat i) . take (coerce prefSize)
            . evens
            . List.tails
            . fst
            . interleaveListRepresentation

    -- search the suffix in each subtree of root
    -- Note: a first node never have a breakpoint
    aux :: Tmin -> (Idx, [ByteString]) -> Tmin
    aux (Tmin ptype sign root@HSRoot {..}) (i, s) =
      Tmin ptype sign $ root {rChildren = map aux2 rChildren}
      where
        aux2 t@(Tree.Node l@HSLabel {..} children) =
          if hsPref `List.isPrefixOf` s
            then snd . go (makeUpdateInfo (drop (coerce hsSize) s) i) $ t
            else t

    -- return the updated node
    go :: UpdateInfo -> HSSubTree -> (Maybe UpdateInfo, HSSubTree)
    go info0@Update {..} t@(Tree.Node l@HSLabel {..} children) =
      case goDown (Just info0) t of
        [] -> case hsInfo of
          LeafInfo {..} ->
            let (info', t'@(Tree.Node l' _)) = updateNode info0 t
             in if inG == isInG && null currentStr && ori == isRev
                  then (Just info', t')
                  else (Nothing, t)
          NodeInfo {..} -> (Nothing, t)
        infosTrees ->
          let (infos_, children') =
                unzip $
                  map
                    ( \(mInfo, t) ->
                        case mInfo of
                          Nothing -> (Nothing, t)
                          Just info -> go info t
                    )
                    infosTrees
              infos = catMaybes infos_
              info' = head infos
              (info'', t') = updateNode info' (Tree.Node l children')
           in if null infos
                then (Nothing, t)
                else assert (length infos == 1) (Just info'', t')

sizeTmin :: Genome -> Genome -> Tmin -> Size
sizeTmin g h = Size . length . listGenomes g h

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
