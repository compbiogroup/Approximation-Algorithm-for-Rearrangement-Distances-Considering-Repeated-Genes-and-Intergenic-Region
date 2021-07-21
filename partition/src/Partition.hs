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
    makeTmin,
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
import Data.Tree (Forest, Tree)
import qualified Data.Tree as Tree
import Data.Vector (Vector, (!))
import qualified Data.Vector as Vec
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
        Just (x, gp) -> case ptype of
          MCISP -> assert (tmin' /= tmin) go part' tmin' breaks'
          RMCISP -> assert (tmin'' /= tmin) go part' tmin'' breaks'
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
data Tmin = Tmin PartitionType Sign HSSTree deriving (Show, Eq)

-- | SuffixTree from the Hitting Set Algorithm, we use it to construct the set Tmin.
--  Leafs correspondent to suffixes started in G are marked with G and Leafs correspondent
--  to suffixes started in H are marked with H.
--  Each node has SumG with the number of leafs of the subtree marked as G and
--  SumH with the number of leafs of the subtree marked as H.
--  Size with the genome size
--  does not contain the characters marking the end of strings.
data HSSTree = HSRoot
  { rSize :: Size,
    rChildren :: [Tree HSLabel]
  }
  deriving (Eq)

instance Show HSSTree where
  show HSRoot {..} =
    "Genome Size:"
      ++ show rSize
      ++ "\n"
      ++ Tree.drawForest (map (fmap show) rChildren)

type HSSubTree = Tree HSLabel

type SumG = Int

type SumH = Int

data HSLabel = HSLabel
  { hsInfo :: Info,
    hsPref :: IdxPair
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

instance Show HSLabel where
  show HSLabel {..} =
    "Prefix:" ++ show hsPref
      ++ ( case hsInfo of
             LeafInfo {..} ->
               (if isInG == InG then " - G" else " - H")
                 ++ ", up:"
                 ++ show lastBreakDistance
                 ++ (if isRev == LR then "" else " R")
                 ++ (if available then "" else " X")
             NodeInfo {..} -> " - sumG: " ++ show sumG ++ ", sumH:" ++ show sumH
         )

data IdxPair = IdxPair (Vector ByteString) Int Int

instance Show IdxPair where
  show idxPair = show $ ipSlice idxPair

instance Eq IdxPair where
  (IdxPair _ lidx1 ridx1) == (IdxPair _ lidx2 ridx2) = lidx1 == lidx1 && ridx2 == ridx2

addHead :: IdxPair -> IdxPair
addHead (IdxPair v lidx ridx) = IdxPair v (lidx - 1) ridx

getHead :: IdxPair -> ByteString
getHead (IdxPair v lidx _) = v ! lidx

getLast :: IdxPair -> ByteString
getLast (IdxPair v _ ridx) = v ! (ridx - 1)

dropHead :: IdxPair -> IdxPair
dropHead (IdxPair v lidx ridx) = IdxPair v (lidx + 1) ridx

dropLast :: IdxPair -> IdxPair
dropLast (IdxPair v lidx ridx) = IdxPair v lidx (ridx - 1)

getHeadPair :: IdxPair -> IdxPair
getHeadPair (IdxPair v lidx ridx) = IdxPair v lidx (lidx + 1)

ipSize :: IdxPair -> Size
ipSize (IdxPair v lidx ridx) = Size (ridx - lidx)

ipSlice :: IdxPair -> Vector ByteString
ipSlice (IdxPair v lidx ridx) = Vec.slice lidx (ridx - lidx) v

ipTakePrefix :: IdxPair -> Int -> IdxPair
ipTakePrefix (IdxPair v lidx ridx) prf = IdxPair v lidx (lidx + prf)

ipDropPrefix :: IdxPair -> Int -> IdxPair
ipDropPrefix (IdxPair v lidx ridx) prf = IdxPair v (lidx + prf) ridx

ipCombine :: IdxPair -> IdxPair -> IdxPair
ipCombine (IdxPair v lidx1 ridx1) (IdxPair _ lidx2 ridx2) =
  assert (ridx1 == lidx2) (IdxPair v lidx1 ridx2)

ipIsPrefixOf :: IdxPair -> IdxPair -> Bool
ipIsPrefixOf idxPair1 idxPair2 = Vec.length v1 <= Vec.length v2 && v1 == Vec.unsafeTake (Vec.length v1) v2
  where
    v1 = ipSlice idxPair1
    v2 = ipSlice idxPair2

ipSplitAt :: IdxPair -> Int -> (IdxPair, IdxPair)
ipSplitAt idxPair i = (ipTakePrefix idxPair i, ipDropPrefix idxPair i)

data STree
  = Node [(IdxPair, STree)]
  | Leaf
  deriving (Show)

makeTmin :: PartitionType -> Genome -> Genome -> Tmin
makeTmin ptype g h = Tmin ptype (genomeIsSigned g) $ makeHSSTree sTree
  where
    n = genomeSize g
    (l1, _) = interleaveListRepresentation g
    (l2, _) = interleaveListRepresentation h
    (rl1, _) = interleaveListRepresentation (invOri g)
    (rl2, _) = interleaveListRepresentation (invOri h)
    str = Vec.fromList $ case ptype of
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

    -- Construct the Suffix Tree
    sTree = go [0 .. length str - 1]
      where
        go :: [Int] -> STree
        go [] = Leaf
        go ss =
          Node
            [ (IdxPair str (begin_idx -1) end_idx, go ssr)
              | (_, sufs) <- suffixMap ss,
                (begin_idx, end_idx, ssr) <- [findEdge sufs]
            ]

        -- input is a list of suffixes and output the begin and end of the edge infix
        -- (first element and last element + 1) and a list of the subtree suffixes.
        findEdge :: [Int] -> (Int, Int, [Int])
        findEdge [] = error patternError
        findEdge [s] = (s, length str, [])
        findEdge sss@(a_idx : ss)
          | null [c_idx | c_idx <- ss, str ! a_idx /= str ! c_idx] =
            let (_, end_idx, ss') =
                  findEdge ((a_idx + 1) : [c_idx + 1 | c_idx <- filter (/= length str - 1) ss])
             in (a_idx, end_idx, ss')
          | otherwise = (a_idx, a_idx, sss)

        suffixMap :: [Int] -> [(ByteString, [Int])]
        suffixMap = map (second reverse) . Map.toList . List.foldl' step Map.empty
          where
            step m suf_idx = Map.alter (f (suf_idx + 1)) (str ! suf_idx) m
            f i Nothing = Just [i]
            f i (Just is) = Just (i : is)

    -- Convert the suffix tree to HSSTree
    makeHSSTree Leaf = error patternError
    makeHSSTree (Node edges) =
      HSRoot n ((\ts -> assert (fmap getG ts == fmap getH ts) ts) subTrees)
      where
        subTrees = concat $ mapMaybe makeSubHSSTree edges
        makeSubHSSTree e@(idxPair, t) =
          if
              | not $ validBeginIR (getHead idxPair) -> Nothing
              | ipSize idxPair == 1 -> Just $ go False e
              | otherwise ->
                Just $
                  go False (getHeadPair idxPair, Node [(dropHead idxPair, t)])

    go withHead (idxPair_, st) =
      case st of
        Leaf ->
          pure $
            flip Tree.Node [] $
              if
                  | "$" `elem` p__ ->
                    let p = Vec.takeWhile (/= "$") p__
                     in HSLabel
                          (LeafInfo 0 InG LR (getHead idxPair__ /= "$"))
                          (ipTakePrefix idxPair__ (length p))
                  | "#" `elem` p__ ->
                    let p = Vec.takeWhile (/= "#") p__
                     in HSLabel
                          (LeafInfo 0 InH LR (getHead idxPair__ /= "#"))
                          (ipTakePrefix idxPair__ (length p))
                  | "%" `elem` p__ ->
                    let p = Vec.takeWhile (/= "%") p__
                     in HSLabel
                          (LeafInfo 0 InG RL (getHead idxPair__ /= "%"))
                          (ipTakePrefix idxPair__ (length p))
                  | "&" `elem` p__ ->
                    let p = Vec.takeWhile (/= "&") p__
                     in HSLabel
                          (LeafInfo 0 InH RL (getHead idxPair__ /= "&"))
                          (ipTakePrefix idxPair__ (length p))
                  | otherwise -> error patternError
        (Node edges) ->
          if size == 0
            then subTrees
            else [Tree.Node (HSLabel (NodeInfo sG sH) idxPair) subTrees]
          where
            size = ipSize idxPair
            sG = sum . map getG $ subTrees
            sH = sum . map getH $ subTrees
            subTrees = concatMap makeSubHSSTree edges
            makeSubHSSTree = go withHead'
            (withHead', idxPair) =
              ( \h ->
                  if validEndIR h
                    then (False, idxPair__)
                    else (True, dropLast idxPair__)
              )
                $ getLast idxPair__
      where
        p__ = ipSlice idxPair__
        idxPair__ = if withHead then addHead idxPair_ else idxPair_

class WalkDownInfo w where
  goDown :: w -> HSSubTree -> [(w, HSSubTree)]

-- information stored while walking on the tree
data SearchDownInfo = SearchDown
  { subStr :: [IdxPair], -- node's substring
    accPref :: [IdxPair] -- accummulated prefix
  }

idxsToVector :: [IdxPair] -> Vector ByteString
idxsToVector = Vec.concat . map ipSlice . reverse

makeSearchDownInfo :: HSSubTree -> SearchDownInfo
makeSearchDownInfo t@(Tree.Node HSLabel {..} children) =
  if ipSize hsPref >= 2
    then SearchDown [ipTakePrefix hsPref 2] [hsPref]
    else SearchDown [ipTakePrefix hsPref 1] [hsPref]

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
            (ipTakePrefix hsPref 2 : accPref)
            (hsPref : accPref)

data UpdateInfo = Update
  { currentStr :: IdxPair, -- remaining of selected suffix
    breakDistance :: Idx, -- number of characters until the breakpoint is reached
    updateDistance :: Idx, -- number of characters until update starts
    leafInG :: IsInG -- whether the suffix's leaf is in G
  }
  deriving (Show)

makeUpdateInfo :: IdxPair -> Idx -> UpdateInfo
makeUpdateInfo idxPair bDist = Update idxPair bDist 0 InG

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
        if hsPref `ipIsPrefixOf` currentStr
          then (Just $ up {currentStr = ipDropPrefix currentStr (coerce . ipSize $ hsPref)}, t)
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
        )
        [Tree.Node (HSLabel hsInfo {available = False} pr) children]
    NodeInfo {..} ->
      Tree.Node
        (HSLabel (NodeInfo sumG sumH) pl)
        [ Tree.Node
            ( if inG == InG
                then HSLabel (NodeInfo (sumG - 1) sumH) pr
                else HSLabel (NodeInfo sumG (sumH - 1)) pr
            )
            children
        ]
  where
    (pl, pr) = ipSplitAt hsPref (coerce idx - 1)

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
    breakDistance' = breakDistance - (coerce . ipSize $ hsPref)
    (info', node') =
      case hsInfo of
        LeafInfo {..} ->
          let lastBreakDistance' = max (coerce breakDistance) lastBreakDistance
           in ( up
                  { updateDistance = coerce lastBreakDistance - (coerce . ipSize $ hsPref),
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
              ( up {updateDistance = updateDistance - (coerce . ipSize $ hsPref), breakDistance = breakDistance'},
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
    toGenome :: (Vector ByteString, Size, Size, IsInG, Ori) -> (Genome, GenomePosition)
    -- Size is the size of the suffix after the sub-genome occurrence
    -- and the boolean indicate whether the occurrence is in A
    toGenome (bs, _, suf_size, inG, ori) =
      (g,) $
        if inG == InG
          then G idx (rSize - genomeSize g) LR
          else H idx (rSize - genomeSize g) LR
      where
        g_ = interleaveListToGenome (toList bs) sign
        (idx, g) = case ori of
          LR -> (coerce $ rSize - suf_size `div` 2 - genomeSize g + 1, g_)
          RL -> (coerce $ suf_size `div` 2 + 1, invOri g_)

    -- Find minimum element of T'.
    -- Note: a first node is never chosen,
    -- because it can not have a breakpoint (so is never unbalanced)
    walkDown :: SearchDownInfo -> HSSubTree -> Maybe (Vector ByteString, Size, Size, IsInG, Ori)
    walkDown gd@SearchDown {..} currentNode =
      let (Tree.Node HSLabel {..} children) = currentNode
       in case hsInfo of
            LeafInfo {..} ->
              if available
                then
                  let v = idxsToVector subStr
                   in Just (v, Size $ Vec.length v, ipSize hsPref - 2, isInG, isRev)
                else Nothing
            NodeInfo {..} ->
              if sumG /= sumH
                then
                  let inG = if sumG > sumH then InG else InH
                      v = idxsToVector subStr
                   in (\(ori, sp) -> (v, Size $ Vec.length v, sp, inG, ori))
                        . second (subtract 2)
                        <$> sum_pfx inG 0 currentNode
                else safeMinimum . mapMaybe (uncurry walkDown) $ goDown gd currentNode
    sum_pfx :: IsInG -> Size -> HSSubTree -> Maybe (Ori, Size)
    sum_pfx inG acc (Tree.Node HSLabel {..} children) =
      case hsInfo of
        LeafInfo {..} ->
          let acc' = acc + ipSize hsPref
           in if inG == isInG && lastBreakDistance < coerce acc' - 1
                then Just (isRev, acc')
                else Nothing
        NodeInfo {..} -> listToMaybe $ mapMaybe (sum_pfx inG (ipSize hsPref + acc)) children
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
    select gidx x =
      zip (repeat i) . map (\i -> IdxPair v i (Vec.length v)) . take (coerce prefSize) . evens $ [0 ..]
      where
        prefSize = gidx + duoIdx duo - 1
        i = 2 * (coerce n - prefSize) - 1
        v = Vec.fromList . fst . interleaveListRepresentation $ x

    -- search the suffix in each subtree of root
    -- Note: a first node never have a breakpoint
    aux :: Tmin -> (Idx, IdxPair) -> Tmin
    aux (Tmin ptype sign root@HSRoot {..}) (i, s) =
      Tmin ptype sign $ root {rChildren = map aux2 rChildren}
      where
        aux2 t@(Tree.Node l@HSLabel {..} children) =
          if hsPref `ipIsPrefixOf` s
            then snd . go (makeUpdateInfo (ipDropPrefix s (coerce $ ipSize hsPref)) i) $ t
            else t

    -- return the updated node
    go :: UpdateInfo -> HSSubTree -> (Maybe UpdateInfo, HSSubTree)
    go info0@Update {..} t@(Tree.Node l@HSLabel {..} children) =
      case goDown (Just info0) t of
        [] -> case hsInfo of
          LeafInfo {..} ->
            let (info', t'@(Tree.Node l' _)) = updateNode info0 t
             in if inG == isInG && ipSize currentStr == 0 && ori == isRev
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
