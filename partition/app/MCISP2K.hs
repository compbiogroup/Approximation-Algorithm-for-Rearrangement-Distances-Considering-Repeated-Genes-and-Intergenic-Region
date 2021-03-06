{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      : MCISP2K
-- Description : Algorithm for Genome partition problems with approximation factor 2k (k in the highest occurrence of a gene in the genome).
-- Copyright   : (c) Gabriel Siqueira, 2021
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module MCISP2K where

import Control.Concurrent.ParallelIO.Global (parallel, stopGlobalPool)
import Control.DeepSeq (force)
import qualified Data.ByteString.Char8 as BS
import Data.Time
import Genomes (Genome, Sign(..), readGenome, writeGenome)
import LocalBase
import Options.Applicative
import Partition (getPartition, reduced, PartitionType(..))

data Args = Args
  { input :: String,
    output :: String,
    noParallel :: Bool,
    signed :: Sign,
    partType :: PartitionType
  }

argsParser :: Parser Args
argsParser =
  Args
    <$> strOption
      ( long "input"
          <> short 'i'
          <> metavar "IFILE"
          <> help "Input file. Each 4 lines of the input file correspond to a instance, each line has a list of comma or space separated values, and represent in order the origin string, the origin intergenic region list, the target string, and the target intergenic region list."
      )
      <*> strOption
        ( long "outfile"
            <> short 'o'
            <> metavar "OFILE"
            <> help "Output file. For each instance five lines are produces in the file, the first four lines correspond to two reduced genomes produced with the partition algorithm (each character of the string correspond to a block of the partition and each integer of the intergenic regions list correspond to a breakpoint). The last line shows the wall clock time required to produce the partition."
        )
      <*> switch
        ( long "no-par"
            <> help "Do not process the genomes in parallel."
        )
      <*> flag Unsigned Signed 
        ( long "signed"
            <> short 's'
            <> help "Whether the input Strings are signed."
        )
      <*> flag MCISP RMCISP
        ( long "rev"
            <> help "Whether to use reverse partition."
        )

opts :: ParserInfo Args
opts =
  info
    (argsParser <**> helper)
    ( fullDesc
        <> progDesc
          "Algorithm for genome partition problems."
    )

main :: IO ()
main = do
  args <- execParser opts
  contents <- BS.readFile (input args)
  let pairs = toQuadruples . filter ((/= '#') . BS.head) . BS.lines $ contents
  ans <-
    if noParallel args
      then mapM (runOne args) pairs
      else do
        ans <- parallel $ map (runOne args) pairs
        stopGlobalPool
        return ans
  BS.writeFile (output args) . BS.unlines . fromAns $ ans
  where
    runOne args bstrs = do
      start <- getCurrentTime
      let !bstrs' = force $ simplifyGenomes (partType args) (signed args) bstrs
      end <- getCurrentTime
      let time = BS.pack . show . realToFrac $ diffUTCTime end start
      return (bstrs', "# Time: " <> (BS.pack . show $ time))

    toQuadruples (s1 : i1 : s2 : i2 : ss) = (s1, i1, s2, i2) : toQuadruples ss
    toQuadruples [] = []
    toQuadruples _ = error "Incorrect number of lines."

    fromAns (((s1, i1, s2, i2), time) : ss) = s1 : i1 : s2 : i2 : time : fromAns ss
    fromAns [] = []

simplifyGenomes :: PartitionType -> Sign -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString) -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString)
simplifyGenomes ptype signed (s1, i1, s2, i2) = (s1', i1', s2', i2')
  where
    (s1', i1') = writeGenome False g'
    (s2', i2') = writeGenome False h'
    (g', h') = reduced part
    part = getPartition ptype g h
    g = readGenome True signed s1 i1
    h = readGenome True signed s2 i2
