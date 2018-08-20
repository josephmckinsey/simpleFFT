{-# LANGUAGE ParallelListComp #-}

module Lib
    ( showC,
      printC,
      bitReverseMap,
      generationFFT,
      generationIFFT,
      decimationFFT,
      decimationIFFT
    ) where

import Data.Complex

-- |Shows the list of complex numbers in a nice way, without reimplementing show itself.
showC :: Show a => [Complex a] -> String
showC [] = "[]"
showC xs = "[\n" ++ (unlines $ map (\i -> "  " ++ show i ++ ",") xs) ++ "]\n"

-- |Print complex numbers easily.
printC :: Show a => [Complex a] -> IO ()
printC xs = putStrLn $ showC xs

-- |Moves the last item to the front. O(n)
-- (a, b, c, d) => (d, a, b, c)
translateBack :: Num a => [Complex a] -> [Complex a]
translateBack [] = []
translateBack xs = last xs : init xs

-- |Fourier transform of translating back. It's just multiplying by e^{-2 \pi i
-- x / N} This will be the meat of the forward calculation.
--
-- (a, b, c, d) => (a, omega * b, omega^2 * c, omega^3 * d)
modulateForward :: RealFloat a => Int -> [Complex a] -> [Complex a]
modulateForward n xs =
    let twiddles = map (\i -> cis (-2 * pi * fromIntegral i / fromIntegral n)) [0..] in
        [w * i | i <- xs | w <- twiddles]

-- |Fourier transform of translating forward. This will be the meat of the
-- inverse calculation. Uses +2 instead of -2 as the coefficient. Equivalent to
-- moving forward one.
--
-- (a, b, c, d) => (a, omega * b, omega^2 * c, omega^3 * d)
modulateReverse :: RealFloat a => Int -> [Complex a] -> [Complex a]
modulateReverse n xs =
    let twiddles = map (\i -> cis (2 * pi * fromIntegral i / fromIntegral n)) [0..] in
        [w * i | i <- xs | w <- twiddles]

-- |Pad adds a zero between every other zero. It's Fourier transform is decompress.
--
-- (a, b, c, d) => (a, 0, b, 0, c, 0, d, 0)
pad :: Num a => [Complex a] -> [Complex a]
pad [] = []
pad (x:xs) = x:(0 :+ 0):(pad xs)

-- |Opposite of pad. Takes adjacent numbers and adds them together.
--
-- (a, b, c, d) => (1/2) * (a + b, c + d)
unpad :: RealFloat a => [Complex a] -> [Complex a]
unpad [] = []
unpad (x:[]) = [x]
unpad (x:y:xs) = (x + y):xs

-- |Fourier transform of pad. Unfolds the list, and divides by two.
--
-- (A, B, C, D) => (1/2) * (A, B, C, D, A, B, C, D)
decompress :: RealFloat a => [Complex a] -> [Complex a]
decompress xs = map (/ 2) (xs ++ xs)

-- |Fourier transform of the unpad, where you remove "zeros", by
-- squishing adjacent values. But fold would conflict with the standard Haskell
-- syntax.
--
-- (A, B, C, D) => (1/2) * (A + C, B + D)
compress :: RealFloat a => Int -> [Complex a] -> [Complex a]
compress n xs =
    let (a, b) = splitAt (div n 2) xs in
        map (/2) $ zipWith (+) a b

-- | Inverse version without 1/2 so that it adds to 1/N at the end.
--
-- (A, B, C, D) => (A + C, B + D)
compress' :: RealFloat a => Int -> [Complex a] -> [Complex a]
compress' n xs =
    let (a, b) = splitAt (div n 2) xs in
        zipWith (+) a b

-- |Small helper to simplify generator function.
--
-- (A, B, C) + (a, b, c) => (A + a, B + b, C + c)
add :: Num a => [a] -> [a] -> [a]
add x y = zipWith (+) x y

-- |One unit of a generative fast Fourier transform. See book for more details
-- and explanation.
generationUnit :: RealFloat a => Int -> [Complex a] -> [Complex a] -> [Complex a]
generationUnit n x y =
    (decompress x) `add` (modulateForward n . decompress $ y)


-- |Inverse (note that this version is basically off by a sign and a constant)
-- |Note that (x ++ x) is used instead of decompress since it only lacks a 2.
generationIUnit :: RealFloat a => Int -> [Complex a] -> [Complex a] -> [Complex a]
generationIUnit n x y =
    (x ++ x) `add` (modulateReverse n $ y ++ y) -- uses ++ so it avoids 1/2 on
                                                -- fold

-- |Decimation version. See book for more details and explanation.
decimationUnit :: RealFloat a => Int -> [Complex a] -> ([Complex a], [Complex a])
decimationUnit n x = (compress n x, compress n . modulateForward n $ x)

-- |Inverse operation off by a factor and a sign. Uses compress' which doesn't
-- have the 2.
decimationIUnit :: RealFloat a => Int -> [Complex a] -> ([Complex a], [Complex a])
decimationIUnit n x = (compress' n x, compress' n . modulateReverse n $ x)

-- |Nice little feature in case you try to make this imperative (or remove the
-- alternating costs more likely). Bit reverse doesn't improve the algorithmic
-- complexity in this case, so I did not include it in the base program.
--
-- It is based on the "Bracewell-Buneman algorithm", but Haskell's linked lists
-- may make it algorithmically worse than an imperative version.
--
-- 4 => (0, 3, 1, 2)
bitReverseMap :: Int -> [Int]
bitReverseMap 1 = [0]
bitReverseMap n =
    let x = map (*2) $ bitReverseMap (div n 2) in
        x ++ map (+1) x

-- |Alternate is the little cousin of a bit reverse. You take too lists and
-- alternate between the elements. This is used in the decimationFFT to combine
-- the two results into one.
--
-- (a, b), (c, d) => (a, c, b, d)
alternate :: [a] -> [a] -> [a]
alternate [] _ = []
alternate _ [] = []
alternate (x:xs) (y:ys) = x:y:(alternate xs ys)

-- |Unalternate is the inverse of the alternate, where it splits up the list the
-- exact same way. This is used in the generationFFT to feed the recursive FFTs.
--
-- (a, b, c, d) => ((a, c), (b, d))
unalternate :: [a] -> ([a], [a])
unalternate [] = ([], [])
unalternate (x:[]) = ([x], [])
unalternate (x:y:xs) = (x:a, y:b)
    where (a, b) = unalternate xs

-- |This puts all the above functions together in a assembly-based Fourier
-- transform. Very elegant! Note that I don't do any checks on n-size, and it
-- fail unexpectedly and silently for those values, since the Fourier transforms
-- depend on it being divisible by 2.
generationFFT :: RealFloat a => Int -> [Complex a] -> [Complex a]
generationFFT 1 x = x
generationFFT n x =
    generationUnit n (generationFFT halfn a) (generationFFT halfn b)
    where (a, b) = unalternate x
          halfn = div n 2

-- |Inverse of the above (really just the same with name changes). (off by a
-- factor and a sign in the process).
generationIFFT :: RealFloat a => Int -> [Complex a] -> [Complex a]
generationIFFT 1 x = x
generationIFFT n x =
    generationIUnit n (generationIFFT halfn a) (generationIFFT halfn b)
    where (a, b) = unalternate x
          halfn = div n 2

-- |Disassembly version of the above. Note how similar it is to the version
-- above, but with a small difference in how it proceeds in computation. This
-- one breaks up lists into smaller values, while the generative goes
-- recursively down first.
decimationFFT :: RealFloat a => Int -> [Complex a] -> [Complex a]
decimationFFT 1 x = x
decimationFFT n x =
    alternate (decimationFFT halfn a) (decimationFFT halfn b)
    where (a, b) = decimationUnit n x
          halfn = div n 2

-- |Decimation inverse. Off by a constant and sign in the process.
decimationIFFT :: RealFloat a => Int -> [Complex a] -> [Complex a]
decimationIFFT 1 x = x
decimationIFFT n x =
    alternate (decimationIFFT halfn a) (decimationIFFT halfn b)
    where (a, b) = decimationIUnit n x
          halfn = div n 2

-- TODO: Add testing features to show that it works.
