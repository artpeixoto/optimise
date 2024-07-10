{-# LANGUAGE NamedFieldPuns, ImportQualifiedPost #-}
module Optimise where

import Data.Vector qualified as V 
import Linear.Algebra
import Linear.Metric(norm)

pontoAPonto :: (a -> b -> c) -> V.Vector a -> V.Vector b -> V.Vector c
pontoAPonto = V.zipWith 

type ParamDerivNum = Double 
type Vetor = V.Vector Double
type Matriz = V.Vector (V.Vector Double)
type Escalar = Double
        
-- essa funcao a derivada de uma funcao f(x) em que x é um vetor e y é escalar. A derivada é calculada numericamente. 
derivNumFnVet ::  ParamDerivNum -> (Vetor -> Double) -> (Vetor -> Vetor)
derivNumFnVet delta f = 
    let 
        obtenhaDeltas :: Vetor -> [(Vetor, Vetor)]
        obtenhaDeltas x0 = 
            let 
                dimX :: Int
                dimX = V.length x0

                facaParXDelta :: Int -> (Vetor, Vetor) 
                facaParXDelta n = 
                    let 
                        vec_with_delta_at_n = 
                            V.fromList [ (\i -> if i == n then delta / 2 else 0) i  |  i <- [0..dimX] ]
                    in 
                        (  pontoAPonto (+)  x0 vec_with_delta_at_n 
                        ,  pontoAPonto (-)  x0 vec_with_delta_at_n
                        )
            in 
                [ facaParXDelta x | x <- [0..dimX] ]

        calculeDySobreDx :: (Vetor, Vetor) -> Double
        calculeDySobreDx (x_plus_delta, x_minus_delta) = 
            ((f x_plus_delta) -  (f x_minus_delta)) / delta 
        in  
            \x -> x |> obtenhaDeltas |> map clculeDySobreDx |> V.fromList 

derivNumFnMat :: ParamDerivNum -> (Vetor -> Vetor) -> (Vetor -> Matriz)
derivNumFnMat delta fVet xVet = 
    let 
        obtenhaDeltas :: Vetor -> [(Vetor, Vetor)]
        obtenhaDeltas x0 = 
            let 
                dimX :: Int
                dimX = V.length x0

                facaParXDelta :: Int -> (Vetor, Vetor) 
                facaParXDelta n = 
                    let 
                        vec_with_delta_at_n = 
                            V.fromList [ (\i -> if i == n then delta / 2 else 0) i  |  i <- [0..dimX] ]
                    in 
                        (  pontoAPonto (+)  x0 vec_with_delta_at_n 
                        ,  pontoAPonto (-)  x0 vec_with_delta_at_n
                        )
            in 
                [ facaParXDelta x | x <- [0..dimX] ]

        calculeDySobreDx :: (Vetor, Vetor) -> Vetor
        calculeDySobreDx (xMaisDelta, xMenosDelta) = 
            (fVet xMaisDelta !-! fVet xMenosDelta) 
            |> fmap (/delta)
    in 
        \xVet -> V.fromList . map calculeDySobreDx . obtenhaDeltas $ xVet

single :: V.Vector a -> a
single v = if V.length v == 1 then v V.! 0 else error "Vetor deve ter apenas um elemento"

derivNumFn :: ParamDerivNum -> (Double -> Double) -> (Double -> Double) 
derivNumFn params f = 
    let f_vet v = f $ single v 
        df_vet = derivNumFnVet params f_vet
    in \x -> (V.! 0) . df_vet . V.singleton $ x

data Passo a b = 
    Passo 
        {   x :: a
        ,   y :: b
        }
    deriving (Show, Eq)

type PassoEscalar = Passo Double Double
type PassoVet     = Passo (Vetor) Double

infixl 0 |>
(|>) :: a -> (a -> b) -> b
x |> f = f x


deveParar :: Double -> Passo Double a  -> Passo Double a -> Bool
deveParar epsilon Passo {x=xAnterior, y=_ } Passo {x=xAtual, y=_} = 
    ( xAnterior - xAtual
    |> abs 
    |> \diff -> diff <= epsilon
    )

devePararVet :: Double -> Passo Vetor a -> PassoVetor a -> Bool

devePararVet epsilon Passo{x=xAnterior, y=_} Passo{x=xAtual, y=_} = 
    let distanciaAndada = xAnterior !-! xAtual |> norm 
    in  distanciaAndada <= epsilon
 


executarPassos :: (Passo a b -> Passo a b -> Bool) -> ((Passo a b) -> (Passo a b)) -> Passo a b -> [Passo a b]

executarPassos criterioParada algoritmo passoInicial = 
    let obtenhaPassosSeguintes passoAnterior =
            let passoAtual = algoritmo passoAnterior
            in if criterioParada passoAnterior passoAtual 
                then passoAtual : []
                else passoAtual : (obtenhaPassosSeguintes passoAtual)
    in passoInicial : obtenhaPassosSeguintes passoInicial
