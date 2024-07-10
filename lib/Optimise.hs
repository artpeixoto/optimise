{-# LANGUAGE NamedFieldPuns #-}
module Optimise where
import qualified Data.Vector as V

-- a função abaixo transforma uma função de duas variaveis 
lift :: (a -> b -> c) -> V.Vector a -> V.Vector b -> V.Vector c
lift = V.zipWith 

data ParamDerivNum = 
    ParamDerivNum
    {   delta :: Double 
    }

-- essa funcao a derivada de uma funcao f(x) em que x é um vetor e y é escalar. A derivada é calculada numericamente. 
derivNumFnVet ::  ParamDerivNum -> (V.Vector Double -> Double) -> (V.Vector Double -> V.Vector Double)

derivNumFnVet ParamDerivNum{delta} f = 
    let 
        obtenhaDeltas :: V.Vector Double -> [(V.Vector Double, V.Vector Double)]
        obtenhaDeltas x_0 = 
            let 
                dimX :: Int
                dimX = V.length x_0

                facaParXDelta :: Int -> (V.Vector Double, V.Vector Double) 
                facaParXDelta n = 
                    let 
                        vec_with_delta_at_n = 
                            V.fromList [ (\i -> if i == n then delta / 2 else 0) i  |  i <- [0..dimX] ]
                    in 
                        (  lift (+)  x_0 vec_with_delta_at_n 
                        ,  lift (-)  x_0 vec_with_delta_at_n
                        )
            in 
                [ facaParXDelta x | x <- [0..dimX] ]

        calculeDySobreDx :: (V.Vector Double, V.Vector Double) -> Double
        calculeDySobreDx (x_plus_delta, x_minus_delta) = 
            ((f x_plus_delta) -  (f x_minus_delta)) / delta 
        in  
            \x -> V.fromList . map calculeDySobreDx . obtenhaDeltas $ x  

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

infixl 0 |>
(|>) :: a -> (a -> b) -> b
x |> f = f x

type Epsilon = Double ;

deveParar :: Epsilon -> Passo a Double -> Passo a Double -> Bool
deveParar epsilon Passo {x=_, y=yAnterior} Passo {x=_, y=yAtual} = 
    ( yAtual - yAnterior 
    |> abs 
    |> \diff -> diff <= epsilon
    )


executarPassos :: (Passo a Double -> Passo a Double -> Bool) -> ((Passo a Double) -> (Passo a Double)) -> Passo a Double -> [Passo a Double]

executarPassos criterioParada algoritmo passoInicial = 
    let obtenhaPassosSeguintes passoAnterior =
            let passoAtual = algoritmo passoAnterior
            in if criterioParada passoAnterior passoAtual 
                then passoAtual : []
                else passoAtual : (obtenhaPassosSeguintes passoAtual)
    in passoInicial : obtenhaPassosSeguintes passoInicial
