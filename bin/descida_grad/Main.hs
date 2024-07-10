{-# LANGUAGE NamedFieldPuns #-}
module Main  where

import qualified Data.Vector as V
import Optimise

data ParamsDescGrad = 
    ParamsDescGrad
        { alpha   :: Double
        , epsilon :: Double
        }

descidaGradienteVet :: 
	ParamsDescGrad -> 
    (V.Vector Double -> Double) ->  -- funcao F
    (V.Vector Double -> V.Vector Double) ->  -- funcao dF/dX
    (V.Vector Double) -- primeiro passo
    -> [Passo (V.Vector Double) Double] -- lista de passos

descidaGradienteVet ParamsDescGrad{ epsilon , alpha } f df x_0=
    let 
        calculePassoSeguinte :: Passo (V.Vector Double) Double -> Passo (V.Vector Double) Double
        calculePassoSeguinte Passo { x = x, y = y} = 
            let dfAtX = df x
                step = V.map (* alpha) dfAtX 
                newX = (-) <$> x <*> step
                newY = f newX
            in 
                Passo{x = newX, y = newY}

        passoSeguinteEhFinal :: Passo (V.Vector Double) Double -> Passo (V.Vector Double) Double -> Bool

        passoSeguinteEhFinal passoAtual passoSeguinte  = 
            let diferencaY = abs  ((y passoAtual) - (y passoSeguinte))
            in  diferencaY < epsilon

        gerePassos passoAnterior =
            let passoAtual = calculePassoSeguinte passoAnterior 
                deveFinalizar = passoSeguinteEhFinal passoAnterior passoAtual
            in if not deveFinalizar 
                then passoAtual : gerePassos passoAtual
                else passoAtual : [] 

        estadoInicial = Passo {x = x_0, y = f x_0}
    in 
        gerePassos estadoInicial 

descidaGradiente 
    :: ParamsDescGrad -> (Double -> Double) -> (Double -> Double) -> Double 
    -> [Passo Double Double]

descidaGradiente params f df x0 = 
    let f_vet = f . single
    	df_vet = V.singleton . df . single
        x0_vet = V.singleton x0
        resultados_descida_vet = descidaGradienteVet params f_vet df_vet x0_vet
        passo_vet_para_passo Passo{x=vx, y=y} = Passo{x = single vx, y= y} 
        resultados_descida = fmap passo_vet_para_passo resultados_descida_vet
    in 
        resultados_descida

main = 
	let 
        f :: Double -> Double
        f x = x^2 + 3*x + 4 

        df :: Double -> Double
        df =    
            let params = ParamDerivNum {delta = 10e-4}
            in  derivNumFn params f

        paramsDescida =  ParamsDescGrad {
            epsilon = 10e-4,
            alpha = 0.01
        }

        x0 = 10.2

        ans_list = descidaGradiente paramsDescida f df x0
    in do
        mapM print ans_list
        return ()

 

