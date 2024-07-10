{-# LANGUAGE NamedFieldPuns, ImportQualifiedPost #-}
module Main where

import qualified Data.Vector as V
import Optimise
import Text.Printf (printf)

data ParamsDescGrad = 
    ParamsDescGrad
        { alpha   :: Double
        }

data FuncoesDescidaGrad = 
    FuncoesDescidaGrad 
        { f  :: Double -> Double
        , df :: Double -> Double
        }

data FuncoesDescidaGradVet = 
    FuncoesDescidaGradVet 
        { fVet :: Vetor -> Double
        , dfVet:: Vetor -> Vetor
        }


passoDescidaGradVet  :: ParamsDescGrad -> FuncoesDescidaGradVet -> PassoVet -> PassoVet 
passoDescidaGradVet ParamsDescGrad {alpha} FuncoesDescidaGradVet{ fVet, dfVet } Passo{x=xAnterior, y=yAnterior} = 
    let 
        gradiente = dfVet xAnterior
        deslocamentoX = V.map (*alpha) gradiente
        xAtual = pontoAPonto (-) xAnterior deslocamentoX 
        yAtual = fVet xAtual
    in Passo {
        x = xAtual,
        y = yAtual
    }

passoDescidaGrad :: ParamsDescGrad -> FuncoesDescidaGrad -> PassoEscalar -> PassoEscalar
passoDescidaGrad params FuncoesDescidaGrad{f, df} Passo{x=xAnterior, y=yAnterior} =
    let 
        fVet  = f . single
        dfVet = V.singleton . df . single 
        xVetAnterior  = V.singleton xAnterior
        passoVet = passoDescidaGradVet params FuncoesDescidaGradVet{fVet, dfVet} Passo{x=xVetAnterior, y=yAnterior}
        passo = let Passo{x=xAtualVet, y=yAtual} = passoVet  in Passo{x=single xAtualVet, y=yAtual}
    in 
        passo

    
main = 
	let 
        f :: Double -> Double
        f x = x^2 + 3*x + 4 

        df :: Double -> Double
        df =    
            let params =  10e-4
            in  derivNumFn params f

        paramsDescida =  ParamsDescGrad {
            alpha = 0.01
        }

        passoInicial = let x = 5 in Passo{x=x, y = f x}
        
        criterioParada = deveParar 1e-3

        algoritmo = passoDescidaGrad paramsDescida FuncoesDescidaGrad{ f, df } 

        resultado = executarPassos criterioParada algoritmo passoInicial
    in do
        mapM_ print resultado
        printf "%d passos" (length resultado)
