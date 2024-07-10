{-# LANGUAGE NamedFieldPuns , ImportQualifiedPost #-}

module Main where

import qualified Data.Vector as V
import Text.Printf (printf)
import Linear
import Optimise


data FuncoesMetodoNewton =
	FuncoesMetodoNewton
		{ f  :: Double -> Double
		, df :: Double -> Double
		, d2f :: Double -> Double
		}
data FuncoesMetodoNewtonVet =
	FuncoesMetodoNewtonVet
		{ fVet  :: Vetor -> Double
		, dfVet :: Vetor -> Vetor
		, d2fVet :: Vetor -> Matriz
		}
passoMetodoNewtonVet :: FuncoesMetodoNewtonVet -> PassoVet -> PassoVEt 
passoMetodoNewton (FuncoesMetodoNewtonVet{fVet, dfVet, d2fVet}) (Passo{x=xAnterior, y=yAnterior})

passoMetodoNewton :: FuncoesMetodoNewton -> PassoEscalar -> PassoEscalar 
passoMetodoNewton (FuncoesMetodoNewton { f, df, d2f }) (Passo {x=xAnterior, y=yAnterior}) =
	let 
		xAtual = xAnterior - (df xAnterior) / (d2f xAnterior)
		yAtual = f xAtual 
	in 
		Passo{ x=xAtual, y=yAtual}
main = 
    let
        fTeste :: Double -> Double
        fTeste x = x^2 + x * 3 + 4;
        delta = 10e-4;
        dfTeste = derivNumFn delta fTeste;
        d2fTeste = derivNumFn delta dfTeste;
        funcoes = 
            FuncoesMetodoNewton
                { f = fTeste
                , df = dfTeste
                , d2f = d2fTeste
                };

        criterioParada :: PassoEscalar -> PassoEscalar -> Bool;
        criterioParada = deveParar 0.001;

        passoInicial :: PassoEscalar;
        passoInicial = 
            let 
                x = 14.0;
                y = fTeste x;
            in Passo{x, y};

        algoritmo = passoMetodoNewton funcoes;
        resultado = executarPassos criterioParada algoritmo passoInicial;
    in do 
        mapM_ print resultado
        printf "%d passos" (length resultado)