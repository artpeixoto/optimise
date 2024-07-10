{-# LANGUAGE NamedFieldPuns #-}
module Main where
import Optimise
import qualified Data.Vector as V

data FuncoesMetodoNewton =
	FuncoesMetodoNewton
		{ f  		:: Double -> Double
		, df 		:: Double -> Double
		, d2f		:: Double -> Double
		}

passoMetodoNewton :: FuncoesMetodoNewton -> PassoEscalar -> PassoEscalar 
passoMetodoNewton FuncoesMetodoNewton{ f, df, d2f } passoAnterior =
	let xAnterior 	= x passoAnterior 
		xAtual 		= xAnterior - (df xAnterior) / (d2f xAnterior)
		yAtual 		= 
	in Passo


  
main :: IO ()
main = 
	let f :: Double -> Double
        f x = x^2 + 3*x + 4 

		paramsDeriv = ParamDerivNum {delta = 10e-4}

        df  = derivNumFn params f
		d2f = derivNumFn params df

		epsilon = 10e-3
		criterioParada = deveParar epsilon 

		passoInicial = 
			let xInicial = 5
				yInicial = f 5
			in Passo{x=xInicial, y=yInicial}

		passos = executarPassos criterioParada passoMetodoNewton passoInicial  
	in 
		






