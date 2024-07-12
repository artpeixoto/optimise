use std::iter;
use nalgebra::DMatrix;
use replace_with::replace_with_or_abort_and_return;

pub type Vetor = nalgebra::DVector<f64>;
pub type Matriz = nalgebra::DMatrix<f64>;

#[derive(Clone, Copy, PartialEq, Eq, Debug,Hash)]
pub struct Passo<X, Y> {
	pub x: X,
	pub y: Y
}

pub type PassoEscalar 				= Passo<f64, f64>;
pub type PassoVet = Passo<Vetor, f64>;   

pub trait CriterioParada<X, Y>{
	fn deve_parar(&self, passo_anterior: &Passo<X, Y> , passo_seguinte: &Passo<X, Y>) -> bool;
}

pub struct CriterioDeltaXMenorQueEpsilon(pub f64);
impl CriterioParada<f64, f64> for CriterioDeltaXMenorQueEpsilon		{
	fn deve_parar(&self, passo_anterior: &Passo<f64, f64> , passo_seguinte: &Passo<f64, f64>) -> bool {
		(passo_anterior.x - passo_seguinte.x).abs() <= self.0
	}
}
impl CriterioParada<Vetor, f64> for CriterioDeltaXMenorQueEpsilon{
	fn deve_parar(&self, passo_anterior: &Passo<Vetor, f64> , passo_seguinte: &Passo<Vetor, f64>) -> bool {
		let delta_x = passo_seguinte.x.clone() - passo_anterior.x.clone();
		delta_x.norm() <= self.0
	}
}

pub trait AlgoritmoOtimizacao<X,Y>{
	fn proximo_passo(&mut self, passo: &Passo<X, Y>) -> Passo<X, Y>;
}

pub fn monte_otimizacao<X, Y, Algo, Crit>
	(mut algo: Algo, crit: Crit, passo_inicial: Passo<X, Y>) -> impl Iterator<Item=Passo<X, Y>> 
	where 
		Algo: AlgoritmoOtimizacao<X,Y>,
		Crit: CriterioParada<X, Y>
{
	#[derive(Clone, PartialEq, Eq)]
	enum EstadoExecOtimizacao<X, Y>{
		Executando(Passo<X,Y>),
		Finalizando(Passo<X,Y>),
		Acabado,
	}
	let mut estado_exec = EstadoExecOtimizacao::Executando(passo_inicial);
	
	iter::from_fn(move || replace_with_or_abort_and_return(&mut estado_exec, |estado_exec: EstadoExecOtimizacao<X, Y>| 
		match estado_exec{
			EstadoExecOtimizacao::Executando(passo_guardado) => {
				let novo_passo = algo.proximo_passo(&passo_guardado);  //calculamos o proximo
				let proximo_estado_exec = if crit.deve_parar(&passo_guardado, &novo_passo){
					EstadoExecOtimizacao::Finalizando(novo_passo)
				} else {
					EstadoExecOtimizacao::Executando(novo_passo)
				};
				(Some(passo_guardado), proximo_estado_exec)
			}
			EstadoExecOtimizacao::Finalizando(passo_guardado) => (Some(passo_guardado), EstadoExecOtimizacao::Acabado),
			EstadoExecOtimizacao::Acabado => (None, EstadoExecOtimizacao::Acabado),
		}
	))
}