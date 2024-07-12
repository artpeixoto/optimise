use nalgebra::DMatrix;

use crate::{ferramentas::{calculos::diff::{self, escalar::derivada_numerica}, exemplos}, monte_otimizacao, AlgoritmoOtimizacao, CriterioDeltaXMenorQueEpsilon, Matriz, Passo, PassoEscalar, Vetor};


#[test]
fn exemplo_vet(){
	let epsilon = 10e-5;
	let criterio_parada = CriterioDeltaXMenorQueEpsilon(epsilon);

	let (f_exemplo, passo_inicial) = exemplos::vetorial::obtenha_f_e_passo_inicial();
	let df_exemplo = diff::vetorial::derivada_numerica(10e-6, &f_exemplo);
	let primeira_hessiana  = diff::matriz::derivada_numerica(10e-6, &df_exemplo)(&passo_inicial.x);

	let mut algoritmo = MetodoQuaseNewton{
		f: &f_exemplo, 
		df: &df_exemplo,
		alfa: 0.5,	
		inv_hess: primeira_hessiana.try_inverse().unwrap(),
		df_x_anterior : df_exemplo(&passo_inicial.x)
	};
	let otimizacao = monte_otimizacao(algoritmo, criterio_parada, passo_inicial);
	otimizacao.enumerate().for_each(|(n_passo, passo)| {println!("{:000} - {:?}", n_passo, passo); });
}
pub struct MetodoQuaseNewton< F, Df> 
	where F: Fn(&Vetor) -> f64, Df: Fn(&Vetor) -> Vetor 
{
	pub f				: F,
	pub df				: Df,
	pub alfa			: f64,
	pub inv_hess 		: Matriz,
	pub df_x_anterior	: Vetor,
}

impl<F, DF> AlgoritmoOtimizacao<Vetor, f64> for MetodoQuaseNewton< F, DF> 
where 
	F : Fn(&Vetor) -> f64,
	DF: Fn(&Vetor) -> Vetor,
{
	fn proximo_passo(&mut self, passo_anterior: &Passo<Vetor, f64>) -> Passo<Vetor, f64> {
		let delta_x = - self.alfa * self.inv_hess.clone() * (self.df)(&passo_anterior.x);	
		let x = passo_anterior.x.clone() + delta_x.clone();
		let y = (self.f)(&x);
		let df_x = (self.df)(&x);
		self.inv_hess = { // calculamos uma estimativa da hessiana usando o metodo de broyden
			let yk: Vetor =  df_x.clone() - self.df_x_anterior.clone() ;
			let ident_matrix = DMatrix::<f64>::identity(yk.len(), yk.len());
			let yk_t_delta_x = (yk.transpose() * delta_x.clone())[0];
			let inv_hess = 
				( 	( 	( ident_matrix.clone() - ( (delta_x.clone() * yk.transpose()) / (yk_t_delta_x)) ) 
					* 	self.inv_hess.clone()
					* 	( ident_matrix.clone() - ( (yk * delta_x.clone().transpose()) / (yk_t_delta_x) ) ) 
					)
					+ ((delta_x.clone() * delta_x.clone().transpose()) / yk_t_delta_x )
				);
			inv_hess
		};
		self.df_x_anterior = df_x;
		Passo{x, y}
	}
}


