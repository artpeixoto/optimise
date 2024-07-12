use crate::{ferramentas::{calculos::{self, diff}, exemplos}, monte_otimizacao, AlgoritmoOtimizacao, CriterioDeltaXMenorQueEpsilon, Passo, PassoEscalar, Vetor};

pub struct DescidaGrad<F, Df>{
	pub f : F,
	pub df : Df,
	pub alpha: f64,
}


impl<F, Df> AlgoritmoOtimizacao<f64, f64> for DescidaGrad<F, Df> 
	where F: Fn(&f64) -> f64, Df: Fn(&f64) -> f64
{
	fn proximo_passo(&mut self, passo_anterior: &Passo<f64, f64>) -> Passo<f64, f64> {
		let novo_x = passo_anterior.x - (self.alpha *  (self.df)(&passo_anterior.x));
		let novo_y = (self.f)(&novo_x);
		Passo{x: novo_x, y: novo_y}
	}
}

impl<F, Df> AlgoritmoOtimizacao<Vetor, f64> for DescidaGrad<F, Df> where F: Fn(&Vetor) -> f64, Df: Fn(&Vetor) -> Vetor{
	fn proximo_passo(&mut self, passo_anterior: &Passo<Vetor, f64>) -> Passo<Vetor, f64> {
		let novo_x = passo_anterior.x.clone() - (self.alpha * (self.df)(&passo_anterior.x));
		let novo_y = (self.f)(&novo_x);
		Passo{x: novo_x, y: novo_y}
	}
}

#[test]
fn exemplo_escalar(){
	let epsilon = 10e-5;
	let criterio_parada = CriterioDeltaXMenorQueEpsilon(epsilon);
	let f_exemplo = |x:&f64| x.powf(2.0) + 3.0*x + 4.0 ;
	let df_exemplo = calculos::diff::escalar::derivada_numerica(10e-5, f_exemplo );
	let mut algoritmo = DescidaGrad{
		f: f_exemplo.clone(), 
		df: df_exemplo,
		alpha: 0.01,	
	};
	let x_inicial = 15.0;
	let y_inicial = f_exemplo(&x_inicial);
	let passo_inicial = PassoEscalar{x: x_inicial, y: y_inicial};
	let otimizacao = monte_otimizacao(algoritmo, criterio_parada, passo_inicial);
	otimizacao.enumerate().for_each(|(n_passo, passo)| {println!("{:000} - {:?}", n_passo, passo); });
}
#[test]
fn exemplo_vet(){
	let epsilon = 10e-5;
	let criterio_parada = CriterioDeltaXMenorQueEpsilon(epsilon);

	let (f_exemplo, passo_inicial) = exemplos::vetorial::obtenha_f_e_passo_inicial();

	let df_exemplo = diff::vetorial::derivada_numerica(10e-6, &f_exemplo);
	let mut algoritmo = DescidaGrad{
		f: &f_exemplo, 
		df: &df_exemplo,
		alpha: 0.01,	
	};
	let otimizacao = monte_otimizacao(algoritmo, criterio_parada, passo_inicial);
	otimizacao.enumerate().for_each(|(n_passo, passo)| {println!("{:000} - {:?}", n_passo, passo); });
}
