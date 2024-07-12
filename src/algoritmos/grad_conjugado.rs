use crate::{ferramentas::{calculos::diff, exemplos}, monte_otimizacao, AlgoritmoOtimizacao, CriterioDeltaXMenorQueEpsilon, Matriz, Passo, Vetor};


pub struct GradConjugado<F, Df>{
	pub f : F,
	pub df: Df,
	pub Q : Matriz,
	pub d : Vetor
}

impl<F, Df> AlgoritmoOtimizacao<Vetor, f64> for GradConjugado<F, Df> where F: Fn(&Vetor) -> f64, Df: Fn(&Vetor) -> Vetor{
	fn proximo_passo(&mut self, passo_anterior: &Passo<Vetor, f64>) -> Passo<Vetor, f64> {
		let q_d = self.Q.clone() * self.d.clone();
		let q_d_dot_d = q_d.dot(&self.d);
		let alfa =     q_d_dot_d / q_d.dot(&q_d);
		let x_novo = passo_anterior.x.clone() + alfa * self.d.clone();
		let dy = (self.df)(&x_novo);
		let beta = q_d.dot(&dy) / q_d_dot_d;
		let d_novo = -dy + beta * self.d.clone();
		self.d = d_novo;
		let y_novo = (self.f)(&x_novo);
	 	Passo{x:x_novo, y:y_novo}
	}
}


#[test]
fn exemplo_vet(){
	let epsilon = 10e-5;
	let criterio_parada = CriterioDeltaXMenorQueEpsilon(epsilon);

	let (f_exemplo, passo_inicial) = exemplos::vetorial::obtenha_f_e_passo_inicial();
	let df_exemplo = diff::vetorial::derivada_numerica(10e-6, &f_exemplo);
	let matriz_hessiana_diagonal =  {
		let hessiana: Matriz = (diff::matriz::derivada_numerica(10e-6, &df_exemplo)(&passo_inicial.x));
		hessiana
	};
	let algoritmo = GradConjugado{
		f: &f_exemplo, 
		df: &df_exemplo,
		Q: matriz_hessiana_diagonal,	
		d: df_exemplo(&passo_inicial.x)
	};
	let otimizacao = monte_otimizacao(algoritmo, criterio_parada, passo_inicial);
	otimizacao.enumerate().for_each(|(n_passo, passo)| {println!("{:000} - {:?}", n_passo, passo); });
}