use crate::{ferramentas::{calculos::diff::{escalar, matriz, vetorial}, exemplos}, monte_otimizacao, AlgoritmoOtimizacao, CriterioDeltaXMenorQueEpsilon, Matriz, Passo, Vetor};




pub struct MetodoNewton<F, DF, D2F>{
	f: F,
	df: DF,
	d2f: D2F 	
}

impl<F, DF, D2F> AlgoritmoOtimizacao<f64, f64> for MetodoNewton<F, DF, D2F> 
where 
	F	: Fn(&f64) -> f64,
	DF	: Fn(&f64) -> f64,
	D2F	: Fn(&f64) -> f64,
{
	fn proximo_passo(&mut self, passo: &Passo<f64, f64>) -> Passo<f64, f64> {
		let proximo_x = passo.x - ((self.df)(&passo.x) / (self.d2f)(&passo.x));
		let proximo_y = (self.f)(&proximo_x);
		Passo{x:proximo_x, y: proximo_y}
	}
}
impl< F, DF, D2F> AlgoritmoOtimizacao<Vetor, f64> for MetodoNewton<F, DF, D2F> 
where 
	F :  Fn(&Vetor) -> f64,
	DF:  Fn(&Vetor) -> Vetor,
	D2F: Fn(&Vetor) -> Matriz,
{
	fn proximo_passo(&mut self, passo: &Passo<Vetor, f64>) -> Passo<Vetor, f64> {
		let hessiana = (self.d2f)(&passo.x);
		let inverso_hessiana = hessiana.try_inverse().expect("a matriz hessiana não pode ser invertida. Isso é bem calvo da sua parte");
		let proximo_x = passo.x.clone() - (inverso_hessiana * (self.df)(&passo.x));
		let proximo_y = (self.f)(&proximo_x);
		Passo{x:proximo_x, y: proximo_y}
	}
}


	#[test]
	fn exemplo_vet() {
		let (f_exemplo, passo_inicial) = exemplos::vetorial::obtenha_f_e_passo_inicial();
		let df_exemplo = vetorial::derivada_numerica(10e-4, &f_exemplo);
		let d2f_exemplo = matriz::derivada_numerica(10e-4,  &df_exemplo);

		let epsilon = 10e-5;
		let criterio_parada = CriterioDeltaXMenorQueEpsilon(epsilon);

		let mut algoritmo = MetodoNewton{
			f: &f_exemplo, 
			df: &df_exemplo,
			d2f: &d2f_exemplo,
		};

		let otimizacao = monte_otimizacao(algoritmo, criterio_parada, passo_inicial);
		for passo_otim in otimizacao{
			println!("passo: {passo_otim:?}")
		}
	}
	#[test]
	fn exemplo_escalar() {
		let (f_exemplo, passo_inicial) = exemplos::escalar::obtenha_f_e_passo_inicial();
		let df_exemplo  = escalar::derivada_numerica(10e-4, &f_exemplo);
		let d2f_exemplo = escalar::derivada_numerica(10e-4, &df_exemplo );

		let epsilon = 10e-5;
		let criterio_parada = CriterioDeltaXMenorQueEpsilon(epsilon);

		let mut algoritmo = MetodoNewton{
			f: &f_exemplo, 
			df: &df_exemplo,
			d2f: &d2f_exemplo,
		};

		let otimizacao = monte_otimizacao(algoritmo, criterio_parada, passo_inicial);
		for passo_otim in otimizacao{
			println!("passo: {passo_otim:?}")
		}
	}
