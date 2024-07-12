pub mod calculos{
		// conjunto de funcionalidades 
		pub mod linear_func{
				#[derive(PartialEq, Clone)]
				pub struct LinearFunc{
					pub a: f64,
					pub b: f64,
				}

				impl LinearFunc{
					pub fn calcule(&self, x: &f64) -> f64{
						self.a * x + self.b
					}
					pub fn cruzamento_com_eixo_x(&self) -> f64{
						-self.b / self.a
					}
					pub fn cruzamento_com_eixo_y(&self) -> f64{
						todo!()				
					}
					pub fn novo_de_ponto_e_derivada(point: (&f64, &f64), deriv: &f64) -> LinearFunc{
						let a = *deriv;	
						let b = -(a * point.0) + point.1;
						LinearFunc{a , b}
					}
				}
				
			}
		
		pub mod diff{
use nalgebra::Vector;

use crate::{ferramentas::graficos::linspace, Vetor};

			pub mod escalar{
				pub fn derivada_numerica(delta: f64, mut f: impl Fn(&f64) -> f64 + Clone) -> impl (Fn(&f64) -> f64) + Clone{
					move |x|{
						let x_plus_delta = x + delta/2.0;
						let x_minus_delta  = x - delta/2.0;
						let y_of_plus_delta = f(&x_plus_delta);
						let y_of_minus_delta = f(&x_minus_delta);
						let df = (y_of_plus_delta - y_of_minus_delta)/delta; 
						df
					}
				}
			}
		
		
			pub mod vetorial{
				use nalgebra::{Const, SMatrix};
				use crate::Vetor;

					
					/// Essa função recebe uma função e um parametro delta e retorna uma outra função a qual retorna uma estimativa da derivada no ponto. Essa estimativa é calculada numericamente.
					pub fn derivada_numerica
						(delta: f64, mut f: impl Fn(&Vetor) -> f64) -> impl Fn(&Vetor) -> Vetor
					{ 
						move |x| {
							let x_dim = x.len();
							let mut ans = Vetor::zeros(x_dim);
							for ix in 0..x_dim {
								let delta_vec = {
									let mut delta_vec = Vetor::zeros(x_dim);
									delta_vec[ix] = delta;
									delta_vec
								};
								let x_mais_delta    = x + delta_vec.clone();
								let x_menos_delta   = x - delta_vec.clone();

								let y_de_mais_delta = f(&x_mais_delta); 
								let y_de_menos_delta = f(&x_menos_delta); 

								let delta_y = y_de_mais_delta - y_de_menos_delta; 
								let df = delta_y / (delta * 2.0);
								ans[ix] = df;
							}
							ans
						}
					}
			}

			pub mod matriz{
				use nalgebra::SMatrix;
			use crate::{Matriz, Vetor};
				

				pub fn derivada_numerica
					(   delta: f64,
						mut f: impl Fn(&Vetor) -> Vetor
					) 
					-> impl Fn(&Vetor) -> Matriz
				{
					move |x| {
						let in_dim  = x.len();
						let out_dim = (f(x)).len(); 

						let mut ans = Matriz::zeros(in_dim, out_dim);
						for ix in 0..in_dim{
							let delta_vec = {
								let mut delta_vec = Vetor::zeros(in_dim);
								delta_vec[ix] = delta  ;
								delta_vec
							};
							let x_mais_delta    = x + delta_vec.clone();
							let x_menos_delta   = x - delta_vec.clone();

							let y_de_mais_delta = f(&x_mais_delta); 
							let y_de_menos_delta = f(&x_menos_delta); 

							let delta_y = y_de_mais_delta - y_de_menos_delta; 
							let df = delta_y / (delta * 2.0);
							for iy in 0..out_dim{
								ans[iy * in_dim + ix] = df[iy]
							}
						}
						ans
					}		
				}

			}

			#[test]
			fn teste_derivada_matriz(){

				let f= |v: &Vetor| -> Vetor{
					let x1 = v[0].clone();
					let x2 = v[1].clone();
					let y1 = x1.powi(2) * 3.0 + -2.0 * x1 + 10.0*x2 + 2.0*x1*x2;
					let y2 = x1 + -2.0 * x1 + 10.0*x2 + 2.0*x1*x2;
					let y3 = x2*1.2;

					Vetor::from_row_slice(&[y1, y2, y3])
				};

				let delta = 10e-4;
				let df = matriz::derivada_numerica(delta, &f); 

				let exemplo_x = Vetor::zeros(2);
				let exemplo_y = df(&exemplo_x);
				println!("{:?}", &exemplo_y);

			}
		}
	}

	pub mod graficos{
	use std::ops::Range;
		use proptest::proptest;
		pub fn linspace(range: Range<f64>, num_pontos: u32) -> impl ExactSizeIterator<Item=f64> + Clone{
			(0..num_pontos)
			.into_iter()
			.map({
				let passo = (range.end - range.start) / ((num_pontos - 1) as f64);
				move |n| {
					return passo * (n as f64) + range.start;
				}
			})
		}

		proptest!{ 
			#[test]
			fn teste_linspace(start:f64, end:f64, num_pontos in 0..100_u32) {
				let space = linspace(start..end, num_pontos).collect::<Vec<_>>();
				println!("([{start}..{end};{num_pontos}]) -> {space:?}");
			}
		}
	}

	pub mod exemplos{
		pub mod escalar{
		use crate::PassoEscalar;
			pub fn obtenha_f_e_passo_inicial() -> (impl Fn(&f64) -> f64, PassoEscalar){
				let f = |x:&f64| x.powf(2.0) + 3.0*x + 4.0 ;

				let x_inicial = 15.0;
				let y_inicial = f(&x_inicial);
				let passo_inicial = PassoEscalar{x: x_inicial, y: y_inicial};

				(f, passo_inicial)
			}
		}
		pub mod vetorial{ 
			use crate::{Passo, PassoVet, Vetor};
			pub fn obtenha_f_e_passo_inicial() -> (impl Fn(&Vetor) -> f64 + Clone, PassoVet){
				let f = |x: &Vetor| -> f64{
					let x_0 = x[0]; 
					let x_1 = x[1];
					x_0.powf(2.0) + x_0*3.0 + 4.0 + x_1.powi(2)*2.0 
				};

				let x_inicial = Vetor::from_row_slice(&[10.0, 15.0]);
				let y_inicial = f(&x_inicial);
				let passo_inicial = Passo{
					x: x_inicial, 
					y: y_inicial
				};

				(f, passo_inicial)
			}
		}
	}
