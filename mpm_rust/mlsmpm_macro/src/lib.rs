extern crate proc_macro;
use proc_macro::TokenStream;
use quote::{quote, format_ident};
use syn::{parse::Parser, punctuated::Punctuated, Expr, ExprLit, Lit, Token};

fn multi_index(a: usize) -> (Vec<usize>, usize) {
    let mut res = vec![];
    for i in 0..=a {
        for dy in 0..=i {
            for dx in 0..=i {
                if dx + dy == i {
                    res.push(dx);
                    res.push(dy);
                }
            }
        }
    }
    let size = res.len() / 2;
    (res, size)
}

fn multi_index_factorial(a: usize) -> (Vec<f64>, usize) {
    fn factorial(num: usize) -> f64 {
        match num {
            0 | 1 => 1.,
            _ => factorial(num - 1) * num as f64,
        }
    }

    let mut res = vec![];
    for i in 0..=a {
        for dy in 0..=i {
            for dx in 0..=i {
                if dx + dy == i {
                    res.push(factorial(dx) * factorial(dy));
                    res.push((dx + dy) as f64);
                }
            }
        }
    }
    let size = res.len() / 2;
    (res, size)
}

#[proc_macro]
pub fn lsmps_poly(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (mut res, size) = multi_index(input);
    quote! {
        fn poly(r: Vector2<f64>) -> SVector<f64, #size> {
            let d = [#(#res,)*];
            let mut res = SVector::<f64, #size>::zeros();
            for i in 0..#size {
                res[i] = r[0].powi(d[i*2] as i32) * r[1].powi(d[i*2 + 1] as i32);
            }

            res
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_scale(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (mut res, size) = multi_index_factorial(input);

    quote! {
        fn scale(rs: f64) -> SMatrix::<f64, #size, #size> {
            let d = [#(#res,)*];
            let mut res = SVector::<f64, #size>::zeros();
            for i in 0..#size {
                res[i] = d[i*2] * rs.powf(-d[i*2 + 1]);
            }
            SMatrix::<f64, #size, #size>::from_diagonal(&res)
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_params(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (mut res, size) = multi_index(input);

    quote! {
        struct LsmpsParams {
            m: SMatrix<f64, #size, #size>,
            f_vel: SMatrix<f64, #size, 2>,
            f_stress: SMatrix<f64, #size, 3>,
            f_pressure: SVector<f64, #size>,
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_params_g2p(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (mut res, size) = multi_index(input);

    quote! {
        struct LsmpsParams {
            m: SMatrix<f64, #size, #size>,
            f_vel: SMatrix<f64, #size, 2>,
        }
    }
    .into()
}

fn parse_one_usize(input: TokenStream) -> usize {
    let input: Expr = syn::parse(input).unwrap();
    let p = if let Expr::Lit(expr_lit) = input {
        if let Lit::Int(litint) = &expr_lit.lit {
            litint.base10_parse::<usize>().unwrap()
        } else {
            panic!();
        }
    } else {
        panic!();
    };

    p
}

fn parse_two_usize(input: TokenStream) -> (usize, usize) {
    let parser = Punctuated::<Expr, Token![,]>::parse_separated_nonempty;
    let input = parser.parse(input).unwrap();
    let p = if let Expr::Lit(expr_lit) = input.first().unwrap() {
        if let Lit::Int(litint) = &expr_lit.lit {
            litint.base10_parse::<usize>().unwrap()
        } else {
            panic!();
        }
    } else {
        panic!();
    };
    let q = if let Expr::Lit(expr_lit) = input.last().unwrap() {
        if let Lit::Int(litint) = &expr_lit.lit {
            litint.base10_parse::<usize>().unwrap()
        } else {
            panic!();
        }
    } else {
        panic!();
    };

    (p, q)
}

#[proc_macro]
pub fn compact_lsmps_func(input: TokenStream) -> TokenStream {
    let (p, q) = parse_two_usize(input);

    let c_func_name = format_ident!("c_{}_{}", p, q);

    quote! {
        fn #c_func_name(bx: usize, by: usize) -> f64 {
            fn factorial(num: usize) -> f64 {
                match num {
                    0 | 1 => 1.,
                    _ => factorial(num - 1) * num as f64,
                }
            }
            
            let b = bx + by;
            if b == 0 {
                1.
            } else if b == #q {
                (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by)) * factorial(#p)
                    / factorial(#p + #q)
            } else {
                (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by))
                    * #q as f64
                    * factorial(#p + #q - b)
                    / factorial(#p + #q)
            }
        }
    }.into()
}

#[proc_macro]
pub fn compact_lsmps_scale(input: TokenStream) -> TokenStream {
    let (p, q) = parse_two_usize(input);

    let (p_res, p_size) = multi_index(p);
    let (q_res, q_size) = multi_index(q);

    let func_name = format_ident!("scale_{}_{}", p, q);
    let c_func_name = format_ident!("c_{}_{}", p, q);

    quote! {
        fn #func_name(rs: f64) -> SMatrix::<f64, #p_size, #p_size> {
            fn factorial(num: usize) -> f64 {
                match num {
                    0 | 1 => 1.,
                    _ => factorial(num - 1) * num as f64,
                }
            }
            
            let alpha = [#(#p_res,)*];
            let beta = [#(#q_res,)*];
            let mut res = SVector::<f64, #p_size>::zeros();
            for i in 0..#p_size {
                res[i] = {
                    let (ax, ay) = (alpha[i*2], alpha[i*2+1]);
                    let a = ax + ay;
                    let sum = (0..#q_size)
                        .map(|j| {
                            let (bx, by) = (beta[j*2], beta[j*2+1]);
                            let b = bx + by;
                            if b > #q_size || bx > ax || by > ay {
                                0.
                            } else {
                                #c_func_name(bx, by) / (factorial(ax - bx) * factorial(ay - by))
                            }
                        })
                        .sum::<f64>();
                    1. / (sum * rs.powi(a as i32))
                }
            }
            SMatrix::<f64, #p_size, #p_size>::from_diagonal(&res)
        }
    }
    .into()
}

#[proc_macro]
pub fn test_lsmps(input: TokenStream) -> TokenStream {
    quote! {
        fn lsmps_4th(settings: &Settings, space: &mut Space) {
        }
    }.into()
}

#[proc_macro]
pub fn lsmps_p2g_func(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let func_name = format_ident!("lsmps_{}th", p);
    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            mlsmpm_macro::lsmps_poly!(#p);
        
            let re = settings.cell_width() * settings.effect_radius as f64;
            let rs = settings.cell_width();
            mlsmpm_macro::lsmps_scale!(#p);
            let scale = scale(rs);
        
            mlsmpm_macro::lsmps_params!(#p);
        
            for p in space.particles.iter() {
                for node in NodeMutIterator::new(
                    settings,
                    &mut space.grid,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    let mass_contrib = node.weight * p.mass;
                    node.node.mass += mass_contrib;
                }
            }
        
            let mut nodes = HashMap::new();
        
            for p in space.particles.iter_mut() {
                let (stress, pressure) = {
                    let (density, volume) = calc_density_and_volume(
                        settings,
                        p,
                        &space.grid,
                        &space.period_bounds,
                        &space.period_bound_rect,
                    );
        
                    let pressure = match settings.pressure {
                        Some(pressure) => pressure(p, space.steps as f64 * settings.dt),
                        None => {
                            let pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                                * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                            if pressure < 0. {
                                0.
                            } else {
                                pressure
                            }
                        }
                    };
        
                    p.pressure = pressure;
        
                    let dudv = p.c;
                    let strain = dudv;
                    let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
        
                    (-pressure * Matrix2f::identity() + viscosity_term, pressure)
                };
        
                for node in NodeIterator::new(
                    settings,
                    &space.grid,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    if node.dist.norm() > re {
                        continue;
                    }
        
                    let params = {
                        let index = node.node.index;
                        if !nodes.contains_key(&index) {
                            let params = LsmpsParams {
                                m: SMatrix::zeros(),
                                f_vel: SMatrix::zeros(),
                                f_stress: SMatrix::zeros(),
                                f_pressure: SVector::zeros(),
                            };
                            nodes.insert(index, params);
                        }
        
                        nodes.get_mut(&node.node.index).unwrap()
                    };
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    // let weight = n.weight;
                    let weight = (1. - (node.dist / re).norm()).powi(2);
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();
                    params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose());
                    let stress = vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                    params.f_stress += weight * poly_r_ij.kronecker(&stress.transpose());
        
                    params.f_pressure += weight * poly_r_ij.kronecker(&Matrix1::new(pressure));
                }
            }
        
            space.grid.par_iter_mut().for_each(|node| {
                if !nodes.contains_key(&node.index) {
                    return;
                }
                let params = nodes.get(&node.index).unwrap();
                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    {
                        let res = scale * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
        
                        let pressure_res = scale * m_inverse * params.f_pressure;
                        node.force = -pressure_res.fixed_slice::<2, 1>(1, 0)
                            + settings.dynamic_viscosity
                                * settings.rho_0
                                * vector![res[(3, 0)] + res[(5, 0)], res[(3, 1)] + res[(5, 1)]];
                    }
        
                    // {
                    //     let res = scale * m_inverse * params.f_stress;
                    //     node.force[0] = res[(1, 0)] + res[(2, 1)];
                    //     node.force[1] = res[(1, 1)] + res[(2, 2)];
                    // }
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_g2p_func(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let func_name = format_ident!("lsmps_{}th", p);
    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            space.particles.par_iter_mut().for_each(|p| {
                p.v = Vector2f::zeros();
                p.c = Matrix2f::zeros();
        
                mlsmpm_macro::lsmps_poly!(#p);
        
                let re = settings.cell_width() * settings.effect_radius as f64;
                let rs = settings.cell_width();
        
                mlsmpm_macro::lsmps_scale!(#p);
                let scale = scale(rs);
        
                mlsmpm_macro::lsmps_params_g2p!(#p);
        
                let mut params = LsmpsParams {
                    m: SMatrix::zeros(),
                    f_vel: SMatrix::zeros(),
                };
        
                for n in NodeIterator::new(
                    settings,
                    &space.grid,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    if n.dist.norm() > re {
                        continue;
                    }
        
                    let r_ij = n.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    // let weight = n.weight;
                    let weight = (1. - (n.dist / re).norm()).powi(2);
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();
                    params.f_vel += weight * poly_r_ij.kronecker(&n.node.v_star.transpose());
                }
        
                if let Some(m_inverted) = params.m.try_inverse() {
                    let res = scale * m_inverted * params.f_vel;
                    p.v = res.row(0).transpose();
        
                    if settings.vx_zero {
                        p.v.x = 0.;
                    }
                    p.x += res.row(0).transpose() * settings.dt;
                    p.c = res.fixed_slice::<2, 2>(1, 0).transpose().into();
                    //p.c = Matrix2::new(0., 0., -0.1, 0.);
                }
            });
        }
    }
    .into()
}