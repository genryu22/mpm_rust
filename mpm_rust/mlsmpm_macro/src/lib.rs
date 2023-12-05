extern crate proc_macro;
use proc_macro::TokenStream;
use quote::{quote, format_ident};
use syn::{parse::Parser, punctuated::Punctuated, Expr, Lit, Token};

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

fn multi_index_vec(a: usize) -> (Vec<usize>, Vec<usize>) {
    let mut res_x = vec![];
    let mut res_y = vec![];
    for i in 0..=a {
        for dy in 0..=i {
            for dx in 0..=i {
                if dx + dy == i {
                    res_x.push(dx);
                    res_y.push(dy);
                }
            }
        }
    }
    (res_x, res_y)
}

fn multi_index_1d(a: usize) -> (Vec<usize>, usize) {
    let mut res = vec![];
    for i in 0..=a {
         res.push(i);
    }
    let size = res.len();
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

fn multi_index_factorial_1d(a: usize) -> (Vec<f64>, usize) {
    fn factorial(num: usize) -> f64 {
        match num {
            0 | 1 => 1.,
            _ => factorial(num - 1) * num as f64,
        }
    }

    let mut res = vec![];
    for i in 0..=a {
        res.push(factorial(i));
        res.push(i as f64);
    }
    let size = res.len() / 2;
    (res, size)
}

#[proc_macro]
pub fn lsmps_poly(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);

    let (dx, dy) = multi_index_vec(input);
    let size = dx.len();

    println!("p={}, n of alpha(|alpha|=p) = {}", input, size);

    quote! {
        fn poly(r: Vector2<f64>) -> SVector<f64, #size> {
            vector![#(r[0].powi(#dx as i32) * r[1].powi(#dy as i32)),*]
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_poly_1d(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (res, size) = multi_index_1d(input);

    quote! {
        fn poly(r: f64) -> SVector<f64, #size> {
            let d = [#(#res),*];

            let mut res = SVector::<f64, #size>::zeros();
            for i in 0..#size {
                res[i] = r.powi(d[i] as i32);
            }

            res
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_scale(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (res, size) = multi_index_factorial(input);

    quote! {
        fn scale(rs: f64) -> SMatrix::<f64, #size, #size> {
            let d = [#(#res),*];
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
pub fn lsmps_scale_1d(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (res, size) = multi_index_factorial_1d(input);

    quote! {
        fn scale(rs: f64) -> SMatrix::<f64, #size, #size> {
            let d = [#(#res),*];
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
    let (_, size) = multi_index(input);

    quote! {
        struct LsmpsParams {
            m: SMatrix<f64, #size, #size>,
            f_vel: SMatrix<f64, #size, 2>,
            f_stress: SMatrix<f64, #size, 3>,
            f_pressure: SVector<f64, #size>,
            count: usize,
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_params_vel(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (_, size) = multi_index(input);

    quote! {
        struct LsmpsParams {
            m: SMatrix<f64, #size, #size>,
            f_vel: SMatrix<f64, #size, 2>,
            count: usize,
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_params_empty(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (_, size) = multi_index(input);

    quote! {
        (SMatrix::<f64, #size, #size>::zeros(), SMatrix::<f64, #size, 2>::zeros())
    }
    .into()
}

#[proc_macro]
pub fn lsmps_params_1d(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (_, size) = multi_index_1d(input);

    quote! {
        struct LsmpsParams {
            m: SMatrix<f64, #size, #size>,
            f_vel: SVector<f64, #size>,
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_params_g2p(input: TokenStream) -> TokenStream {
    let input: usize = parse_one_usize(input);
    let (_, size) = multi_index(input);

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
    let p = match input {
        Expr::Lit(expr_lit) => {
            if let Lit::Int(litint) = &expr_lit.lit {
                litint.base10_parse::<usize>().unwrap()
            } else {
                panic!("int でない。");
            }
        },
        _ => {
            panic!("expr でない。");
        }
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

            // (-1. as f64).powi(b as i32) / (factorial(bx) * factorial(by)) * factorial(#q) / factorial(#q - b) * factorial(#p + #q - b) / factorial(#p + #q)
        }
    }.into()
}

#[proc_macro]
pub fn compact_scale(input: TokenStream) -> TokenStream {
    let (p, q) = parse_two_usize(input);

    let (ax, ay) = multi_index_vec(p);
    let (bx, by) = multi_index_vec(q);
    let (p_size, q_size) = (ax.len(), bx.len());

    let scale_p_q = format_ident!("scale_{}_{}", p, q);

    quote! {
        fn #scale_p_q(rs: f64) -> SMatrix::<f64, #p_size, #p_size> {
            fn factorial(num: usize) -> f64 {
                match num {
                    0 | 1 => 1.,
                    _ => factorial(num - 1) * num as f64,
                }
            }

            fn c(bx: usize, by: usize) -> f64 {
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
    
                // (-1. as f64).powi(b as i32) / (factorial(bx) * factorial(by)) * factorial(#q) / factorial(#q - b) * factorial(#p + #q - b) / factorial(#p + #q)
            }

            let (ax, ay) = ([#(#ax),*], [#(#ay),*]);
            let (bx, by) = ([#(#bx),*], [#(#by),*]);
            let mut diag = SVector::<f64, #p_size>::zeros();
            for i in 0..#p_size {
                diag[i] = 1. / (0..#q_size).map(|j| {
                    if bx[j] > ax[i] || by[j] > ay[i] {
                        0.
                    } else {
                        c(bx[j], by[j]) / (factorial(ax[i] - bx[j]) * factorial(ay[i] - by[j]))
                    }
                }).sum::<f64>() / rs.powi((ax[i] + ay[i]) as i32);
            }

            SMatrix::<f64, #p_size, #p_size>::from_diagonal(&diag)
        }
    }.into()
}

#[proc_macro]
pub fn lsmps_p2g_func(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let func_name = format_ident!("lsmps_{}", p);
    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_scale!(#p);
            mlsmpm_macro::lsmps_params!(#p);
        
            let rs = settings.cell_width();
            let scale = scale(rs);

            parallel!(settings, space.particles, |p| {
                for node in
                    NodeIndexIterator::new(settings, p, &space.period_bounds, &space.period_bound_rect)
                {
                    let mass_contrib = node.weight * p.mass;
                    {
                        let mut node = space.grid[node.index].lock().unwrap();
                        node.mass += mass_contrib;
                    }
                }
            });

            let nodes = (0..space.grid.len()).into_par_iter().map(|_| {
                std::sync::Mutex::new(LsmpsParams {
                    m: SMatrix::zeros(),
                    f_vel: SMatrix::zeros(),
                    f_stress: SMatrix::zeros(),
                    f_pressure: SVector::zeros(),
                    count: 0,
            })}).collect::<Vec<_>>();

            parallel!(settings, space.particles, |p| {
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
        
                    (match settings.pressure_grad {
                        Some(_) => viscosity_term,
                        None => (-pressure * Matrix2f::identity() + viscosity_term)
                    }, pressure)
                };
        
                for node in NodeIndexIterator::new(
                    settings,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();
                    params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose());
                    let stress = vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                    params.f_stress += weight * poly_r_ij.kronecker(&stress.transpose());
                    params.count += 1;
                }
            });

            parallel!(settings, space.grid, |node| {
                let mut node = node.lock().unwrap();
                let index = node.index.0 + node.index.1 * (settings.grid_width + 1);
                let params = nodes[index].lock().unwrap();

                if params.count == 0 {
                    return;
                }

                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    {
                        let res = scale * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                    }

                    {
                        let res = scale * m_inverse * params.f_stress;
                        node.force = match settings.pressure_grad {
                            Some(pressure_grad) => {
                                let x = node.index.0 as f64 * settings.cell_width();
                                let y = node.index.1 as f64 * settings.cell_width();
    
                                -settings.pressure_grad.unwrap()(x, y, space.steps as f64 * settings.dt)
                            },
                            None => Vector2f::zeros(),
                        } + vector![res[(1, 0)] + res[(2, 1)], res[(1, 1)] + res[(2, 2)]];
                    }
                } else {
                    let sing_values = params.m.singular_values();
                    println!("モーメント行列の逆行列が求まりませんでした。 max={} min={} cond={} (x, y) = ({} {}), steps={}",
                        sing_values.max(),
                        sing_values.min(),
                        sing_values.max() / sing_values.min(),
                        node.index.0 as f64 * settings.cell_width(),
                        node.index.1 as f64 * settings.cell_width(),
                        space.steps
                    );
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_p2g_laplacian_func(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let func_name = format_ident!("lsmps_laplacian_{}", p);
    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_scale!(#p);
            mlsmpm_macro::lsmps_params!(#p);
        
            let rs = settings.cell_width();
            let scale = scale(rs);

            let nodes = (0..space.grid.len()).into_par_iter().map(|_| {
                std::sync::Mutex::new(LsmpsParams {
                    m: SMatrix::zeros(),
                    f_vel: SMatrix::zeros(),
                    f_stress: SMatrix::zeros(),
                    f_pressure: SVector::zeros(),
                    count: 0,
            })}).collect::<Vec<_>>();

            parallel!(settings, space.particles, |p| {
                for node in NodeIndexIterator::new(
                    settings,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();
                    params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose());
                    params.count += 1;
                }
            });

            parallel!(settings, space.grid, |node| {
                let mut node = node.lock().unwrap();
                let index = node.index.0 + node.index.1 * (settings.grid_width + 1);
                let params = nodes[index].lock().unwrap();

                if params.count == 0 {
                    return;
                }

                if let Some(m_inverse) = (params.m).try_inverse() {
                    {
                        let res = scale * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();

                        {
                            let x = node.index.0 * settings.cell_width();
                            let y = node.index.1 * settings.cell_width();
                            node.force = -settings.pressure_grad.unwrap()(x, y, space.steps as f64 * settings.dt) + settings.dynamic_viscosity
                                * settings.rho_0
                                * vector![res[(3, 0)] + res[(5, 0)], res[(3, 1)] + res[(5, 1)]];
                        }
                    }
                } else {
                    let sing_values = params.m.singular_values();
                    println!("モーメント行列の逆行列が求まりませんでした。 max={} min={} cond={} (x, y) = ({} {}), steps={}",
                        sing_values.max(),
                        sing_values.min(),
                        sing_values.max() / sing_values.min(),
                        node.index.0 as f64 * settings.cell_width(),
                        node.index.1 as f64 * settings.cell_width(),
                        space.steps
                    );
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn lsmps_g2p_func(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let func_name = format_ident!("lsmps_{}", p);
    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            parallel!(settings, space.particles, |p| {
                p.v = Vector2f::zeros();
                p.c = Matrix2f::zeros();
        
                mlsmpm_macro::lsmps_poly!(#p);
                mlsmpm_macro::lsmps_scale!(#p);
                mlsmpm_macro::lsmps_params_g2p!(#p);
        
                // let re = settings.cell_width() * settings.effect_radius as f64;
                let rs = settings.cell_width();
                let scale = scale(rs);
        
                let mut params = LsmpsParams {
                    m: SMatrix::zeros(),
                    f_vel: SMatrix::zeros(),
                };
        
                for n in NodeIndexIterator::new(
                    settings,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    let r_ij = n.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = n.weight;
                    // let weight = (1. - (n.dist / re).norm()).powi(2);
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();

                    {
                        let node = space.grid[n.index].lock().unwrap();
                        params.f_vel += weight * poly_r_ij.kronecker(&node.v_star.transpose());
                    }
                }
        
                if let Some(m_inverted) = params.m.try_inverse() {
                    let res = scale * m_inverted * params.f_vel;
                    p.v = res.row(0).transpose();
        
                    if settings.vx_zero {
                        p.v.x = 0.;
                    }
                    if !settings.calc_convection_term {
                        p.x += p.v * settings.dt;
                    }
                    p.c = res.fixed_view::<2, 2>(1, 0).transpose().into();
                    p.v_lsmps = res.rows(0, res.shape().0).transpose();
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn compact_p2g_func(input: TokenStream) -> TokenStream {
    let (p, q) = parse_two_usize(input);
    let func_name = format_ident!("compact_{}_{}", p, q);
    let scale_p_0 = format_ident!("scale_{}_0", p);
    let scale_p_q = format_ident!("scale_{}_{}", p, q);
    let c_p_q = format_ident!("c_{}_{}", p, q);
    let (dx, dy) = multi_index_vec(q);
    let i = 0..(dx.len());

    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            parallel!(settings, space.particles, |p| {
                for node in
                    NodeIndexIterator::new(settings, p, &space.period_bounds, &space.period_bound_rect)
                {
                    let mass_contrib = node.weight * p.mass;
                    {
                        let mut node = space.grid[node.index].lock().unwrap();
                        node.mass += mass_contrib;
                    }
                }
            });
        
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_params!(#p);
            mlsmpm_macro::compact_lsmps_func!(#p, 0);
            mlsmpm_macro::compact_lsmps_func!(#p, #q);
            mlsmpm_macro::compact_scale!(#p, 0);
            mlsmpm_macro::compact_scale!(#p, #q);
        
            let rs = settings.cell_width();
            let scale_stress = #scale_p_0(rs);
            let scale_vel = #scale_p_q(rs);

            let nodes = (0..space.grid.len()).into_par_iter().map(|_| {
                std::sync::Mutex::new(LsmpsParams {
                    m: SMatrix::zeros(),
                    f_vel: SMatrix::zeros(),
                    f_stress: SMatrix::zeros(),
                    f_pressure: SVector::zeros(),
                    count: 0,
            })}).collect::<Vec<_>>();
        
            parallel!(settings, space.particles, |p| {
                let stress = {
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
        
                    match settings.pressure_grad {
                        Some(_) => viscosity_term,
                        None => (-pressure * Matrix2f::identity() + viscosity_term)
                    }
                };
        
                for node in NodeIndexIterator::new(
                    settings,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();

                    #({
                        let (i, bx, by) = (#i, #dx, #dy);
                        // きれいに書く
                        if p.v_lsmps.shape().1 > 1 {
                            params.f_vel += weight * poly_r_ij.kronecker(&p.v_lsmps.column(i).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                        } else {
                            // step = 
                            if i == 0 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            } else if i == 1 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.c.column(0).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            } else if i == 2 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.c.column(1).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            }
                        }
                    })*
        
                    let stress = vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                    params.f_stress += weight * poly_r_ij.kronecker(&stress.transpose());
                    params.count += 1;
                }
            });

            parallel!(settings, space.grid, |node| {
                let mut node = node.lock().unwrap();
                let index = node.index.0 + node.index.1 * (settings.grid_width + 1);
                let params = nodes[index].lock().unwrap();

                if params.count == 0 {
                    return;
                }
            
                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    {
                        let res = scale_vel * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                        if settings.calc_convection_term {
                            node.force = -settings.rho_0 * Vector2::<f64>::new(
                                res[(0, 0)] * res[(1, 0)] + res[(0, 1)] * res[(2, 0)],
                                res[(0, 0)] * res[(1, 1)] + res[(0, 1)] * res[(2, 1)],
                            );
                        }
                    }

                    {
                        let res = scale_stress * m_inverse * params.f_stress;
                        let force = node.force;
                        node.force = force + match settings.pressure_grad {
                            Some(pressure_grad) => {
                                let x = node.index.0 as f64 * settings.cell_width();
                                let y = node.index.1 as f64 * settings.cell_width();
    
                                -settings.pressure_grad.unwrap()(x, y, space.steps as f64 * settings.dt)
                            },
                            None => Vector2f::zeros(),
                        } + vector![res[(1, 0)] + res[(2, 1)], res[(1, 1)] + res[(2, 2)]];
                    }
                } else {
                    let sing_values = params.m.singular_values();
                    println!("モーメント行列の逆行列が求まりませんでした。 max={} min={} cond={}", sing_values.max(), sing_values.min(), sing_values.max() / sing_values.min());
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn compact_p2g_laplacian_func(input: TokenStream) -> TokenStream {
    let (p, q) = parse_two_usize(input);
    let func_name = format_ident!("compact_laplacian_{}_{}", p, q);
    let scale_p_q = format_ident!("scale_{}_{}", p, q);
    let c_p_q = format_ident!("c_{}_{}", p, q);
    let (dx, dy) = multi_index_vec(q);
    let i = 0..(dx.len());

    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_params_vel!(#p);
            mlsmpm_macro::compact_lsmps_func!(#p, #q);
            mlsmpm_macro::compact_scale!(#p, #q);
        
            let rs = settings.cell_width();
            let scale_vel = #scale_p_q(rs);

            let nodes = (0..space.grid.len()).into_par_iter().map(|_| {
                std::sync::Mutex::new(LsmpsParams {
                    m: SMatrix::zeros(),
                    f_vel: SMatrix::zeros(),
                    count: 0,
            })}).collect::<Vec<_>>();
        
            parallel!(settings, space.particles, |p| {
                for node in NodeIndexIterator::new(
                    settings,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();

                    #({
                        let (i, bx, by) = (#i, #dx, #dy);
                        // きれいに書く
                        if p.v_lsmps.shape().1 > 1 {
                            params.f_vel += weight * poly_r_ij.kronecker(&p.v_lsmps.column(i).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                        } else {
                            // step = 0
                            if i == 0 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            } else if i == 1 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.c.column(0).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            } else if i == 2 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.c.column(1).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            }
                        }
                    })*

                    params.count += 1;
                }
            });

            parallel!(settings, space.grid, |node| {
                let mut node = node.lock().unwrap();
                let index = node.index.0 + node.index.1 * (settings.grid_width + 1);
                let params = nodes[index].lock().unwrap();

                if params.count == 0 {
                    return;
                }
            
                if let Some(m_inverse) = params.m.try_inverse() {
                    {
                        let res = scale_vel * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                        if settings.calc_convection_term {
                            node.force = -settings.rho_0 * Vector2::<f64>::new(
                                res[(0, 0)] * res[(1, 0)] + res[(0, 1)] * res[(2, 0)],
                                res[(0, 0)] * res[(1, 1)] + res[(0, 1)] * res[(2, 1)],
                            );
                        }

                        {
                            let x = node.index.0 as f64 * settings.cell_width();
                            let y = node.index.1 as f64 * settings.cell_width();
                            let force = node.force;
                            node.force = force + -settings.pressure_grad.unwrap()(x, y, space.steps as f64 * settings.dt) + settings.dynamic_viscosity
                                * settings.rho_0
                                * vector![res[(3, 0)] + res[(5, 0)], res[(3, 1)] + res[(5, 1)]];
                        }
                    }
                } else {
                    let sing_values = params.m.singular_values();
                    println!("モーメント行列の逆行列が求まりませんでした。 max={} min={} cond={}", sing_values.max(), sing_values.min(), sing_values.max() / sing_values.min());
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn compact_v_p2g_func(input: TokenStream) -> TokenStream {
    let (p, q) = parse_two_usize(input);
    let func_name = format_ident!("compact_v_{}_{}", p, q);
    let scale_p_q = format_ident!("scale_{}_{}", p, q);
    let c_p_q = format_ident!("c_{}_{}", p, q);
    let (dx, dy) = multi_index_vec(q);
    let i = 0..(dx.len());

    quote! {
        fn #func_name(settings: &Settings, space: &mut Space) {
            mlsmpm(settings, space);
        
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_params!(#p);
            mlsmpm_macro::compact_lsmps_func!(#p, #q);
            mlsmpm_macro::compact_scale!(#p, #q);
        
            let rs = settings.cell_width();
            let scale_vel = #scale_p_q(rs);

            let nodes = (0..space.grid.len()).into_par_iter().map(|_| {
                std::sync::Mutex::new(LsmpsParams {
                    m: SMatrix::zeros(),
                    f_vel: SMatrix::zeros(),
                    f_stress: SMatrix::zeros(),
                    f_pressure: SVector::zeros(),
                    count: 0,
            })}).collect::<Vec<_>>();
        
            parallel!(settings, space.particles, |p| {
                for node in NodeIndexIterator::new(
                    settings,
                    p,
                    &space.period_bounds,
                    &space.period_bound_rect,
                ) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();

                    #({
                        let (i, bx, by) = (#i, #dx, #dy);
                        if p.v_lsmps.shape().1 > 1 {
                            params.f_vel += weight * poly_r_ij.kronecker(&p.v_lsmps.column(i).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                        } else {
                            // step = 0
                            if i == 0 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            } else if i == 1 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.c.column(0).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            } else if i == 2 {
                                params.f_vel += weight * poly_r_ij.kronecker(&p.c.column(1).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                            }
                        }
                    })*

                    params.count += 1;
                }
            });

            parallel!(settings, space.grid, |node| {
                let mut node = node.lock().unwrap();
                let index = node.index.0 + node.index.1 * (settings.grid_width + 1);
                let params = nodes[index].lock().unwrap();

                if params.count == 0 {
                    return;
                }

                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    {
                        let res = scale_vel * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                        if settings.calc_convection_term {
                            let force = node.force;
                            node.force = force - settings.rho_0 * Vector2::<f64>::new(
                                res[(0, 0)] * res[(1, 0)] + res[(0, 1)] * res[(2, 0)],
                                res[(0, 0)] * res[(1, 1)] + res[(0, 1)] * res[(2, 1)],
                            );
                        }
                    }
                } else {
                    let sing_values = params.m.singular_values();
                    println!("モーメント行列の逆行列が求まりませんでした。 {}", sing_values.max() / sing_values.min());
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn test_dissipation_lsmps_func(input: TokenStream) -> TokenStream {
    let (p, e) = parse_two_usize(input);
    let p2g_func_name = format_ident!("lsmps_p2g_{}", p);
    let g2p_func_name = format_ident!("lsmps_g2p_{}", p);
    quote! {
        fn #p2g_func_name(particles: &Vec<Particle>, grid: &mut Vec<Node>, cell_width: f64, grid_size: usize) {
            mlsmpm_macro::lsmps_poly_1d!(#p);
            mlsmpm_macro::lsmps_scale_1d!(#p);
            mlsmpm_macro::lsmps_params_1d!(#p);
        
            let rs = cell_width;
            let scale = scale(rs);
        
            let params = particles
                .par_iter()
                .flat_map(|p| {
                    let e = #e as i32;
                    (-e..=e)
                        .map(|gx| {
                            let base_node = (p.x / cell_width).floor() as i32;
                            (
                                (base_node + gx) as f64 * cell_width,
                                (base_node + gx).rem_euclid(grid_size as i32 + 1) as usize,
                            )
                        })
                        .filter(|(n_pos, _)| (n_pos - p.x).abs() <= e as f64 * cell_width)
                        .map(|(n_pos, n_index)| {
                            let dist = n_pos - p.x;
                            let r_ij = -dist / rs;
                            let poly_r_ij = poly(r_ij);
                            let weight = (1. - (dist / (e as f64 * cell_width)).abs()).powi(2);
        
                            (
                                n_index,
                                LsmpsParams {
                                    m: weight * poly_r_ij * poly_r_ij.transpose(),
                                    f_vel: weight * poly_r_ij.kronecker(&Vector1::new(p.v)),
                                },
                            )
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
        
            let mut params_sums = {
                let mut params_sums = vec![];
                for _ in 0..grid.len() {
                    params_sums.push(LsmpsParams {
                        m: SMatrix::zeros(),
                        f_vel: SVector::zeros(),
                    });
                }
                params_sums
            };
            for (n_index, params) in params {
                params_sums[n_index].m += params.m;
                params_sums[n_index].f_vel += params.f_vel;
            }
        
            for (index, vel) in params_sums
                .par_iter_mut()
                .map(|params| {
                    if let Some(m_inverse) = params.m.try_inverse() {
                        let res = scale * m_inverse * params.f_vel;
                        res.row(0).transpose().x
                    } else {
                        0.
                    }
                })
                .enumerate()
                .collect::<Vec<_>>()
            {
                grid[index].v = vel;
            }
        }

        fn #g2p_func_name(particles: &mut Vec<Particle>, grid: &Vec<Node>, cell_width: f64, grid_size: usize) {
            mlsmpm_macro::lsmps_poly_1d!(#p);
            mlsmpm_macro::lsmps_scale_1d!(#p);
            mlsmpm_macro::lsmps_params_1d!(#p);
        
            let rs = cell_width;
            let scale = scale(rs);
        
            impl std::iter::Sum for LsmpsParams {
                fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                    iter.reduce(|mut acc, e| {
                        acc.m += e.m;
                        acc.f_vel += e.f_vel;
                        acc
                    })
                    .unwrap_or(LsmpsParams {
                        m: SMatrix::zeros(),
                        f_vel: SVector::zeros(),
                    })
                }
            }
        
            particles.iter_mut().for_each(|p| {
                let e = #e as i32;
                let params = (-e..=e)
                    .map(|gx| {
                        let base_node = (p.x / cell_width).floor() as i32;
                        (
                            (base_node + gx) as f64 * cell_width,
                            (base_node + gx).rem_euclid(grid_size as i32 + 1) as usize,
                        )
                    })
                    .filter(|(n_pos, _)| (n_pos - p.x).abs() <= e as f64 * cell_width)
                    .map(|(n_pos, n_index)| {
                        let dist = n_pos - p.x;
                        let r_ij = dist / rs;
                        let poly_r_ij = poly(r_ij);
                        let weight = (1. - (dist / (e as f64 * cell_width)).abs()).powi(2);
        
                        LsmpsParams {
                            m: weight * poly_r_ij * poly_r_ij.transpose(),
                            f_vel: weight * poly_r_ij.kronecker(&Vector1::new(grid[n_index].v)),
                        }
                    })
                    .sum::<LsmpsParams>();
        
                if let Some(m_inverted) = params.m.try_inverse() {
                    let res = scale * m_inverted * params.f_vel;
                    p.v = res.row(0).x;
                    p.c = res.row(1).x;
                }
            });
        }
    }
    .into()
}

#[proc_macro]
pub fn scheme_div_force(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let func_name = format_ident!("scheme_div_force_{}", p);

    quote! {
        fn #func_name(settings: &Settings) -> Vec<(f64, f64, SVector<f64, 2>)> {
            let (mut particles, mut grid, periodic_boundary_rect) = new_for_taylor_green(settings);

            mlsmpm::g2p_lsmps_2nd(
                &mut particles,
                &grid,
                settings,
                periodic_boundary_rect.clone(),
            );

            grid.par_iter_mut().for_each(|node| node.reset());

            let period_bounds = vec![];
            let periodic_boundary_rect = Some(periodic_boundary_rect);
        
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_scale!(#p);
            mlsmpm_macro::lsmps_params!(#p);
        
            let rs = settings.cell_width();
            let scale = scale(rs);
        
            let nodes = {
                let mut nodes = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    nodes.push(std::sync::Mutex::new(LsmpsParams {
                        m: SMatrix::zeros(),
                        f_vel: SMatrix::zeros(),
                        f_stress: SMatrix::zeros(),
                        f_pressure: SVector::zeros(),
                        count: 0,
                    }));
                }
                nodes
            };
        
            parallel!(settings, particles, |p| {
                let viscosity_term = {
                    let dudv = p.c();
                    let strain = dudv;
                    let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
        
                    viscosity_term
                };
        
                for node in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();
        
                    let stress = vector![
                        viscosity_term[(0, 0)],
                        viscosity_term[(0, 1)],
                        viscosity_term[(1, 1)]
                    ];
                    params.f_stress += weight * poly_r_ij.kronecker(&stress.transpose());
                    params.count += 1;
                }
            });
        
            let result = {
                let mut result = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    result.push(Mutex::new([0.; 2]));
                }
                result
            };

            nodes.par_iter().enumerate().for_each(|(index, params)| {
                let params = params.lock().unwrap();

                if params.count == 0 {
                    return;
                }

                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    let res = scale * m_inverse * params.f_stress;
        
                    let node_force = SVector::<f64, 2>::new(res[(1, 0)] + res[(2, 1)], res[(1, 1)] + res[(2, 2)]);
                    let mut result = result[index].lock().unwrap();
                    result[0] += node_force.x;
                    result[1] += node_force.y;
                }
            });
        
            grid.par_iter()
                .enumerate()
                .map(|(index, _)| {
                    (
                        (index % (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        (index / (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        {
                            let result = result[index].lock().unwrap();
                            Vector2::<f64>::new(result[0], result[1])
                        }
                    )
                })
                .filter(|(x, y, _)| 4. <= *x && *x < 6. && 4. <= *y && *y < 6.)
                .collect::<Vec<_>>()
        }
    }.into()
}
#[proc_macro]
pub fn scheme_div_force_1st(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let func_name = format_ident!("scheme_div_force_{}_1st", p);

    quote! {
        fn #func_name(settings: &Settings) -> Vec<(f64, f64, SVector<f64, 2>)> {
            let (mut particles, mut grid, periodic_boundary_rect) = new_for_taylor_green(settings);

            mlsmpm::g2p_lsmps_1st(
                &mut particles,
                &grid,
                settings,
                periodic_boundary_rect.clone(),
            );

            grid.par_iter_mut().for_each(|node| node.reset());

            let period_bounds = vec![];
            let periodic_boundary_rect = Some(periodic_boundary_rect);
        
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_scale!(#p);
            mlsmpm_macro::lsmps_params!(#p);
        
            let rs = settings.cell_width();
            let scale = scale(rs);
        
            let nodes = {
                let mut nodes = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    nodes.push(std::sync::Mutex::new(LsmpsParams {
                        m: SMatrix::zeros(),
                        f_vel: SMatrix::zeros(),
                        f_stress: SMatrix::zeros(),
                        f_pressure: SVector::zeros(),
                        count: 0,
                    }));
                }
                nodes
            };
        
            parallel!(settings, particles, |p| {
                let viscosity_term = {
                    let dudv = p.c();
                    let strain = dudv;
                    let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
        
                    viscosity_term
                };
        
                for node in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();
        
                    let stress = vector![
                        viscosity_term[(0, 0)],
                        viscosity_term[(0, 1)],
                        viscosity_term[(1, 1)]
                    ];
                    params.f_stress += weight * poly_r_ij.kronecker(&stress.transpose());
                    params.count += 1;
                }
            });
        
            let result = {
                let mut result = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    result.push(Mutex::new([0.; 2]));
                }
                result
            };

            nodes.par_iter().enumerate().for_each(|(index, params)| {
                let params = params.lock().unwrap();

                if params.count == 0 {
                    return;
                }

                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    let res = scale * m_inverse * params.f_stress;
        
                    let node_force = SVector::<f64, 2>::new(res[(1, 0)] + res[(2, 1)], res[(1, 1)] + res[(2, 2)]);
                    let mut result = result[index].lock().unwrap();
                    result[0] += node_force.x;
                    result[1] += node_force.y;
                }
            });
        
            grid.par_iter()
                .enumerate()
                .map(|(index, _)| {
                    (
                        (index % (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        (index / (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        {
                            let result = result[index].lock().unwrap();
                            Vector2::<f64>::new(result[0], result[1])
                        }
                    )
                })
                .filter(|(x, y, _)| 4. <= *x && *x < 6. && 4. <= *y && *y < 6.)
                .collect::<Vec<_>>()
        }
    }.into()
}

#[proc_macro]
pub fn scheme_laplacian_velocity(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let q = 1;
    let func_name = format_ident!("scheme_laplacian_velocity_{}", p);
    let scale_p_q = format_ident!("scale_{}_{}", p, q);
    let c_p_q = format_ident!("c_{}_{}", p, q);
    let (dx, dy) = multi_index_vec(q);
    let i = 0..(dx.len());

    quote! {
        fn #func_name(settings: &Settings) -> Vec<(f64, f64, SVector<f64, 2>)> {
            let (mut particles, mut grid, periodic_boundary_rect) = new_for_taylor_green(settings);

            mlsmpm::g2p_lsmps_1st(
                &mut particles,
                &grid,
                settings,
                periodic_boundary_rect.clone(),
            );

            grid.par_iter_mut().for_each(|node| node.reset());

            let period_bounds = vec![];
            let periodic_boundary_rect = Some(periodic_boundary_rect);
        
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_params!(#p);
            mlsmpm_macro::compact_lsmps_func!(#p, #q);
            mlsmpm_macro::compact_scale!(#p, #q);
        
            let rs = settings.cell_width();
            let scale = #scale_p_q(rs);
        
            let nodes = {
                let mut nodes = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    nodes.push(std::sync::Mutex::new(LsmpsParams {
                        m: SMatrix::zeros(),
                        f_vel: SMatrix::zeros(),
                        f_stress: SMatrix::zeros(),
                        f_pressure: SVector::zeros(),
                        count: 0,
                    }));
                }
                nodes
            };
        
            parallel!(settings, particles, |p| {
                for node in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();

                    #({
                        let (i, bx, by) = (#i, #dx, #dy);
                        params.f_vel += weight * poly_r_ij.kronecker(&p.x_lsmps().column(i).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                    })*

                    params.count += 1;
                }
            });
        
            let result = {
                let mut result = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    result.push(Mutex::new([0.; 2]));
                }
                result
            };

            nodes.par_iter().enumerate().for_each(|(index, params)| {
                let params = params.lock().unwrap();

                if params.count == 0 {
                    return;
                }

                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    let res = scale * m_inverse * params.f_vel;
        
                    let node_force = settings.dynamic_viscosity
                        * settings.rho_0
                        * vector![res[(3, 0)] + res[(5, 0)], res[(3, 1)] + res[(5, 1)]];
                    let mut result = result[index].lock().unwrap();
                    result[0] += node_force.x;
                    result[1] += node_force.y;
                }
            });
        
            grid.par_iter()
                .enumerate()
                .map(|(index, _)| {
                    (
                        (index % (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        (index / (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        {
                            let result = result[index].lock().unwrap();
                            Vector2::<f64>::new(result[0], result[1])
                        }
                    )
                })
                .filter(|(x, y, _)| 4. <= *x && *x < 6. && 4. <= *y && *y < 6.)
                .collect::<Vec<_>>()
        }        
    }.into()
}

#[proc_macro]
pub fn scheme_laplacian_velocity_g2p_2nd(input: TokenStream) -> TokenStream {
    let p = parse_one_usize(input);
    let q = 2;
    let func_name = format_ident!("scheme_laplacian_velocity_{}_g2p_2nd", p);
    let scale_p_q = format_ident!("scale_{}_{}", p, q);
    let c_p_q = format_ident!("c_{}_{}", p, q);
    let (dx, dy) = multi_index_vec(q);
    let i = 0..(dx.len());

    quote! {
        fn #func_name(settings: &Settings) -> Vec<(f64, f64, SVector<f64, 2>)> {
            let (mut particles, mut grid, periodic_boundary_rect) = new_for_taylor_green(settings);

            mlsmpm::g2p_lsmps_2nd(
                &mut particles,
                &grid,
                settings,
                periodic_boundary_rect.clone(),
            );

            grid.par_iter_mut().for_each(|node| node.reset());

            let period_bounds = vec![];
            let periodic_boundary_rect = Some(periodic_boundary_rect);
        
            mlsmpm_macro::lsmps_poly!(#p);
            mlsmpm_macro::lsmps_params!(#p);
            mlsmpm_macro::compact_lsmps_func!(#p, #q);
            mlsmpm_macro::compact_scale!(#p, #q);
        
            let rs = settings.cell_width();
            let scale = #scale_p_q(rs);
        
            let nodes = {
                let mut nodes = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    nodes.push(std::sync::Mutex::new(LsmpsParams {
                        m: SMatrix::zeros(),
                        f_vel: SMatrix::zeros(),
                        f_stress: SMatrix::zeros(),
                        f_pressure: SVector::zeros(),
                        count: 0,
                    }));
                }
                nodes
            };
        
            parallel!(settings, particles, |p| {
                for node in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
                    let mut params = nodes[node.index].lock().unwrap();
        
                    let r_ij = -node.dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = node.weight;
        
                    params.m += weight * poly_r_ij * poly_r_ij.transpose();

                    #({
                        let (i, bx, by) = (#i, #dx, #dy);
                        params.f_vel += weight * poly_r_ij.kronecker(&p.x_lsmps().column(i).transpose()) * #c_p_q(bx, by) * (-node.dist.x).powi(bx as i32) * (-node.dist.y).powi(by as i32);
                    })*

                    params.count += 1;
                }
            });
        
            let result = {
                let mut result = Vec::with_capacity(grid.len());
                for _ in 0..grid.len() {
                    result.push(Mutex::new([0.; 2]));
                }
                result
            };

            nodes.par_iter().enumerate().for_each(|(index, params)| {
                let params = params.lock().unwrap();

                if params.count == 0 {
                    return;
                }

                if let Some(m_inverse) = (params.m + SMatrix::identity() * 0.).try_inverse() {
                    let res = scale * m_inverse * params.f_vel;
        
                    let node_force = settings.dynamic_viscosity
                        * settings.rho_0
                        * vector![res[(3, 0)] + res[(5, 0)], res[(3, 1)] + res[(5, 1)]];
                    let mut result = result[index].lock().unwrap();
                    result[0] += node_force.x;
                    result[1] += node_force.y;
                }
            });
        
            grid.par_iter()
                .enumerate()
                .map(|(index, _)| {
                    (
                        (index % (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        (index / (settings.grid_width + 1)) as f64 * settings.cell_width(),
                        {
                            let result = result[index].lock().unwrap();
                            Vector2::<f64>::new(result[0], result[1])
                        }
                    )
                })
                .filter(|(x, y, _)| 4. <= *x && *x < 6. && 4. <= *y && *y < 6.)
                .collect::<Vec<_>>()
        }        
    }.into()
}