extern crate proc_macro;
use proc_macro::TokenStream;
use quote::quote;
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
    let input: usize = input.to_string().parse().unwrap();
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
    let input: usize = input.to_string().parse().unwrap();
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
    let input: usize = input.to_string().parse().unwrap();
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
    let input: usize = input.to_string().parse().unwrap();
    let (mut res, size) = multi_index(input);

    quote! {
        struct LsmpsParams {
            m: SMatrix<f64, #size, #size>,
            f_vel: SMatrix<f64, #size, 2>,
        }
    }
    .into()
}

fn factorial(num: usize) -> f64 {
    match num {
        0 | 1 => 1.,
        _ => factorial(num - 1) * num as f64,
    }
}

fn c(bx: usize, by: usize, p: usize, q: usize) -> f64 {
    let b = bx + by;
    if b == 0 {
        1.
    } else if b == q {
        (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by)) * factorial(p)
            / factorial(p + q)
    } else {
        (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by))
            * q as f64
            * factorial(p + q - b)
            / factorial(p + q)
    }
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
pub fn compact_lsmps_s(input: TokenStream) -> TokenStream {
    let (p, q) = parse_two_usize(input);
    println!("{}, {}", p, q);

    let (mut res, size) = multi_index(p);

    quote! {
        struct LsmpsParams {
            m: SMatrix<f64, #size, #size>,
            f_vel: SMatrix<f64, #size, 2>,
        }
    }
    .into()
}
