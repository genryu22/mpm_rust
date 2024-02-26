fn main() {
    fn test(kernel_func: fn(f64) -> f64) {
        let split = 50;
        let dx = 3. / split as f64;
        for i in 0..=split {
            let x = -1.5 + i as f64 * dx;
            println!("{}, {}", x, kernel_func(x));
        }
    }

    test(quadratic_kernel);
    test(left_half_quadratic_kernel);
    test(left_wall_quadratic_kernel);
    test(right_half_quadratic_kernel);
    test(right_wall_quadratic_kernel);
}

fn quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1.5, -0.5, 0.5, 1.5], x)
}

fn left_half_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1., -0.5, 0.5, 1.5], x)
}

fn left_wall_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![0., 0., 0.5, 1.5], x)
        + b_spline_basis(0, 2, &vec![0., 0., 0., 0.5], x)
}

fn right_half_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1.5, -0.5, 0.5, 1.], x)
}

fn right_wall_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1.5, -0.5, 0., 0.], x)
        + b_spline_basis(0, 2, &vec![-0.5, 0., 0., 0.], x)
}

fn b_spline_basis(j: usize, k: usize, knot: &Vec<f64>, t: f64) -> f64 {
    if knot.len() as i32 - k as i32 - 2 < 0 || j > knot.len() - 2 {
        panic!();
    }

    if k == 0 {
        if knot[j] <= t && t < knot[j + 1] {
            1.
        } else {
            0.
        }
    } else {
        fn w(j: usize, k: usize, knot: &Vec<f64>, t: f64) -> f64 {
            if knot[j + k] == knot[j] {
                0.
            } else {
                (t - knot[j]) / (knot[j + k] - knot[j])
            }
        }

        w(j, k, knot, t) * b_spline_basis(j, k - 1, knot, t)
            + (1. - w(j + 1, k, knot, t)) * b_spline_basis(j + 1, k - 1, knot, t)
    }
}
