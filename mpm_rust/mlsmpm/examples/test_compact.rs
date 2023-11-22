fn main() {
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
            (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by))
                * factorial(p)
                / factorial(p + q)
        } else {
            (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by))
                * q as f64
                * factorial(p + q - b)
                / factorial(p + q)
        }
    }

    fn s(ax: usize, ay: usize, rs: f64, p: usize, q: usize) -> f64 {
        let a = ax + ay;
        let sum = [(0, 0), (1, 0), (0, 1)]
            .iter()
            .filter(|(bx, by)| bx + by <= q && *bx <= ax && *by <= ay)
            .map(|&(bx, by)| c(bx, by, p, q) / (factorial(ax - bx) * factorial(ay - by)))
            .sum::<f64>();
        1. / (sum * rs.powi(a as i32))
    }

    println!("{}", s(1, 0, 1., 2, 1));
    println!("{}", s(0, 1, 1., 2, 1));
    println!("{}", s(2, 0, 1., 2, 1));
    println!("{}", s(1, 1, 1., 2, 1));
    println!("{}", s(0, 2, 1., 2, 1));
}
