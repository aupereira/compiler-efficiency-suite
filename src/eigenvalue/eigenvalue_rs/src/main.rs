use rand::Rng;
use std::env;
use std::mem;
use std::thread;

fn sgn(x: f64) -> f64 {
    return if x >= 0.0 { 1.0 } else { -1.0 };
}

fn rand_symmetric(n: usize, rng: &mut rand::rngs::ThreadRng, a: &mut Vec<Vec<f64>>) {
    for i in 0..n {
        for j in i + 1..n {
            a[i][j] = rng.gen();
            a[j][i] = a[i][j];
        }
    }
}

fn matrix_multiply(
    n: usize,
    a: &Vec<Vec<f64>>,
    b: &Vec<Vec<f64>>,
    c: &mut Vec<Vec<f64>>,
    btmp: &mut Vec<Vec<f64>>,
) {
    let mut sum: f64;

    for i in 0..n {
        for j in 0..n {
            btmp[i][j] = b[j][i];
            // unsafe {btmp[i][j] = *b.get_unchecked(j).get_unchecked(i);}
        }
    }

    for i in 0..n {
        for j in 0..n {
            sum = 0.0;
            for k in 0..n {
                sum += a[i][k] * btmp[j][k];
                // unsafe {sum += *a.get_unchecked(i).get_unchecked(k) * *btmp.get_unchecked(j).get_unchecked(k);}
            }
            c[i][j] = sum;
            // unsafe {*c.get_unchecked_mut(i).get_unchecked_mut(j) = sum;}
        }
    }
}

fn norm2(n: usize, v: &Vec<f64>) -> f64 {
    let mut sum: f64 = 0.0;

    for i in 0..n {
        sum += v[i] * v[i];
    }

    return sum.sqrt();
}

fn s2opai(len: usize, v: &Vec<f64>, a: &mut Vec<Vec<f64>>) {
    for i in 0..len {
        for j in 0..len {
            a[i][j] = -2.0 * v[i] * v[j];
        }
    }

    for i in 0..len {
        a[i][i] += 1.0;
    }
}

fn v_div_s(len: usize, v: &mut Vec<f64>, s: f64) {
    for i in 0..len {
        v[i] /= s;
    }
}

fn tridiagonalize(
    n: usize,
    a: &mut Vec<Vec<f64>>,
    p: &mut Vec<Vec<f64>>,
    temp: &mut Vec<Vec<f64>>,
    tpose: &mut Vec<Vec<f64>>,
    v: &mut Vec<f64>,
) {
    let mut alpha: f64;
    let mut r: f64;

    for k in 0..n - 2 {
        for j in 0..k + 1 {
            v[j] = 0.0;
        }
        for j in k + 1..n {
            v[j] = a[j][k];
        }

        alpha = -sgn(v[k + 1]) * norm2(n, &v);
        r = -2.0 * (0.5 * ((alpha * alpha) - (alpha * v[k + 1]))).sqrt();

        v[k + 1] = (v[k + 1] - alpha) / r;
        for j in k + 2..n {
            v[j] /= r;
        }

        s2opai(n, &v, p);
        matrix_multiply(n, p, a, temp, tpose);
        matrix_multiply(n, temp, p, a, tpose);
    }
}

fn householder(
    len: usize,
    a: &mut Vec<Vec<f64>>,
    a1: &mut Vec<Vec<f64>>,
    q: &mut Vec<Vec<f64>>,
    h: &mut Vec<Vec<f64>>,
    temp: &mut Vec<Vec<f64>>,
    tpose: &mut Vec<Vec<f64>>,
    v: &mut Vec<f64>,
) {
    *a1 = a.clone();

    let s: f64 = a1[len - 1][len - 1];
    let t: f64 = a1[len - 2][len - 2];
    let x: f64 = a1[len - 2][len - 1];
    let d: f64 = (t - s) / 2.0;
    let shift: f64 = s - (x * x) / (d + sgn(d) * (d * d + x * x).sqrt());

    for i in 0..len {
        a1[i][i] -= shift;
    }

    for k in 0..len - 1 {
        for j in 0..k {
            v[j] = 0.0;
        }
        for j in k..len {
            v[j] = a1[j][k];
        }

        v[k] += norm2(len, v) * sgn(v[k + 1]);
        v_div_s(len, v, norm2(len, v));

        if k == 0 {
            s2opai(len, &v, q);
            matrix_multiply(len, q, a1, temp, tpose);
            mem::swap(a1, temp);
        } else {
            s2opai(len, &v, h);
            matrix_multiply(len, h, a1, temp, tpose);
            mem::swap(a1, temp);
            matrix_multiply(len, q, h, temp, tpose);
            mem::swap(q, temp);
        }
    }

    matrix_multiply(len, a1, q, a, tpose);

    for i in 0..len {
        a[i][i] += shift;
    }
}

fn eig_qr(
    iter: usize,
    a: &mut Vec<Vec<f64>>,
    eig_vals: &mut Vec<f64>,
    t1: &mut Vec<Vec<f64>>,
    t2: &mut Vec<Vec<f64>>,
    t3: &mut Vec<Vec<f64>>,
    t4: &mut Vec<Vec<f64>>,
    t5: &mut Vec<Vec<f64>>,
    v: &mut Vec<f64>,
) {
    let mut ind: usize = a.len();

    while ind > 1 {
        for _ in 0..iter {
            householder(ind, a, t1, t2, t3, t4, t5, v);
        }
        eig_vals[ind - 1] = a[ind - 1][ind - 1];
        ind -= 1;
    }

    eig_vals[0] = a[0][0];
}

fn eigen_loop(size: usize, iter: usize, loops: usize) {
    let mut rng = rand::thread_rng();

    let mut a: Vec<Vec<f64>> = vec![vec![0.0; size]; size];

    let mut t1: Vec<Vec<f64>> = vec![vec![0.0; size]; size];
    let mut t2: Vec<Vec<f64>> = vec![vec![0.0; size]; size];
    let mut t3: Vec<Vec<f64>> = vec![vec![0.0; size]; size];
    let mut t4: Vec<Vec<f64>> = vec![vec![0.0; size]; size];
    let mut t5: Vec<Vec<f64>> = vec![vec![0.0; size]; size];
    let mut v: Vec<f64> = vec![0.0; size];

    let mut eig_vals: Vec<f64> = vec![0.0; size];

    for _ in 0..loops {
        rand_symmetric(size, &mut rng, &mut a);
        tridiagonalize(size, &mut a, &mut t1, &mut t2, &mut t3, &mut v);
        eig_qr(
            iter,
            &mut a,
            &mut eig_vals,
            &mut t1,
            &mut t2,
            &mut t3,
            &mut t4,
            &mut t5,
            &mut v,
        );
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 5 {
        println!("Usage: ./eigenvalue_rs <Matrix Size> <Convergence Iterations> <Loops> <Threads>");
        return;
    }

    let m_size: usize = args[1].parse::<usize>().unwrap();
    let c_iter: usize = args[2].parse::<usize>().unwrap();
    let num_loops: usize = args[3].parse::<usize>().unwrap();
    let num_threads: usize = args[4].parse::<usize>().unwrap();

    if num_threads == 1 {
        eigen_loop(m_size, c_iter, num_loops);
        std::process::exit(0);
    }

    let handles: Vec<_> = (0..num_threads)
        .map(|_| {
            let size: usize = m_size.clone();
            let iter: usize = c_iter.clone();
            let loops: usize = num_loops.clone();
            thread::spawn(move || {
                eigen_loop(size, iter, loops);
            })
        })
        .collect();

    for handle in handles {
        handle.join().unwrap();
    }
}
