use rand::Rng;
use std::env;
use std::f64::consts::PI;
use std::thread;

fn fft(x: &mut [f64]) {
    let (mut real, mut imag, mut t_real, mut t_imag): (f64, f64, f64, f64);
    
    let n: usize = x.len();

    if n == 2 {
        return;
    }

    let mut even: Vec<f64> = vec![0.0; n / 2];
    let mut odd: Vec<f64> = vec![0.0; n / 2];

    for i in (0..n / 2).step_by(2) {
        even[i] = x[2 * i];
        even[i + 1] = x[2 * i + 1];
        odd[i] = x[2 * i + 2];
        odd[i + 1] = x[2 * i + 3];
    }

    fft(&mut even);
    fft(&mut odd);

    for k in (0..n / 2).step_by(2) {
        imag = -2.0 * PI * (k as f64) / (n as f64);
        real = f64::cos(imag);
        imag = f64::sin(imag);

        t_real = real * odd[k] - imag * odd[k + 1];
        t_imag = real * odd[k + 1] + imag * odd[k];

        x[k] = even[k] + t_real;
        x[k + 1] = even[k + 1] + t_imag;

        x[k + n / 2] = even[k] - t_real;
        x[k + n / 2 + 1] = even[k + 1] - t_imag;
    }
}

fn fft_loop(size: usize, loops: usize) {
    let mut rng = rand::thread_rng();

    let mut x: Vec<f64>;

    for _ in 0..loops {
        x = (0..2 * size).map(|_| rng.gen_range(0.0..1.0)).collect();
        fft(&mut x);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 4 {
        println!("Usage: ./fft_rs <Input Size> <Loops> <Threads>");
        return;
    }

    let fft_size: usize = 1 << args[1].parse::<usize>().unwrap();
    let num_loops: usize = args[2].parse::<usize>().unwrap();
    let num_threads: usize = args[3].parse::<usize>().unwrap();

    if num_threads == 1 {
        fft_loop(fft_size, num_loops);
        std::process::exit(0);
    }

    let handles: Vec<_> = (0..num_threads)
        .map(|_| {
            let fs: usize = fft_size.clone();
            let nl: usize = num_loops.clone();
            thread::spawn(move || {
                fft_loop(fs, nl);
            })
        })
        .collect();

    for handle in handles {
        handle.join().unwrap();
    }
}
