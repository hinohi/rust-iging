extern crate rand;

use rand::prelude::*;
use std::f64;
use std::string::String;
use std::fmt::Write;

const L: usize = 32;
const N: usize = L * L;

struct System {
    cell: [i32; N],
    temperature: f64,
}

impl System {
    fn new() -> System {
        let t = 10.0;
        System {
            cell: [1; N],
            temperature: t,
        }
    }

    fn update(&mut self, rng: &mut ThreadRng) {
        let mut order: Vec<usize> = (1..N).collect();
        let beta = 1.0 / self.temperature;
        rng.shuffle(&mut order);
        for i in order {
            let y = i / L;
            let x = i % L;
            let y1 = (y + L - 1) % L;
            let y2 = (y + 1) % L;
            let x1 = (x + L - 1) % L;
            let x2 = (x + 1) % L;
            let local_mag = self.cell[y1 * L + x]
                + self.cell[y * L + x1]
                + self.cell[y * L + x2]
                + self.cell[y2 + L + x];
            let de = 2 * local_mag * self.cell[y * L + x];
            // Metropolis method
            if de < 0 || rng.gen::<f64>() < (-de as f64 * beta).exp() {
                self.cell[y * L + x] = -self.cell[y * L + x];
            }
        }
    }

    fn set_temperature(&mut self, temperature: f64) {
        self.temperature = temperature;
    }

    fn calc_quantity(&self) -> (f64, f64) {
        let mut enery = 0;
        let mut mag = 0;
        for i in 0..N {
            let y = i / L;
            let x = i % L;
            let y2 = (y + 1) % L;
            let x2 = (x + 1) % L;
            enery += self.cell[i] * (self.cell[y * L + x2] + self.cell[y2 * L + x]);
            mag += self.cell[i];
        }
        (enery as f64 / N as f64, mag as f64 / N as f64)
    }
}

struct Quantity {
    temperature: f64,
    sample: u32,
    energy: f64,
    energy2: f64,
    magnetic_charge: f64,
    magnetic_charge2: f64,
    magnetic_charge4: f64,
}

impl Quantity {
    fn new(temperature: f64) -> Quantity {
        Quantity {
            temperature,
            sample: 0,
            energy: 0.0,
            energy2: 0.0,
            magnetic_charge: 0.0,
            magnetic_charge2: 0.0,
            magnetic_charge4: 0.0,
        }
    }

    fn add(&mut self, enery: f64, magnetic_charge: f64) {
        self.sample += 1;
        self.energy += enery;
        self.energy2 += enery * enery;
        self.magnetic_charge += magnetic_charge;
        let mag2 = magnetic_charge * magnetic_charge;
        self.magnetic_charge2 += mag2;
        self.magnetic_charge4 += mag2 * mag2;
    }

    fn dump(&self) -> String {
        let mut s = String::new();
        let energy = self.energy / self.sample as f64;
        let energy2 = self.energy2 / self.sample as f64;
        let energy_error = if energy2 - energy * energy > 0.0 {
            ((energy2 - energy * energy) / (self.sample - 1) as f64).sqrt()
        } else {
            0.0
        };
        let mag = self.magnetic_charge / self.sample as f64;
        let mag2 = self.magnetic_charge2 / self.sample as f64;
        let mag4 = self.magnetic_charge4 / self.sample as f64;
        let mag_error = if mag2 - mag * mag > 0.0 {
            ((mag2 - mag * mag) / (self.sample - 1) as f64).sqrt()
        } else {
            0.0
        };
        let binder = (3.0 - mag4 / (mag2 * mag2)) * 0.5;
        write!(
            s,
            "{} {} {} {} {} {}",
            self.temperature, energy, energy_error, mag, mag_error, binder
        ).unwrap();
        s
    }
}

fn main() {
    let mut rng = thread_rng();
    let mut system = System::new();
    let max_temp = 5.0;
    let min_temp = 0.1;
    let temp_step = 100;
    let n_sample = 10_000;
    let sample_step = 10;

    system.set_temperature(max_temp);
    for _ in 0..N * 10 {
        system.update(&mut rng);
    }

    for step in 0..temp_step + 1 {
        let temperature = max_temp - (max_temp - min_temp) / temp_step as f64 * step as f64;
        system.set_temperature(temperature);
        let mut q = Quantity::new(temperature);
        for _ in 0..n_sample {
            for _ in 0..sample_step {
                system.update(&mut rng);
            }
            let (e, m) = system.calc_quantity();
            q.add(e, m);
        }
        println!("{}", q.dump());
    }
}
