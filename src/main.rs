extern crate nalgebra as na;
use na::Vector2;
use itertools::izip;

fn evolve_cell(f_in: &[f64], f_out_field: &mut [f64], f_out_idx: &[usize],
    weights: &[f64], velocities: &[Vector2<f64>], inv_tau: f64)
{
    let density: f64 = f_in.iter().sum();
    let velocity: Vector2<f64> = izip!(f_in.iter(), velocities.iter())
        .map(|(f, c)| {*f / density * *c}).sum();

    let vv = velocity.dot(&velocity);
    for (i, f, w, c) in izip!(
        f_out_idx.iter(), f_in.iter(), weights.iter(), velocities.iter())
    {
        let cv = c.dot(&velocity);
        f_out_field[*i] = -inv_tau * (f - density*w *
            (1.0 + 3.0*cv + 4.5*cv.powi(2) - 1.5*vv));
    }
}

const W: usize = 400;
const H: usize = 100;

const VELOCITIES: [Vector2<f64>; 9] = [
    Vector2::new(0., 0.),
    Vector2::new(0., 1.),
    Vector2::new(0.,-1.),
    Vector2::new(1., 0.),
    Vector2::new(-1., 0.),
    Vector2::new(1., 1.),
    Vector2::new(1., -1.),
    Vector2::new(-1., 1.),
    Vector2::new(-1., -1.)
];

const WEIGHTS: [f64;9] = [4./9.,1./9.,1./36.,1./9.,1./36.,1./9.,1./36.,1./9.,1./36.];

fn create_grid_connections() -> Vec<usize>
{
    let mut ret = Vec::<usize>::with_capacity(W*H*VELOCITIES.len());

    for y in 0..H {
        for x in 0..W {
            for i in 0..VELOCITIES.len() {
                let v = VELOCITIES[i];
                let y = (y + v.y as usize) % H;
                let x = (x + v.x as usize) % W;
                ret.push(y*W*VELOCITIES.len() + x*VELOCITIES.len() + i);
            }
        }
    }

    ret
}

fn main()
{
    let grid = create_grid_connections();

    let mut input = vec![0.0f64; grid.len()];
    let mut output = vec![0.0f64; grid.len()];
    loop {
        for (f_in, f_out_idx) in izip!(
            input.chunks(VELOCITIES.len()),
            grid.chunks(VELOCITIES.len())
        ) {
            evolve_cell(f_in, &mut output[..],
                f_out_idx, &WEIGHTS[..], &VELOCITIES[..], 1./0.6);
        }

        // Double buffer swap:
        std::mem::swap(&mut input, &mut output);
    }
}
