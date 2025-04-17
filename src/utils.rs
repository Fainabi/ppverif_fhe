pub fn normalize(fs: &mut [f32]) {
    let sqrt_norm: f32 = fs.iter().map(|&x| x * x).sum::<f32>().sqrt();
    fs.iter_mut().for_each(|v| *v /= sqrt_norm);
}
