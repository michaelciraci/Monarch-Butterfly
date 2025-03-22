use num_complex::Complex;

monarch_derive::generate_powers_of_two!();

#[inline]
pub fn butterfly5(x: [Complex<f32>; 5]) -> [Complex<f32>; 5] {
    let twiddle1 = compute_twiddle(1, 5);
    let twiddle2 = compute_twiddle(2, 5);
    let x14p = x[1] + x[4];
    let x14n = x[1] - x[4];
    let x23p = x[2] + x[3];
    let x23n = x[2] - x[3];
    let sum = x[0] + x14p + x23p;
    let b14re_a = x[0].re + twiddle1.re * x14p.re + twiddle2.re * x23p.re;
    let b14re_b = twiddle1.im * x14n.im + twiddle2.im * x23n.im;
    let b23re_a = x[0].re + twiddle2.re * x14p.re + twiddle1.re * x23p.re;
    let b23re_b = twiddle2.im * x14n.im + -twiddle1.im * x23n.im;
    let b14im_a = x[0].im + twiddle1.re * x14p.im + twiddle2.re * x23p.im;
    let b14im_b = twiddle1.im * x14n.re + twiddle2.im * x23n.re;
    let b23im_a = x[0].im + twiddle2.re * x14p.im + twiddle1.re * x23p.im;
    let b23im_b = twiddle2.im * x14n.re + -twiddle1.im * x23n.re;
    let out1re = b14re_a - b14re_b;
    let out1im = b14im_a + b14im_b;
    let out2re = b23re_a - b23re_b;
    let out2im = b23im_a + b23im_b;
    let out3re = b23re_a + b23re_b;
    let out3im = b23im_a - b23im_b;
    let out4re = b14re_a + b14re_b;
    let out4im = b14im_a - b14im_b;

    [
        sum,
        Complex::new(out1re, out1im),
        Complex::new(out2re, out2im),
        Complex::new(out3re, out3im),
        Complex::new(out4re, out4im),
    ]
}

#[inline]
fn compute_twiddle(index: usize, fft_len: usize) -> Complex<f32> {
    let constant = -2.0 * std::f32::consts::PI / fft_len as f32;
    let angle = constant * index as f32;

    Complex::new(angle.cos(), angle.sin())
}

#[cfg(test)]
mod tests {
    use num_complex::Complex;

    use crate::*;

    fn butterfly(x: Vec<Complex<f32>>) -> Vec<Complex<f32>> {
        let n = x.len();
        if n <= 1 {
            x
        } else {
            let even: Vec<_> = butterfly(x.iter().step_by(2).cloned().collect());
            let odd: Vec<_> = butterfly(x.iter().skip(1).step_by(2).cloned().collect());

            let t: Vec<_> = (0..(n / 2))
                .map(|k| {
                    Complex::exp(
                        -2.0 * Complex::<f32>::i() * std::f32::consts::PI * k as f32 / n as f32,
                    ) * odd[k]
                })
                .collect();

            [
                (0..(n / 2)).map(|k| even[k] + t[k]).collect::<Vec<_>>(),
                (0..(n / 2)).map(|k| even[k] - t[k]).collect::<Vec<_>>(),
            ]
            .concat()
        }
    }

    #[test]
    fn test_8() {
        let a = butterfly((1..9).map(|i| Complex::new(i as f32, 0.0)).collect());
        assert_eq!(a[0], Complex::new(36.0, 0.0));
        assert_eq!(a[6], Complex::new(-4.0, -4.0));
    }

    #[test]
    fn test_butterfly_2() {
        assert_eq!(
            butterfly2([Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)]).to_vec(),
            butterfly(vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)])
        );
    }

    #[test]
    fn test_butterfly_4() {
        assert_eq!(
            butterfly4([
                Complex::new(1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(3.0, 0.0),
                Complex::new(4.0, 0.0)
            ])
            .to_vec(),
            butterfly(vec![
                Complex::new(1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(3.0, 0.0),
                Complex::new(4.0, 0.0),
            ])
        );
    }

    #[test]
    fn test_butterfly_1024() {
        let v: Vec<_> = (0..1024)
            .map(|i: i32| Complex::new(i as f32, i as f32))
            .collect();
        let a: [Complex<f32>; 1024] = (0..1024)
            .map(|i| Complex::new(i as f32, i as f32))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();

        let ba = butterfly1024(a);
        let bv = butterfly(v);

        assert_eq!(&ba, &*bv);
    }
}
