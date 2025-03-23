use num_complex::Complex;

monarch_derive::generate_powers_of_two!();

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
