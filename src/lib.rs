use std::f32::consts::PI;

use num_complex::Complex;

pub fn butterfly(x: Vec<Complex<f32>>) -> Vec<Complex<f32>> {
    let n = x.len();
    if n <= 1 {
        x
    } else {
        let even: Vec<_> = butterfly(x.iter().step_by(2).cloned().collect());
        let odd: Vec<_> = butterfly(x.iter().skip(1).step_by(2).cloned().collect());

        let t: Vec<_> = (0..(n / 2))
            .map(|k| Complex::exp(-2.0 * Complex::<f32>::i() * PI * k as f32 / n as f32) * odd[k])
            .collect();

        [
            (0..(n / 2)).map(|k| even[k] + t[k]).collect::<Vec<_>>(),
            (0..(n / 2)).map(|k| even[k] - t[k]).collect::<Vec<_>>(),
        ]
        .concat()
    }
}

pub fn butterfly1(x: [Complex<f32>; 1]) -> [Complex<f32>; 1] {
    x
}

pub fn butterfly2(x: [Complex<f32>; 2]) -> [Complex<f32>; 2] {
    let n = 2;

    let even: [Complex<f32>; 1] = butterfly1([x[0]]);
    let odd: [Complex<f32>; 1] = butterfly1([x[1]]);

    let t: [Complex<f32>; 1] =
        [Complex::exp(-2.0 * Complex::<f32>::i() * PI * 0.0 / n as f32) * odd[0]];

    [even[0] + t[0], even[0] - t[0]]
}

pub fn butterfly4(x: [Complex<f32>; 4]) -> [Complex<f32>; 4] {
    let n = 4;

    let even: [Complex<f32>; 2] = butterfly2([x[0], x[2]]);
    let odd: [Complex<f32>; 2] = butterfly2([x[1], x[3]]);

    let t: [Complex<f32>; 2] = [
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 0.0 / n as f32) * odd[0],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 1.0 / n as f32) * odd[1],
    ];

    [
        even[0] + t[0],
        even[1] + t[1],
        even[0] - t[0],
        even[1] - t[1],
    ]
}

pub fn butterfly8(x: [Complex<f32>; 8]) -> [Complex<f32>; 8] {
    let n = 8;

    let even: [Complex<f32>; 4] = butterfly4([x[0], x[2], x[4], x[6]]);
    let odd: [Complex<f32>; 4] = butterfly4([x[1], x[3], x[5], x[7]]);

    let t: [Complex<f32>; 4] = [
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 0.0 / n as f32) * odd[0],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 1.0 / n as f32) * odd[1],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 2.0 / n as f32) * odd[2],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 3.0 / n as f32) * odd[3],
    ];

    [
        even[0] + t[0],
        even[1] + t[1],
        even[2] + t[2],
        even[3] + t[3],
        even[0] - t[0],
        even[1] - t[1],
        even[2] - t[2],
        even[3] - t[3],
    ]
}

pub fn butterfly16(x: [Complex<f32>; 16]) -> [Complex<f32>; 16] {
    let n = 16;

    let even: [Complex<f32>; 8] = butterfly8([x[0], x[2], x[4], x[6], x[8], x[10], x[12], x[14]]);
    let odd: [Complex<f32>; 8] = butterfly8([x[1], x[3], x[5], x[7], x[9], x[11], x[13], x[15]]);

    let t: [Complex<f32>; 8] = [
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 0.0 / n as f32) * odd[0],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 1.0 / n as f32) * odd[1],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 2.0 / n as f32) * odd[2],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 3.0 / n as f32) * odd[3],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 4.0 / n as f32) * odd[4],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 5.0 / n as f32) * odd[5],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 6.0 / n as f32) * odd[6],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 7.0 / n as f32) * odd[7],
    ];

    [
        even[0] + t[0],
        even[1] + t[1],
        even[2] + t[2],
        even[3] + t[3],
        even[4] + t[4],
        even[5] + t[5],
        even[6] + t[6],
        even[7] + t[7],
        even[0] - t[0],
        even[1] - t[1],
        even[2] - t[2],
        even[3] - t[3],
        even[4] - t[4],
        even[5] - t[5],
        even[6] - t[6],
        even[7] - t[7],
    ]
}

#[cfg(test)]
mod tests {
    use num_complex::Complex;

    use crate::*;

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
}
