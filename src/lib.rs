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

#[inline]
pub fn butterfly1(x: [Complex<f32>; 1]) -> [Complex<f32>; 1] {
    x
}

#[inline]
pub fn butterfly2(x: [Complex<f32>; 2]) -> [Complex<f32>; 2] {
    let n = 2;

    let even: [Complex<f32>; 1] = butterfly1([x[0]]);
    let odd: [Complex<f32>; 1] = butterfly1([x[1]]);

    let t: [Complex<f32>; 1] =
        [Complex::exp(-2.0 * Complex::<f32>::i() * PI * 0.0 / n as f32) * odd[0]];

    [even[0] + t[0], even[0] - t[0]]
}

#[inline]
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

#[inline]
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

#[inline]
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

#[inline]
pub fn butterfly32(x: [Complex<f32>; 32]) -> [Complex<f32>; 32] {
    let n = 32;

    let even: [Complex<f32>; 16] = butterfly16([
        x[0], x[2], x[4], x[6], x[8], x[10], x[12], x[14], x[16], x[18], x[20], x[22], x[24],
        x[26], x[28], x[30],
    ]);
    let odd: [Complex<f32>; 16] = butterfly16([
        x[1], x[3], x[5], x[7], x[9], x[11], x[13], x[15], x[17], x[19], x[21], x[23], x[25],
        x[27], x[29], x[31],
    ]);

    let t: [Complex<f32>; 16] = [
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 0.0 / n as f32) * odd[0],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 1.0 / n as f32) * odd[1],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 2.0 / n as f32) * odd[2],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 3.0 / n as f32) * odd[3],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 4.0 / n as f32) * odd[4],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 5.0 / n as f32) * odd[5],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 6.0 / n as f32) * odd[6],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 7.0 / n as f32) * odd[7],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 8.0 / n as f32) * odd[8],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 9.0 / n as f32) * odd[9],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 10.0 / n as f32) * odd[10],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 11.0 / n as f32) * odd[11],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 12.0 / n as f32) * odd[12],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 13.0 / n as f32) * odd[13],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 14.0 / n as f32) * odd[14],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 15.0 / n as f32) * odd[15],
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
        even[8] + t[8],
        even[9] + t[9],
        even[10] + t[10],
        even[11] + t[11],
        even[12] + t[12],
        even[13] + t[13],
        even[14] + t[14],
        even[15] + t[15],
        even[0] - t[0],
        even[1] - t[1],
        even[2] - t[2],
        even[3] - t[3],
        even[4] - t[4],
        even[5] - t[5],
        even[6] - t[6],
        even[7] - t[7],
        even[8] - t[8],
        even[9] - t[9],
        even[10] - t[10],
        even[11] - t[11],
        even[12] - t[12],
        even[13] - t[13],
        even[14] - t[14],
        even[15] - t[15],
    ]
}

#[inline]
pub fn butterfly64(x: [Complex<f32>; 64]) -> [Complex<f32>; 64] {
    let n = 32;

    let even: [Complex<f32>; 32] = butterfly32([
        x[0], x[2], x[4], x[6], x[8], x[10], x[12], x[14], x[16], x[18], x[20], x[22], x[24],
        x[26], x[28], x[30], x[32], x[34], x[36], x[38], x[40], x[42], x[44], x[46], x[48], x[50],
        x[52], x[54], x[56], x[58], x[60], x[62],
    ]);
    let odd: [Complex<f32>; 32] = butterfly32([
        x[1], x[3], x[5], x[7], x[9], x[11], x[13], x[15], x[17], x[19], x[21], x[23], x[25],
        x[27], x[29], x[31], x[33], x[35], x[37], x[39], x[41], x[43], x[45], x[47], x[49], x[51],
        x[53], x[55], x[57], x[59], x[61], x[63],
    ]);

    let t: [Complex<f32>; 32] = [
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 0.0 / n as f32) * odd[0],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 1.0 / n as f32) * odd[1],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 2.0 / n as f32) * odd[2],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 3.0 / n as f32) * odd[3],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 4.0 / n as f32) * odd[4],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 5.0 / n as f32) * odd[5],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 6.0 / n as f32) * odd[6],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 7.0 / n as f32) * odd[7],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 8.0 / n as f32) * odd[8],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 9.0 / n as f32) * odd[9],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 10.0 / n as f32) * odd[10],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 11.0 / n as f32) * odd[11],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 12.0 / n as f32) * odd[12],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 13.0 / n as f32) * odd[13],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 14.0 / n as f32) * odd[14],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 15.0 / n as f32) * odd[15],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 16.0 / n as f32) * odd[16],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 17.0 / n as f32) * odd[17],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 18.0 / n as f32) * odd[18],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 19.0 / n as f32) * odd[19],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 20.0 / n as f32) * odd[20],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 21.0 / n as f32) * odd[21],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 22.0 / n as f32) * odd[22],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 23.0 / n as f32) * odd[23],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 24.0 / n as f32) * odd[24],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 25.0 / n as f32) * odd[25],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 26.0 / n as f32) * odd[26],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 27.0 / n as f32) * odd[27],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 28.0 / n as f32) * odd[28],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 29.0 / n as f32) * odd[29],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 30.0 / n as f32) * odd[30],
        Complex::exp(-2.0 * Complex::<f32>::i() * PI * 31.0 / n as f32) * odd[31],
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
        even[8] + t[8],
        even[9] + t[9],
        even[10] + t[10],
        even[11] + t[11],
        even[12] + t[12],
        even[13] + t[13],
        even[14] + t[14],
        even[15] + t[15],
        even[16] + t[16],
        even[17] + t[17],
        even[18] + t[18],
        even[19] + t[19],
        even[20] + t[20],
        even[21] + t[21],
        even[22] + t[22],
        even[23] + t[23],
        even[24] + t[24],
        even[25] + t[25],
        even[26] + t[26],
        even[27] + t[27],
        even[28] + t[28],
        even[29] + t[29],
        even[30] + t[30],
        even[31] + t[31],
        even[0] - t[0],
        even[1] - t[1],
        even[2] - t[2],
        even[3] - t[3],
        even[4] - t[4],
        even[5] - t[5],
        even[6] - t[6],
        even[7] - t[7],
        even[8] - t[8],
        even[9] - t[9],
        even[10] - t[10],
        even[11] - t[11],
        even[12] - t[12],
        even[13] - t[13],
        even[14] - t[14],
        even[15] - t[15],
        even[0] - t[0],
        even[1] - t[1],
        even[2] - t[2],
        even[3] - t[3],
        even[4] - t[4],
        even[5] - t[5],
        even[6] - t[6],
        even[7] - t[7],
        even[8] - t[8],
        even[9] - t[9],
        even[10] - t[10],
        even[11] - t[11],
        even[12] - t[12],
        even[13] - t[13],
        even[14] - t[14],
        even[15] - t[15],
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
