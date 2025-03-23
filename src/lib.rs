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


fn mixed100(x: [Complex<f32>; 100]) -> [Complex<f32>; 100] {
    let mut cross_fft_len = 5;
    let mut twiddle_factors = vec![];

    for factor in [5, 4] {
        let cross_fft_columns = cross_fft_len;
        cross_fft_len *= factor;

        for i in 0..cross_fft_columns {
            for k in 1..factor {
                let twiddle = compute_twiddle(i * k, cross_fft_len);
                twiddle_factors.push(twiddle);
            }
        }
    }

    let mut base = {
        let b0 = butterfly5([x[0], x[1], x[2], x[3], x[4]]);
        let b1 = butterfly5([x[5], x[6], x[7], x[8], x[9]]);
        let b2 = butterfly5([x[10], x[11], x[12], x[13], x[14]]);
        let b3 = butterfly5([x[15], x[16], x[17], x[18], x[19]]);
        let b4 = butterfly5([x[20], x[21], x[22], x[23], x[24]]);
        let b5 = butterfly5([x[25], x[26], x[27], x[28], x[19]]);
        let b6 = butterfly5([x[30], x[31], x[32], x[33], x[34]]);
        let b7 = butterfly5([x[35], x[36], x[37], x[38], x[39]]);
        let b8 = butterfly5([x[40], x[41], x[42], x[43], x[44]]);
        let b9 = butterfly5([x[45], x[46], x[47], x[48], x[49]]);
        let b10 = butterfly5([x[50], x[51], x[52], x[53], x[54]]);
        let b11 = butterfly5([x[55], x[56], x[57], x[58], x[59]]);
        let b12 = butterfly5([x[60], x[61], x[62], x[63], x[64]]);
        let b13 = butterfly5([x[65], x[66], x[67], x[68], x[69]]);
        let b14 = butterfly5([x[70], x[71], x[72], x[73], x[74]]);
        let b15 = butterfly5([x[75], x[76], x[77], x[78], x[79]]);
        let b16 = butterfly5([x[80], x[81], x[82], x[83], x[84]]);
        let b17 = butterfly5([x[85], x[86], x[87], x[88], x[89]]);
        let b18 = butterfly5([x[90], x[91], x[92], x[93], x[94]]);
        let b19 = butterfly5([x[95], x[96], x[97], x[98], x[99]]);
        [
            b0[0], b0[1], b0[2], b0[3], b0[4], b1[0], b1[1], b1[2], b1[3], b1[4], b2[0], b2[1],
            b2[2], b2[3], b2[4], b3[0], b3[1], b3[2], b3[3], b3[4], b4[0], b4[1], b4[2], b4[3],
            b4[4], b5[0], b5[1], b5[2], b5[3], b5[4], b6[0], b6[1], b6[2], b6[3], b6[4], b7[0],
            b7[1], b7[2], b7[3], b7[4], b8[0], b8[1], b8[2], b8[3], b8[4], b9[0], b9[1], b9[2],
            b9[3], b9[4], b10[0], b10[1], b10[2], b10[3], b10[4], b11[0], b11[1], b11[2], b11[3],
            b11[4], b12[0], b12[1], b12[2], b12[3], b12[4], b13[0], b13[1], b13[2], b13[3], b13[4],
            b14[0], b14[1], b14[2], b14[3], b14[4], b15[0], b15[1], b15[2], b15[3], b15[4], b16[0],
            b16[1], b16[2], b16[3], b16[4], b17[0], b17[1], b17[2], b17[3], b17[4], b18[0], b18[1],
            b18[2], b18[3], b18[4], b19[0], b19[1], b19[2], b19[3], b19[4],
        ]
    };

    let mut cross_fft_len = 5;
    let mut layer_twiddles = twiddle_factors;

    for factor in [5, 4] {
        let cross_fft_columns = cross_fft_len;
        cross_fft_len *= factor;

        // First iteration do FFTs of size 25
        let b0 = butterfly5([
            base[0],
            base[5] * layer_twiddles[0],
            base[10] * layer_twiddles[1],
            base[15] * layer_twiddles[2],
            base[20] * layer_twiddles[3],
        ]);
        let b1 = butterfly5([
            base[1],
            base[6] * layer_twiddles[4],
            base[11] * layer_twiddles[5],
            base[16] * layer_twiddles[6],
            base[21] * layer_twiddles[7],
        ]);
        let b2 = butterfly5([
            base[2],
            base[7] * layer_twiddles[8],
            base[12] * layer_twiddles[9],
            base[17] * layer_twiddles[10],
            base[22] * layer_twiddles[11],
        ]);
        let b3 = butterfly5([
            base[3],
            base[8] * layer_twiddles[12],
            base[13] * layer_twiddles[13],
            base[18] * layer_twiddles[14],
            base[23] * layer_twiddles[15],
        ]);
        let b4 = butterfly5([
            base[4],
            base[9] * layer_twiddles[16],
            base[14] * layer_twiddles[17],
            base[19] * layer_twiddles[18],
            base[24] * layer_twiddles[19],
        ]);
    }

    todo!()
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
