use monarch_butterfly::*;
use num_complex::Complex;

macro_rules! assert_slice_equal {
    ($lhs:expr,$rhs:expr) => {
        assert_slice_equal!($lhs, $rhs, 0.0001);
    };
    ($lhs:expr,$rhs:expr,$tol:expr) => {
        assert_eq!($lhs.len(), $rhs.len());
        for idx in 0..($lhs.len()) {
            assert!(
                ($lhs[idx].re - $rhs[idx].re).abs() < $tol,
                "index {} does not match: ({} - {}).abs() < {}",
                idx,
                $lhs[idx].re,
                $rhs[idx].re,
                $tol
            );
            assert!(
                ($lhs[idx].im - $rhs[idx].im).abs() < $tol,
                "index {} does not match: ({} - {}).abs() < {}",
                idx,
                $lhs[idx].im,
                $rhs[idx].im,
                $tol
            );
        }
    };
}

macro_rules! compare_against_rustfft {
    ($test:ident, $n:expr) => {
        #[test]
        fn $test() {
            let mut p = rustfft::FftPlanner::new();
            let plan = p.plan_fft_forward($n);
            let mut buf: Vec<_> = (0..$n)
                .map(|i| Complex::<f64>::new(i as f64 + 7.0, (i * i % 31) as f64))
                .collect();
            let monarch = fft::<$n, _, _>(&buf);
            plan.process(&mut buf);
            assert_slice_equal!(monarch, buf);
        }
    };
}

compare_against_rustfft!(forward_1, 1);
compare_against_rustfft!(forward_2, 2);
compare_against_rustfft!(forward_3, 3);
compare_against_rustfft!(forward_4, 4);
compare_against_rustfft!(forward_5, 5);
compare_against_rustfft!(forward_6, 6);
compare_against_rustfft!(forward_7, 7);
compare_against_rustfft!(forward_8, 8);
compare_against_rustfft!(forward_9, 9);
compare_against_rustfft!(forward_10, 10);
compare_against_rustfft!(forward_11, 11);
compare_against_rustfft!(forward_12, 12);
compare_against_rustfft!(forward_13, 13);
compare_against_rustfft!(forward_14, 14);
compare_against_rustfft!(forward_15, 15);
compare_against_rustfft!(forward_16, 16);
compare_against_rustfft!(forward_17, 17);
compare_against_rustfft!(forward_18, 18);
compare_against_rustfft!(forward_19, 19);
compare_against_rustfft!(forward_20, 20);
compare_against_rustfft!(forward_21, 21);
compare_against_rustfft!(forward_22, 22);
compare_against_rustfft!(forward_23, 23);
compare_against_rustfft!(forward_24, 24);
compare_against_rustfft!(forward_25, 25);
compare_against_rustfft!(forward_26, 26);
compare_against_rustfft!(forward_27, 27);
compare_against_rustfft!(forward_28, 28);
compare_against_rustfft!(forward_29, 29);
compare_against_rustfft!(forward_30, 30);
compare_against_rustfft!(forward_31, 31);
compare_against_rustfft!(forward_32, 32);
compare_against_rustfft!(forward_33, 33);
compare_against_rustfft!(forward_34, 34);
compare_against_rustfft!(forward_35, 35);
compare_against_rustfft!(forward_36, 36);
compare_against_rustfft!(forward_37, 37);
compare_against_rustfft!(forward_38, 38);
compare_against_rustfft!(forward_39, 39);
compare_against_rustfft!(forward_40, 40);
compare_against_rustfft!(forward_41, 41);
compare_against_rustfft!(forward_42, 42);
compare_against_rustfft!(forward_43, 43);
compare_against_rustfft!(forward_44, 44);
compare_against_rustfft!(forward_45, 45);
compare_against_rustfft!(forward_46, 46);
compare_against_rustfft!(forward_47, 47);
compare_against_rustfft!(forward_48, 48);
compare_against_rustfft!(forward_49, 49);
compare_against_rustfft!(forward_50, 50);
compare_against_rustfft!(forward_51, 51);
compare_against_rustfft!(forward_52, 52);
compare_against_rustfft!(forward_53, 53);
compare_against_rustfft!(forward_54, 54);
compare_against_rustfft!(forward_55, 55);
compare_against_rustfft!(forward_56, 56);
compare_against_rustfft!(forward_57, 57);
compare_against_rustfft!(forward_58, 58);
compare_against_rustfft!(forward_59, 59);
compare_against_rustfft!(forward_60, 60);
compare_against_rustfft!(forward_61, 61);
compare_against_rustfft!(forward_62, 62);
compare_against_rustfft!(forward_63, 63);
compare_against_rustfft!(forward_64, 64);
compare_against_rustfft!(forward_65, 65);
compare_against_rustfft!(forward_66, 66);
compare_against_rustfft!(forward_67, 67);
compare_against_rustfft!(forward_68, 68);
compare_against_rustfft!(forward_69, 69);
compare_against_rustfft!(forward_70, 70);
compare_against_rustfft!(forward_71, 71);
compare_against_rustfft!(forward_72, 72);
compare_against_rustfft!(forward_73, 73);
compare_against_rustfft!(forward_74, 74);
compare_against_rustfft!(forward_75, 75);
compare_against_rustfft!(forward_76, 76);
compare_against_rustfft!(forward_77, 77);
compare_against_rustfft!(forward_78, 78);
compare_against_rustfft!(forward_79, 79);
compare_against_rustfft!(forward_80, 80);
compare_against_rustfft!(forward_81, 81);
compare_against_rustfft!(forward_82, 82);
compare_against_rustfft!(forward_83, 83);
compare_against_rustfft!(forward_84, 84);
compare_against_rustfft!(forward_85, 85);
compare_against_rustfft!(forward_86, 86);
compare_against_rustfft!(forward_87, 87);
compare_against_rustfft!(forward_88, 88);
compare_against_rustfft!(forward_89, 89);
compare_against_rustfft!(forward_90, 90);
compare_against_rustfft!(forward_91, 91);
compare_against_rustfft!(forward_92, 92);
compare_against_rustfft!(forward_93, 93);
compare_against_rustfft!(forward_94, 94);
compare_against_rustfft!(forward_95, 95);
compare_against_rustfft!(forward_96, 96);
compare_against_rustfft!(forward_97, 97);
compare_against_rustfft!(forward_98, 98);
compare_against_rustfft!(forward_99, 99);
compare_against_rustfft!(forward_100, 100);
compare_against_rustfft!(forward_101, 101);
compare_against_rustfft!(forward_102, 102);
compare_against_rustfft!(forward_103, 103);
compare_against_rustfft!(forward_104, 104);
compare_against_rustfft!(forward_105, 105);
compare_against_rustfft!(forward_106, 106);
compare_against_rustfft!(forward_107, 107);
compare_against_rustfft!(forward_108, 108);
compare_against_rustfft!(forward_109, 109);
compare_against_rustfft!(forward_110, 110);
compare_against_rustfft!(forward_111, 111);
compare_against_rustfft!(forward_112, 112);
compare_against_rustfft!(forward_113, 113);
compare_against_rustfft!(forward_114, 114);
compare_against_rustfft!(forward_115, 115);
compare_against_rustfft!(forward_116, 116);
compare_against_rustfft!(forward_117, 117);
compare_against_rustfft!(forward_118, 118);
compare_against_rustfft!(forward_119, 119);
compare_against_rustfft!(forward_120, 120);
compare_against_rustfft!(forward_121, 121);
compare_against_rustfft!(forward_122, 122);
compare_against_rustfft!(forward_123, 123);
compare_against_rustfft!(forward_124, 124);
compare_against_rustfft!(forward_125, 125);
compare_against_rustfft!(forward_126, 126);
compare_against_rustfft!(forward_127, 127);
compare_against_rustfft!(forward_128, 128);
compare_against_rustfft!(forward_129, 129);
compare_against_rustfft!(forward_130, 130);
compare_against_rustfft!(forward_131, 131);
compare_against_rustfft!(forward_132, 132);
compare_against_rustfft!(forward_133, 133);
compare_against_rustfft!(forward_134, 134);
compare_against_rustfft!(forward_135, 135);
compare_against_rustfft!(forward_136, 136);
compare_against_rustfft!(forward_137, 137);
compare_against_rustfft!(forward_138, 138);
compare_against_rustfft!(forward_139, 139);
compare_against_rustfft!(forward_140, 140);
compare_against_rustfft!(forward_141, 141);
compare_against_rustfft!(forward_142, 142);
compare_against_rustfft!(forward_143, 143);
compare_against_rustfft!(forward_144, 144);
compare_against_rustfft!(forward_145, 145);
compare_against_rustfft!(forward_146, 146);
compare_against_rustfft!(forward_147, 147);
compare_against_rustfft!(forward_148, 148);
compare_against_rustfft!(forward_149, 149);
compare_against_rustfft!(forward_150, 150);
compare_against_rustfft!(forward_151, 151);
compare_against_rustfft!(forward_152, 152);
compare_against_rustfft!(forward_153, 153);
compare_against_rustfft!(forward_154, 154);
compare_against_rustfft!(forward_155, 155);
compare_against_rustfft!(forward_156, 156);
compare_against_rustfft!(forward_157, 157);
compare_against_rustfft!(forward_158, 158);
compare_against_rustfft!(forward_159, 159);
compare_against_rustfft!(forward_160, 160);
compare_against_rustfft!(forward_161, 161);
compare_against_rustfft!(forward_162, 162);
compare_against_rustfft!(forward_163, 163);
compare_against_rustfft!(forward_164, 164);
compare_against_rustfft!(forward_165, 165);
compare_against_rustfft!(forward_166, 166);
compare_against_rustfft!(forward_167, 167);
compare_against_rustfft!(forward_168, 168);
compare_against_rustfft!(forward_169, 169);
compare_against_rustfft!(forward_170, 170);
compare_against_rustfft!(forward_171, 171);
compare_against_rustfft!(forward_172, 172);
compare_against_rustfft!(forward_173, 173);
compare_against_rustfft!(forward_174, 174);
compare_against_rustfft!(forward_175, 175);
compare_against_rustfft!(forward_176, 176);
compare_against_rustfft!(forward_177, 177);
compare_against_rustfft!(forward_178, 178);
compare_against_rustfft!(forward_179, 179);
compare_against_rustfft!(forward_180, 180);
compare_against_rustfft!(forward_181, 181);
compare_against_rustfft!(forward_182, 182);
compare_against_rustfft!(forward_183, 183);
compare_against_rustfft!(forward_184, 184);
compare_against_rustfft!(forward_185, 185);
compare_against_rustfft!(forward_186, 186);
compare_against_rustfft!(forward_187, 187);
compare_against_rustfft!(forward_188, 188);
compare_against_rustfft!(forward_189, 189);
compare_against_rustfft!(forward_190, 190);
compare_against_rustfft!(forward_191, 191);
compare_against_rustfft!(forward_192, 192);
compare_against_rustfft!(forward_193, 193);
compare_against_rustfft!(forward_194, 194);
compare_against_rustfft!(forward_195, 195);
compare_against_rustfft!(forward_196, 196);
compare_against_rustfft!(forward_197, 197);
compare_against_rustfft!(forward_198, 198);
compare_against_rustfft!(forward_199, 199);
compare_against_rustfft!(forward_200, 200);

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
fn test_butterfly_2() {
    assert_eq!(
        fft::<2, _, _>([Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)]).to_vec(),
        butterfly(vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)])
    );
}

#[test]
fn test_fft3() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(3);
    let mut buf = vec![
        Complex::<f64>::new(0.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(2.0, 0.0),
    ];

    let monarch = fft::<3, _, _>(&buf);

    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn round_trip_complex_f32() {
    let input: [Complex<f32>; 15] =
        core::array::from_fn(|i| Complex::new((i % 7) as f32, (i % 3) as f32));
    let spectrum = fft::<15, _, _>(&input[..]);
    let output = ifft::<15, _, _>(spectrum);
    for (actual, sample) in output.iter().zip(input) {
        assert!((*actual - sample * 15.0).norm() < 1e-4);
    }
}

#[test]
fn power_of_two_128_matches_rustfft_in_both_directions() {
    let input: [Complex<f32>; 128] = core::array::from_fn(|i| {
        Complex::new(
            (i * 11 % 31) as f32 / 16.0 - 1.0,
            (i * i % 23) as f32 / 12.0 - 1.0,
        )
    });
    let mut planner = rustfft::FftPlanner::new();
    for (actual, plan) in [
        (fft::<128, _, _>(input), planner.plan_fft_forward(128)),
        (ifft::<128, _, _>(input), planner.plan_fft_inverse(128)),
    ] {
        let mut expected = input;
        plan.process(&mut expected);
        assert_slice_equal!(actual, expected);
    }
}

#[test]
fn optimized_paths_match_rustfft_in_both_directions() {
    fn check<const N: usize, T: num_traits::Float + num_traits::FloatConst + rustfft::FftNum>(
        tolerance: T,
    ) {
        let dense: [Complex<T>; N] = core::array::from_fn(|i| {
            let re = (i * 13 % 31) as f64 / 31.0 - 0.5;
            let im = (i * i % 17) as f64 / 17.0 - 0.5;
            Complex::new(T::from_f64(re).unwrap(), T::from_f64(im).unwrap())
        });
        let mut impulse = [Complex::new(T::zero(), T::zero()); N];
        impulse[N / 4 - 1] = Complex::new(T::from_f64(0.75).unwrap(), T::from_f64(-0.25).unwrap());

        let mut planner = rustfft::FftPlanner::new();
        let forward = planner.plan_fft_forward(N);
        let inverse = planner.plan_fft_inverse(N);
        for input in [dense, impulse] {
            for (actual, plan) in [
                (fft::<N, _, _>(input), &forward),
                (ifft::<N, _, _>(&input[..]), &inverse),
            ] {
                let mut expected = input;
                plan.process(&mut expected);
                for (bin, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
                    assert!(
                        (actual - expected).norm() < tolerance,
                        "{N}-point {:?} bin {bin}: {actual:?} differs from {expected:?}",
                        plan.fft_direction()
                    );
                }
            }
        }
    }

    check::<8, _>(1e-5_f32);
    check::<16, _>(1e-5_f32);
    check::<18, _>(1e-5_f32);
    check::<24, _>(1e-5_f32);
    check::<27, _>(1e-5_f32);
    check::<32, _>(1e-5_f32);
    check::<48, _>(1e-5_f32);
    check::<48, _>(1e-11_f64);
    check::<64, _>(1e-5_f32);
    check::<96, _>(1e-5_f32);
    check::<96, _>(1e-11_f64);
    check::<162, _>(2e-5_f32);
    check::<192, _>(2e-5_f32);
    check::<193, _>(2e-5_f32);
    check::<198, _>(2e-5_f32);

    check::<29, _>(1e-5_f32);
    check::<97, _>(2e-5_f32);
    check::<179, _>(2e-5_f32);
    check::<199, _>(2e-5_f32);
    check::<199, _>(1e-11_f64);

    check::<47, _>(1e-5_f32);
    check::<47, _>(1e-11_f64);
    check::<94, _>(2e-5_f32);
    check::<167, _>(2e-5_f32);
    check::<167, _>(1e-11_f64);
}

#[test]
fn prime_13_complex_impulse() {
    let mut input = [Complex::new(0.0_f32, 0.0); 13];
    input[5] = Complex::new(0.5, -1.0);
    for (spectrum, sign) in [
        (fft::<13, _, _>(input), -1.0),
        (ifft::<13, _, _>(input), 1.0),
    ] {
        for (k, bin) in spectrum.iter().enumerate() {
            let angle = sign * core::f64::consts::TAU * (5 * k) as f64 / 13.0;
            let expected = Complex::new(0.5, -1.0) * Complex::from_polar(1.0, angle);
            assert!((Complex::new(bin.re as f64, bin.im as f64) - expected).norm() < 2e-7);
        }
    }
}

#[test]
fn test_butterfly_4() {
    assert_eq!(
        fft::<4, _, _>([
            Complex::new(1.0_f32, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0)
        ]),
        [
            Complex::new(10.0, 0.0),
            Complex::new(-2.0, 2.0),
            Complex::new(-2.0, 0.0),
            Complex::new(-2.0, -2.0),
        ]
    );
}
