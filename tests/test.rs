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
                .map(|i| Complex::<f64>::new(i as f64, 0.0))
                .collect();
            let monarch = fft::<$n, _, _>(&buf);
            plan.process(&mut buf);
            assert_slice_equal!(monarch, buf);
        }
    };
}

compare_against_rustfft!(test_1, 1);
compare_against_rustfft!(test_2, 2);
compare_against_rustfft!(test_3, 3);
compare_against_rustfft!(test_4, 4);
compare_against_rustfft!(test_5, 5);
compare_against_rustfft!(test_6, 6);
compare_against_rustfft!(test_7, 7);
compare_against_rustfft!(test_8, 8);
compare_against_rustfft!(test_9, 9);
compare_against_rustfft!(test_10, 10);
compare_against_rustfft!(test_11, 11);
compare_against_rustfft!(test_12, 12);
compare_against_rustfft!(test_13, 13);
compare_against_rustfft!(test_14, 14);
compare_against_rustfft!(test_15, 15);
compare_against_rustfft!(test_16, 16);
compare_against_rustfft!(test_17, 17);
compare_against_rustfft!(test_18, 18);
compare_against_rustfft!(test_19, 19);
compare_against_rustfft!(test_20, 20);
compare_against_rustfft!(test_21, 21);
compare_against_rustfft!(test_22, 22);
compare_against_rustfft!(test_23, 23);
compare_against_rustfft!(test_24, 24);
compare_against_rustfft!(test_25, 25);
compare_against_rustfft!(test_26, 26);
compare_against_rustfft!(test_27, 27);
compare_against_rustfft!(test_28, 28);
compare_against_rustfft!(test_29, 29);
compare_against_rustfft!(test_30, 30);
compare_against_rustfft!(test_31, 31);
compare_against_rustfft!(test_32, 32);
compare_against_rustfft!(test_33, 33);
compare_against_rustfft!(test_34, 34);
compare_against_rustfft!(test_35, 35);
compare_against_rustfft!(test_36, 36);
compare_against_rustfft!(test_37, 37);
compare_against_rustfft!(test_38, 38);
compare_against_rustfft!(test_39, 39);
compare_against_rustfft!(test_40, 40);
compare_against_rustfft!(test_41, 41);
compare_against_rustfft!(test_42, 42);
compare_against_rustfft!(test_43, 43);
compare_against_rustfft!(test_44, 44);
compare_against_rustfft!(test_45, 45);
compare_against_rustfft!(test_46, 46);
compare_against_rustfft!(test_47, 47);
compare_against_rustfft!(test_48, 48);
compare_against_rustfft!(test_49, 49);
compare_against_rustfft!(test_50, 50);
compare_against_rustfft!(test_51, 51);
compare_against_rustfft!(test_52, 52);
compare_against_rustfft!(test_53, 53);
compare_against_rustfft!(test_54, 54);
compare_against_rustfft!(test_55, 55);
compare_against_rustfft!(test_56, 56);
compare_against_rustfft!(test_57, 57);
compare_against_rustfft!(test_58, 58);
compare_against_rustfft!(test_59, 59);
compare_against_rustfft!(test_60, 60);
compare_against_rustfft!(test_61, 61);
compare_against_rustfft!(test_62, 62);
compare_against_rustfft!(test_63, 63);
compare_against_rustfft!(test_64, 64);
compare_against_rustfft!(test_65, 65);
compare_against_rustfft!(test_66, 66);
compare_against_rustfft!(test_67, 67);
compare_against_rustfft!(test_68, 68);
compare_against_rustfft!(test_69, 69);
compare_against_rustfft!(test_70, 70);
compare_against_rustfft!(test_71, 71);
compare_against_rustfft!(test_72, 72);
compare_against_rustfft!(test_73, 73);
compare_against_rustfft!(test_74, 74);
compare_against_rustfft!(test_75, 75);
compare_against_rustfft!(test_76, 76);
compare_against_rustfft!(test_77, 77);
compare_against_rustfft!(test_78, 78);
compare_against_rustfft!(test_79, 79);
compare_against_rustfft!(test_80, 80);
compare_against_rustfft!(test_81, 81);
compare_against_rustfft!(test_82, 82);
compare_against_rustfft!(test_83, 83);
compare_against_rustfft!(test_84, 84);
compare_against_rustfft!(test_85, 85);
compare_against_rustfft!(test_86, 86);
compare_against_rustfft!(test_87, 87);
compare_against_rustfft!(test_88, 88);
compare_against_rustfft!(test_89, 89);
compare_against_rustfft!(test_90, 90);
compare_against_rustfft!(test_91, 91);
compare_against_rustfft!(test_92, 92);
compare_against_rustfft!(test_93, 93);
compare_against_rustfft!(test_94, 94);
compare_against_rustfft!(test_95, 95);
compare_against_rustfft!(test_96, 96);
compare_against_rustfft!(test_97, 97);
compare_against_rustfft!(test_98, 98);
compare_against_rustfft!(test_99, 99);
compare_against_rustfft!(test_100, 100);
compare_against_rustfft!(test_101, 101);
compare_against_rustfft!(test_102, 102);
compare_against_rustfft!(test_103, 103);
compare_against_rustfft!(test_104, 104);
compare_against_rustfft!(test_105, 105);

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
fn test_butterfly_4() {
    assert_eq!(
        fft::<4, _, _>([
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
