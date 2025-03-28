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
