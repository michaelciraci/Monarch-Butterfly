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
        fft2([Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)]).to_vec(),
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

    let monarch = fft3(&buf);

    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_butterfly_4() {
    assert_eq!(
        fft4([
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
fn test_fft5() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(5);
    let mut buf: Vec<_> = (0..5).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

    let monarch = fft5(&buf);

    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft6() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(6);
    let mut buf: Vec<_> = (0..6).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

    let monarch = fft6(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft7() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(7);
    let mut buf: Vec<_> = (0..7).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

    let monarch = fft7(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_8() {
    let a = butterfly((1..9).map(|i| Complex::new(i as f32, 0.0)).collect());
    assert_eq!(a[0], Complex::new(36.0, 0.0));
    assert_eq!(a[6], Complex::new(-4.0, -4.0));
}

#[test]
fn test_fft9() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(9);
    let mut buf: Vec<_> = (0..9).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

    let monarch = fft9(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft10() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(10);
    let mut buf: Vec<_> = (0..10)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft10(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft11() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(11);
    let mut buf: Vec<_> = (0..11)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft11(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft12() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(12);
    let mut buf: Vec<_> = (0..12)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft12(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft13() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(13);
    let mut buf: Vec<_> = (0..13)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft13(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft14() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(14);
    let mut buf: Vec<_> = (0..14)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft14(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft15() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(15);
    let mut buf: Vec<_> = (0..15)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft15(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft17() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(17);
    let mut buf: Vec<_> = (0..17)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft17(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft18() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(18);
    let mut buf: Vec<_> = (0..18)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft18(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft19() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(19);
    let mut buf: Vec<_> = (0..19)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft19(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft20() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(20);
    let mut buf: Vec<_> = (0..20)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft20(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft21() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(21);
    let mut buf: Vec<_> = (0..21)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft21(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft22() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(22);
    let mut buf: Vec<_> = (0..22)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft22(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft23() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(23);
    let mut buf: Vec<_> = (0..23)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft23(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft24() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(24);
    let mut buf: Vec<_> = (0..24)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft24(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft25() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(25);
    let mut buf: Vec<_> = (0..25)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft25(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft26() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(26);
    let mut buf: Vec<_> = (0..26)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft26(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft27() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(27);
    let mut buf: Vec<_> = (0..27)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft27(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft28() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(28);
    let mut buf: Vec<_> = (0..28)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft28(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft29() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(29);
    let mut buf: Vec<_> = (0..29)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft29(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft30() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(30);
    let mut buf: Vec<_> = (0..30)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft30(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft31() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(31);
    let mut buf: Vec<_> = (0..31)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft31(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft33() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(33);
    let mut buf: Vec<_> = (0..33)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft33(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft34() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(34);
    let mut buf: Vec<_> = (0..34)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft34(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft35() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(35);
    let mut buf: Vec<_> = (0..35)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft35(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft36() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(36);
    let mut buf: Vec<_> = (0..36)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft36(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft37() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(37);
    let mut buf: Vec<_> = (0..37)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft37(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft38() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(38);
    let mut buf: Vec<_> = (0..38)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft38(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft39() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(39);
    let mut buf: Vec<_> = (0..39)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft39(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft40() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(40);
    let mut buf: Vec<_> = (0..40)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft40(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft41() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(41);
    let mut buf: Vec<_> = (0..41)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft41(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft42() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(42);
    let mut buf: Vec<_> = (0..42)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft42(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft43() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(43);
    let mut buf: Vec<_> = (0..43)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft43(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft44() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(44);
    let mut buf: Vec<_> = (0..44)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft44(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft45() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(45);
    let mut buf: Vec<_> = (0..45)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft45(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft46() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(46);
    let mut buf: Vec<_> = (0..46)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft46(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft47() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(47);
    let mut buf: Vec<_> = (0..47)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft47(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft48() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(48);
    let mut buf: Vec<_> = (0..48)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft48(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft49() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(49);
    let mut buf: Vec<_> = (0..49)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft49(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft50() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(50);
    let mut buf: Vec<_> = (0..50)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft50(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft51() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(51);
    let mut buf: Vec<_> = (0..51)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft51(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft52() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(52);
    let mut buf: Vec<_> = (0..52)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft52(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft53() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(53);
    let mut buf: Vec<_> = (0..53)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft53(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft54() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(54);
    let mut buf: Vec<_> = (0..54)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft54(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft55() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(55);
    let mut buf: Vec<_> = (0..55)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft55(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft56() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(56);
    let mut buf: Vec<_> = (0..56)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft56(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft57() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(57);
    let mut buf: Vec<_> = (0..57)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft57(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft58() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(58);
    let mut buf: Vec<_> = (0..58)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft58(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
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

    let ba = fft1024(a);
    let bv = butterfly(v);

    assert_slice_equal!(ba, bv, 0.001);
}
