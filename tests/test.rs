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
fn test_fft59() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(59);
    let mut buf: Vec<_> = (0..59)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft59(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft60() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(60);
    let mut buf: Vec<_> = (0..60)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft60(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft61() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(61);
    let mut buf: Vec<_> = (0..61)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft61(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft62() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(62);
    let mut buf: Vec<_> = (0..62)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft62(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft63() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(63);
    let mut buf: Vec<_> = (0..63)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft63(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft64() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(64);
    let mut buf: Vec<_> = (0..64)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft64(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft65() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(65);
    let mut buf: Vec<_> = (0..65)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft65(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft66() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(66);
    let mut buf: Vec<_> = (0..66)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft66(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft67() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(67);
    let mut buf: Vec<_> = (0..67)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft67(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft68() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(68);
    let mut buf: Vec<_> = (0..68)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft68(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft69() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(69);
    let mut buf: Vec<_> = (0..69)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft69(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft70() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(70);
    let mut buf: Vec<_> = (0..70)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft70(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft71() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(71);
    let mut buf: Vec<_> = (0..71)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft71(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft72() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(72);
    let mut buf: Vec<_> = (0..72)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft72(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft73() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(73);
    let mut buf: Vec<_> = (0..73)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft73(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft74() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(74);
    let mut buf: Vec<_> = (0..74)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft74(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft75() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(75);
    let mut buf: Vec<_> = (0..75)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft75(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft76() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(76);
    let mut buf: Vec<_> = (0..76)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft76(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft77() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(77);
    let mut buf: Vec<_> = (0..77)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft77(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft78() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(78);
    let mut buf: Vec<_> = (0..78)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft78(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft79() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(79);
    let mut buf: Vec<_> = (0..79)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft79(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft80() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(80);
    let mut buf: Vec<_> = (0..80)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft80(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft81() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(81);
    let mut buf: Vec<_> = (0..81)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft81(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft82() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(82);
    let mut buf: Vec<_> = (0..82)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft82(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft83() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(83);
    let mut buf: Vec<_> = (0..83)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft83(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft84() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(84);
    let mut buf: Vec<_> = (0..84)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft84(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft85() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(85);
    let mut buf: Vec<_> = (0..85)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft85(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft86() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(86);
    let mut buf: Vec<_> = (0..86)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft86(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft87() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(87);
    let mut buf: Vec<_> = (0..87)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft87(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft88() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(88);
    let mut buf: Vec<_> = (0..88)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft88(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft89() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(89);
    let mut buf: Vec<_> = (0..89)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft89(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft90() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(90);
    let mut buf: Vec<_> = (0..90)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft90(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft91() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(91);
    let mut buf: Vec<_> = (0..91)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft91(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft92() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(92);
    let mut buf: Vec<_> = (0..92)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft92(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft93() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(93);
    let mut buf: Vec<_> = (0..93)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft93(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft94() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(94);
    let mut buf: Vec<_> = (0..94)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft94(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft95() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(95);
    let mut buf: Vec<_> = (0..95)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft95(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft96() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(96);
    let mut buf: Vec<_> = (0..96)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft96(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft97() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(97);
    let mut buf: Vec<_> = (0..97)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft97(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft98() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(98);
    let mut buf: Vec<_> = (0..98)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft98(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft99() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(99);
    let mut buf: Vec<_> = (0..99)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft99(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft100() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(100);
    let mut buf: Vec<_> = (0..100)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft100(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft101() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(101);
    let mut buf: Vec<_> = (0..101)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft101(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft102() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(102);
    let mut buf: Vec<_> = (0..102)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft102(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft103() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(103);
    let mut buf: Vec<_> = (0..103)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft103(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft104() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(104);
    let mut buf: Vec<_> = (0..104)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft104(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft105() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(105);
    let mut buf: Vec<_> = (0..105)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft105(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft106() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(106);
    let mut buf: Vec<_> = (0..106)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft106(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft107() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(107);
    let mut buf: Vec<_> = (0..107)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft107(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft108() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(108);
    let mut buf: Vec<_> = (0..108)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft108(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft109() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(109);
    let mut buf: Vec<_> = (0..109)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft109(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft110() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(110);
    let mut buf: Vec<_> = (0..110)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft110(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft111() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(111);
    let mut buf: Vec<_> = (0..111)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft111(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft112() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(112);
    let mut buf: Vec<_> = (0..112)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft112(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft113() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(113);
    let mut buf: Vec<_> = (0..113)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft113(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft114() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(114);
    let mut buf: Vec<_> = (0..114)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft114(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft115() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(115);
    let mut buf: Vec<_> = (0..115)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft115(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft116() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(116);
    let mut buf: Vec<_> = (0..116)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft116(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft117() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(117);
    let mut buf: Vec<_> = (0..117)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft117(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft118() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(118);
    let mut buf: Vec<_> = (0..118)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft118(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft119() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(119);
    let mut buf: Vec<_> = (0..119)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft119(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft120() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(120);
    let mut buf: Vec<_> = (0..120)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft120(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft121() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(121);
    let mut buf: Vec<_> = (0..121)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft121(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft122() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(122);
    let mut buf: Vec<_> = (0..122)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft122(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft123() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(123);
    let mut buf: Vec<_> = (0..123)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft123(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft124() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(124);
    let mut buf: Vec<_> = (0..124)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft124(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

// #[test]
// fn test_fft125() {
//     let mut p = rustfft::FftPlanner::new();
//     let plan = p.plan_fft_forward(125);
//     let mut buf: Vec<_> = (0..125)
//         .map(|i| Complex::<f64>::new(i as f64, 0.0))
//         .collect();

//     let monarch = fft125(&buf);
//     plan.process(&mut buf);

//     assert_slice_equal!(monarch, buf);
// }

#[test]
fn test_fft126() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(126);
    let mut buf: Vec<_> = (0..126)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft126(&buf);
    plan.process(&mut buf);

    assert_slice_equal!(monarch, buf);
}

#[test]
fn test_fft127() {
    let mut p = rustfft::FftPlanner::new();
    let plan = p.plan_fft_forward(127);
    let mut buf: Vec<_> = (0..127)
        .map(|i| Complex::<f64>::new(i as f64, 0.0))
        .collect();

    let monarch = fft127(&buf);
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
