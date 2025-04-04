#![doc = include_str!("../README.md")]

#![allow(clippy::excessive_precision)]
#![forbid(unsafe_code)]
#![no_std]

use num_complex::Complex;
use num_traits::{Float, FloatConst};

const SQRT_3: f64 = 1.7320508075688772;
const SQRT_3_DIV_2: f64 = SQRT_3 / 2.0;

monarch_derive::generate_switch!();
monarch_derive::generate_powers_of_two!();
monarch_derive::generate_coprimes!();
monarch_derive::generate_mixed_radix!();
monarch_derive::generate_primes!();
monarch_derive::generate_iffts!();

fn _compute_twiddle<T: Float + FloatConst>(index: usize, fft_len: usize) -> Complex<T> {
    let constant = T::from(-2.0).unwrap() * T::PI() / T::from(fft_len).unwrap();
    // index * -2PI / fft_len
    let angle = constant * T::from(index).unwrap();

    Complex::new(angle.cos(), angle.sin())
}

#[doc = concat!("Inner FFT")]
#[inline(always)]
pub fn fft3<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 3] {
    let n = 3;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let twiddle: Complex<T> = Complex::new(T::from(-0.5).unwrap(), -T::from(SQRT_3_DIV_2).unwrap());

    let xp = x[1] + x[2];
    let xn = x[1] - x[2];
    let sum = x[0] + xp;

    let temp_a = x[0]
        + Complex {
            re: twiddle.re * xp.re,
            im: twiddle.re * xp.im,
        };
    let temp_b = Complex {
        re: -twiddle.im * xn.im,
        im: twiddle.im * xn.re,
    };

    [sum, temp_a + temp_b, temp_a - temp_b]
}

#[doc = concat!("Inner FFT")]
#[inline(always)]
pub fn fft9<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 9] {
    let n = 9;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let twiddle1: Complex<T> = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle2: Complex<T> = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle4: Complex<T> = Complex::new(
        T::from(-0.93969262078590832).unwrap(),
        T::from(-0.34202014332566888).unwrap(),
    );

    let first = fft3([x[0], x[3], x[6]]);
    let second = fft3([x[1], x[4], x[7]]);
    let third = fft3([x[2], x[5], x[8]]);

    let row0 = fft3([first[0], second[0], third[0]]);
    let row1 = fft3([first[1], second[1] * twiddle1, third[1] * twiddle2]);
    let row2 = fft3([first[2], second[2] * twiddle2, third[2] * twiddle4]);

    [
        row0[0], row1[0], row2[0], row0[1], row1[1], row2[1], row0[2], row1[2], row2[2],
    ]
}

#[doc = concat!("Inner FFT")]
#[inline(always)]
pub fn fft18<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 18] {
    let n = 18;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let twiddle0 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle1 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle2 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle3 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle4 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle5 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle6 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle7 = Complex::new(
        T::from(0.93969262078590842).unwrap(),
        T::from(-0.34202014332566871).unwrap(),
    );
    let twiddle8 = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle9 = Complex::new(T::from(0.5).unwrap(), T::from(-0.8660254037844386).unwrap());
    let twiddle10 = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle11 = Complex::new(
        T::from(-0.1736481776669303).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle12 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle13 = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle14 = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle15 = Complex::new(
        T::from(-0.5).unwrap(),
        T::from(-0.86602540378443881).unwrap(),
    );
    let twiddle16 = Complex::new(
        T::from(-0.93969262078590832).unwrap(),
        T::from(-0.34202014332566888).unwrap(),
    );
    let twiddle17 = Complex::new(
        T::from(-0.93969262078590842).unwrap(),
        T::from(0.34202014332566866).unwrap(),
    );

    let row0 = fft6([x[0], x[3], x[6], x[9], x[12], x[15]]);
    let row1 = fft6([x[1], x[4], x[7], x[10], x[13], x[16]]);
    let row2 = fft6([x[2], x[5], x[8], x[11], x[14], x[17]]);

    let col0 = fft3([row0[0] * twiddle0, row1[0] * twiddle6, row2[0] * twiddle12]);
    let col1 = fft3([row0[1] * twiddle1, row1[1] * twiddle7, row2[1] * twiddle13]);
    let col2 = fft3([row0[2] * twiddle2, row1[2] * twiddle8, row2[2] * twiddle14]);
    let col3 = fft3([row0[3] * twiddle3, row1[3] * twiddle9, row2[3] * twiddle15]);
    let col4 = fft3([row0[4] * twiddle4, row1[4] * twiddle10, row2[4] * twiddle16]);
    let col5 = fft3([row0[5] * twiddle5, row1[5] * twiddle11, row2[5] * twiddle17]);

    [
        col0[0], col1[0], col2[0], col3[0], col4[0], col5[0], col0[1], col1[1], col2[1], col3[1],
        col4[1], col5[1], col0[2], col1[2], col2[2], col3[2], col4[2], col5[2],
    ]
}

#[doc = concat!("Inner FFT")]
#[inline(always)]
pub fn fft27<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 27] {
    let n = 27;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let row0 = fft9([x[0], x[3], x[6], x[9], x[12], x[15], x[18], x[21], x[24]]);
    let row1 = fft9([x[1], x[4], x[7], x[10], x[13], x[16], x[19], x[22], x[25]]);
    let row2 = fft9([x[2], x[5], x[8], x[11], x[14], x[17], x[20], x[23], x[26]]);

    let twiddle0 = Complex::new(
        T::from(0.97304487057982381).unwrap(),
        T::from(-0.23061587074244017).unwrap(),
    );
    let twiddle1 = Complex::new(
        T::from(0.89363264032341228).unwrap(),
        T::from(-0.44879918020046217).unwrap(),
    );
    let twiddle2 = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle3 = Complex::new(
        T::from(0.59715859170278618).unwrap(),
        T::from(-0.80212319275504373).unwrap(),
    );
    let twiddle4 = Complex::new(
        T::from(0.3960797660391569).unwrap(),
        T::from(-0.918216106880274).unwrap(),
    );
    let twiddle5 = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle6 = Complex::new(
        T::from(-0.058144828910475774).unwrap(),
        T::from(-0.99830815827126817).unwrap(),
    );
    let twiddle7 = Complex::new(
        T::from(-0.28680323271109021).unwrap(),
        T::from(-0.9579895123154889).unwrap(),
    );
    let twiddle8 = Complex::new(
        T::from(-0.68624163786873338).unwrap(),
        T::from(-0.72737364157304896).unwrap(),
    );
    let twiddle9 = Complex::new(
        T::from(-0.93969262078590832).unwrap(),
        T::from(-0.34202014332566888).unwrap(),
    );
    let twiddle10 = Complex::new(
        T::from(-0.99323835774194302).unwrap(),
        T::from(0.11609291412523012).unwrap(),
    );
    let twiddle11 = Complex::new(
        T::from(-0.83548781141293649).unwrap(),
        T::from(0.54950897807080601).unwrap(),
    );

    let col0 = fft3([row0[0], row1[0], row2[0]]);
    let col1 = fft3([row0[1], row1[1] * twiddle0, row2[1] * twiddle1]);
    let col2 = fft3([row0[2], row1[2] * twiddle1, row2[2] * twiddle3]);
    let col3 = fft3([row0[3], row1[3] * twiddle2, row2[3] * twiddle5]);
    let col4 = fft3([row0[4], row1[4] * twiddle3, row2[4] * twiddle7]);
    let col5 = fft3([row0[5], row1[5] * twiddle4, row2[5] * twiddle8]);
    let col6 = fft3([row0[6], row1[6] * twiddle5, row2[6] * twiddle9]);
    let col7 = fft3([row0[7], row1[7] * twiddle6, row2[7] * twiddle10]);
    let col8 = fft3([row0[8], row1[8] * twiddle7, row2[8] * twiddle11]);

    [
        col0[0], col1[0], col2[0], col3[0], col4[0], col5[0], col6[0], col7[0], col8[0], col0[1],
        col1[1], col2[1], col3[1], col4[1], col5[1], col6[1], col7[1], col8[1], col0[2], col1[2],
        col2[2], col3[2], col4[2], col5[2], col6[2], col7[2], col8[2],
    ]
}

#[doc = concat!("Inner FFT")]
#[inline(always)]
pub fn fft125<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 125] {
    let n = 125;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let twiddle0 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle1 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle2 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle3 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle4 = Complex::new(
        T::from(0.9685831611286311).unwrap(),
        T::from(-0.2486898871648548).unwrap(),
    );
    let twiddle5 = Complex::new(
        T::from(0.8763066800438636).unwrap(),
        T::from(-0.4817536741017153).unwrap(),
    );
    let twiddle6 = Complex::new(
        T::from(0.7289686274214116).unwrap(),
        T::from(-0.6845471059286887).unwrap(),
    );
    let twiddle7 = Complex::new(
        T::from(0.5358267949789965).unwrap(),
        T::from(-0.8443279255020151).unwrap(),
    );
    let twiddle8 = Complex::new(
        T::from(0.8763066800438636).unwrap(),
        T::from(-0.4817536741017153).unwrap(),
    );
    let twiddle9 = Complex::new(
        T::from(0.5358267949789965).unwrap(),
        T::from(-0.8443279255020151).unwrap(),
    );
    let twiddle10 = Complex::new(
        T::from(0.0627905195293133).unwrap(),
        T::from(-0.9980267284282716).unwrap(),
    );
    let twiddle11 = Complex::new(
        T::from(-0.4257792915650727).unwrap(),
        T::from(-0.9048270524660195).unwrap(),
    );
    let twiddle12 = Complex::new(
        T::from(0.7289686274214116).unwrap(),
        T::from(-0.6845471059286887).unwrap(),
    );
    let twiddle13 = Complex::new(
        T::from(0.0627905195293133).unwrap(),
        T::from(-0.9980267284282716).unwrap(),
    );
    let twiddle14 = Complex::new(
        T::from(-0.6374239897486897).unwrap(),
        T::from(-0.7705132427757893).unwrap(),
    );
    let twiddle15 = Complex::new(
        T::from(-0.9921147013144779).unwrap(),
        T::from(-0.1253332335643041).unwrap(),
    );
    let twiddle16 = Complex::new(
        T::from(0.5358267949789965).unwrap(),
        T::from(-0.8443279255020151).unwrap(),
    );
    let twiddle17 = Complex::new(
        T::from(-0.4257792915650727).unwrap(),
        T::from(-0.9048270524660195).unwrap(),
    );
    let twiddle18 = Complex::new(
        T::from(-0.9921147013144779).unwrap(),
        T::from(-0.1253332335643041).unwrap(),
    );
    let twiddle19 = Complex::new(
        T::from(-0.6374239897486895).unwrap(),
        T::from(0.7705132427757894).unwrap(),
    );
    let twiddle20 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle21 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle22 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle23 = Complex::new(T::from(1).unwrap(), T::from(-0).unwrap());
    let twiddle24 = Complex::new(
        T::from(0.9987369566060175).unwrap(),
        T::from(-0.050244318179769556).unwrap(),
    );
    let twiddle25 = Complex::new(
        T::from(0.9949510169813002).unwrap(),
        T::from(-0.1003617148512149).unwrap(),
    );
    let twiddle26 = Complex::new(
        T::from(0.9886517447379141).unwrap(),
        T::from(-0.15022558912075706).unwrap(),
    );
    let twiddle27 = Complex::new(
        T::from(0.9798550523842469).unwrap(),
        T::from(-0.19970998051440703).unwrap(),
    );
    let twiddle28 = Complex::new(
        T::from(0.9949510169813002).unwrap(),
        T::from(-0.1003617148512149).unwrap(),
    );
    let twiddle29 = Complex::new(
        T::from(0.9798550523842469).unwrap(),
        T::from(-0.19970998051440703).unwrap(),
    );
    let twiddle30 = Complex::new(
        T::from(0.954864544746643).unwrap(),
        T::from(-0.2970415815770349).unwrap(),
    );
    let twiddle31 = Complex::new(
        T::from(0.9202318473658704).unwrap(),
        T::from(-0.3913736668372024).unwrap(),
    );
    let twiddle32 = Complex::new(
        T::from(0.9886517447379141).unwrap(),
        T::from(-0.15022558912075706).unwrap(),
    );
    let twiddle33 = Complex::new(
        T::from(0.954864544746643).unwrap(),
        T::from(-0.2970415815770349).unwrap(),
    );
    let twiddle34 = Complex::new(
        T::from(0.8994052515663711).unwrap(),
        T::from(-0.4371157666509329).unwrap(),
    );
    let twiddle35 = Complex::new(
        T::from(0.8235325976284275).unwrap(),
        T::from(-0.5672689491267565).unwrap(),
    );
    let twiddle36 = Complex::new(
        T::from(0.9798550523842469).unwrap(),
        T::from(-0.19970998051440703).unwrap(),
    );
    let twiddle37 = Complex::new(
        T::from(0.9202318473658704).unwrap(),
        T::from(-0.3913736668372024).unwrap(),
    );
    let twiddle38 = Complex::new(
        T::from(0.8235325976284275).unwrap(),
        T::from(-0.5672689491267565).unwrap(),
    );
    let twiddle39 = Complex::new(
        T::from(0.6936533058128049).unwrap(),
        T::from(-0.7203090248879068).unwrap(),
    );
    let twiddle40 = Complex::new(
        T::from(0.9685831611286311).unwrap(),
        T::from(-0.2486898871648548).unwrap(),
    );
    let twiddle41 = Complex::new(
        T::from(0.8763066800438636).unwrap(),
        T::from(-0.4817536741017153).unwrap(),
    );
    let twiddle42 = Complex::new(
        T::from(0.7289686274214116).unwrap(),
        T::from(-0.6845471059286887).unwrap(),
    );
    let twiddle43 = Complex::new(
        T::from(0.5358267949789965).unwrap(),
        T::from(-0.8443279255020151).unwrap(),
    );
    let twiddle44 = Complex::new(
        T::from(0.954864544746643).unwrap(),
        T::from(-0.2970415815770349).unwrap(),
    );
    let twiddle45 = Complex::new(
        T::from(0.8235325976284275).unwrap(),
        T::from(-0.5672689491267565).unwrap(),
    );
    let twiddle46 = Complex::new(
        T::from(0.6178596130903343).unwrap(),
        T::from(-0.7862884321366189).unwrap(),
    );
    let twiddle47 = Complex::new(
        T::from(0.35641187871325075).unwrap(),
        T::from(-0.934328942456612).unwrap(),
    );
    let twiddle48 = Complex::new(
        T::from(0.9387338576538741).unwrap(),
        T::from(-0.34464292317451706).unwrap(),
    );
    let twiddle49 = Complex::new(
        T::from(0.7624425110114479).unwrap(),
        T::from(-0.6470559615694443).unwrap(),
    );
    let twiddle50 = Complex::new(
        T::from(0.49272734154829156).unwrap(),
        T::from(-0.8701837546695257).unwrap(),
    );
    let twiddle51 = Complex::new(
        T::from(0.1626371651948835).unwrap(),
        T::from(-0.986685944207868).unwrap(),
    );
    let twiddle52 = Complex::new(
        T::from(0.9202318473658704).unwrap(),
        T::from(-0.3913736668372024).unwrap(),
    );
    let twiddle53 = Complex::new(
        T::from(0.6936533058128049).unwrap(),
        T::from(-0.7203090248879068).unwrap(),
    );
    let twiddle54 = Complex::new(
        T::from(0.35641187871325075).unwrap(),
        T::from(-0.934328942456612).unwrap(),
    );
    let twiddle55 = Complex::new(
        T::from(-0.037690182669934576).unwrap(),
        T::from(-0.9992894726405892).unwrap(),
    );
    let twiddle56 = Complex::new(
        T::from(0.8994052515663711).unwrap(),
        T::from(-0.4371157666509329).unwrap(),
    );
    let twiddle57 = Complex::new(
        T::from(0.6178596130903343).unwrap(),
        T::from(-0.7862884321366189).unwrap(),
    );
    let twiddle58 = Complex::new(
        T::from(0.21200710992205452).unwrap(),
        T::from(-0.9772681235681935).unwrap(),
    );
    let twiddle59 = Complex::new(
        T::from(-0.23649899702372465).unwrap(),
        T::from(-0.971631732914674).unwrap(),
    );
    let twiddle60 = Complex::new(
        T::from(0.8763066800438636).unwrap(),
        T::from(-0.4817536741017153).unwrap(),
    );
    let twiddle61 = Complex::new(
        T::from(0.5358267949789965).unwrap(),
        T::from(-0.8443279255020151).unwrap(),
    );
    let twiddle62 = Complex::new(
        T::from(0.0627905195293133).unwrap(),
        T::from(-0.9980267284282716).unwrap(),
    );
    let twiddle63 = Complex::new(
        T::from(-0.4257792915650727).unwrap(),
        T::from(-0.9048270524660195).unwrap(),
    );
    let twiddle64 = Complex::new(
        T::from(0.8509944817946918).unwrap(),
        T::from(-0.5251746299612957).unwrap(),
    );
    let twiddle65 = Complex::new(
        T::from(0.44838321609003223).unwrap(),
        T::from(-0.8938414241512638).unwrap(),
    );
    let twiddle66 = Complex::new(
        T::from(-0.0878511965507432).unwrap(),
        T::from(-0.9961336091431725).unwrap(),
    );
    let twiddle67 = Complex::new(
        T::from(-0.5979049830575189).unwrap(),
        T::from(-0.8015669848708765).unwrap(),
    );
    let twiddle68 = Complex::new(
        T::from(0.8235325976284275).unwrap(),
        T::from(-0.5672689491267565).unwrap(),
    );
    let twiddle69 = Complex::new(
        T::from(0.35641187871325075).unwrap(),
        T::from(-0.934328942456612).unwrap(),
    );
    let twiddle70 = Complex::new(
        T::from(-0.23649899702372465).unwrap(),
        T::from(-0.971631732914674).unwrap(),
    );
    let twiddle71 = Complex::new(
        T::from(-0.7459411454241821).unwrap(),
        T::from(-0.6660118674342517).unwrap(),
    );
    let twiddle72 = Complex::new(
        T::from(0.7939903986478354).unwrap(),
        T::from(-0.6079302976946054).unwrap(),
    );
    let twiddle73 = Complex::new(
        T::from(0.260841506289897).unwrap(),
        T::from(-0.9653816388332739).unwrap(),
    );
    let twiddle74 = Complex::new(
        T::from(-0.3797790955218012).unwrap(),
        T::from(-0.925077206834458).unwrap(),
    );
    let twiddle75 = Complex::new(
        T::from(-0.8639234171928353).unwrap(),
        T::from(-0.5036232016357609).unwrap(),
    );
    let twiddle76 = Complex::new(
        T::from(0.7624425110114479).unwrap(),
        T::from(-0.6470559615694443).unwrap(),
    );
    let twiddle77 = Complex::new(
        T::from(0.1626371651948835).unwrap(),
        T::from(-0.986685944207868).unwrap(),
    );
    let twiddle78 = Complex::new(
        T::from(-0.5144395337815065).unwrap(),
        T::from(-0.8575266561936522).unwrap(),
    );
    let twiddle79 = Complex::new(
        T::from(-0.9470983049947443).unwrap(),
        T::from(-0.32094360980720926).unwrap(),
    );
    let twiddle80 = Complex::new(
        T::from(0.7289686274214116).unwrap(),
        T::from(-0.6845471059286887).unwrap(),
    );
    let twiddle81 = Complex::new(
        T::from(0.0627905195293133).unwrap(),
        T::from(-0.9980267284282716).unwrap(),
    );
    let twiddle82 = Complex::new(
        T::from(-0.6374239897486897).unwrap(),
        T::from(-0.7705132427757893).unwrap(),
    );
    let twiddle83 = Complex::new(
        T::from(-0.9921147013144779).unwrap(),
        T::from(-0.1253332335643041).unwrap(),
    );
    let twiddle84 = Complex::new(
        T::from(0.6936533058128049).unwrap(),
        T::from(-0.7203090248879068).unwrap(),
    );
    let twiddle85 = Complex::new(
        T::from(-0.037690182669934576).unwrap(),
        T::from(-0.9992894726405892).unwrap(),
    );
    let twiddle86 = Complex::new(
        T::from(-0.7459411454241821).unwrap(),
        T::from(-0.6660118674342517).unwrap(),
    );
    let twiddle87 = Complex::new(
        T::from(-0.9971589002606139).unwrap(),
        T::from(0.07532680552793279).unwrap(),
    );
    let twiddle88 = Complex::new(
        T::from(0.6565857557529564).unwrap(),
        T::from(-0.7542513807361038).unwrap(),
    );
    let twiddle89 = Complex::new(
        T::from(-0.13779029068463805).unwrap(),
        T::from(-0.9904614256966512).unwrap(),
    );
    let twiddle90 = Complex::new(
        T::from(-0.8375280400421417).unwrap(),
        T::from(-0.5463943467342692).unwrap(),
    );
    let twiddle91 = Complex::new(
        T::from(-0.9620276715860859).unwrap(),
        T::from(0.2729519355173252).unwrap(),
    );
    let twiddle92 = Complex::new(
        T::from(0.6178596130903343).unwrap(),
        T::from(-0.7862884321366189).unwrap(),
    );
    let twiddle93 = Complex::new(
        T::from(-0.23649899702372465).unwrap(),
        T::from(-0.971631732914674).unwrap(),
    );
    let twiddle94 = Complex::new(
        T::from(-0.9101059706849958).unwrap(),
        T::from(-0.4143755809932839).unwrap(),
    );
    let twiddle95 = Complex::new(
        T::from(-0.8881364488135446).unwrap(),
        T::from(0.45957986062148776).unwrap(),
    );
    let twiddle96 = Complex::new(
        T::from(0.5775727034222675).unwrap(),
        T::from(-0.816339250717184).unwrap(),
    );
    let twiddle97 = Complex::new(
        T::from(-0.3328195445229868).unwrap(),
        T::from(-0.9429905358928644).unwrap(),
    );
    let twiddle98 = Complex::new(
        T::from(-0.9620276715860859).unwrap(),
        T::from(-0.27295193551732505).unwrap(),
    );
    let twiddle99 = Complex::new(
        T::from(-0.7784623015670232).unwrap(),
        T::from(0.6276913612907007).unwrap(),
    );
    let twiddle100 = Complex::new(
        T::from(0.5358267949789965).unwrap(),
        T::from(-0.8443279255020151).unwrap(),
    );
    let twiddle101 = Complex::new(
        T::from(-0.4257792915650727).unwrap(),
        T::from(-0.9048270524660195).unwrap(),
    );
    let twiddle102 = Complex::new(
        T::from(-0.9921147013144779).unwrap(),
        T::from(-0.1253332335643041).unwrap(),
    );
    let twiddle103 = Complex::new(
        T::from(-0.6374239897486895).unwrap(),
        T::from(0.7705132427757894).unwrap(),
    );
    let twiddle104 = Complex::new(
        T::from(0.49272734154829156).unwrap(),
        T::from(-0.8701837546695257).unwrap(),
    );
    let twiddle105 = Complex::new(
        T::from(-0.5144395337815065).unwrap(),
        T::from(-0.8575266561936522).unwrap(),
    );
    let twiddle106 = Complex::new(
        T::from(-0.9996841892832999).unwrap(),
        T::from(0.02513009544333757).unwrap(),
    );
    let twiddle107 = Complex::new(
        T::from(-0.47070393216533246).unwrap(),
        T::from(0.8822912264349533).unwrap(),
    );
    let twiddle108 = Complex::new(
        T::from(0.44838321609003223).unwrap(),
        T::from(-0.8938414241512638).unwrap(),
    );
    let twiddle109 = Complex::new(
        T::from(-0.5979049830575189).unwrap(),
        T::from(-0.8015669848708765).unwrap(),
    );
    let twiddle110 = Complex::new(
        T::from(-0.9845643345292053).unwrap(),
        T::from(0.1750230589752761).unwrap(),
    );
    let twiddle111 = Complex::new(
        T::from(-0.2850192624699761).unwrap(),
        T::from(0.958521789017376).unwrap(),
    );
    let twiddle112 = Complex::new(
        T::from(0.4029064357136627).unwrap(),
        T::from(-0.9152411726209175).unwrap(),
    );
    let twiddle113 = Complex::new(
        T::from(-0.6753328081210245).unwrap(),
        T::from(-0.7375131173581739).unwrap(),
    );
    let twiddle114 = Complex::new(
        T::from(-0.9470983049947443).unwrap(),
        T::from(0.3209436098072095).unwrap(),
    );
    let twiddle115 = Complex::new(
        T::from(-0.08785119655074321).unwrap(),
        T::from(0.9961336091431725).unwrap(),
    );
    let twiddle116 = Complex::new(
        T::from(0.35641187871325075).unwrap(),
        T::from(-0.934328942456612).unwrap(),
    );
    let twiddle117 = Complex::new(
        T::from(-0.7459411454241821).unwrap(),
        T::from(-0.6660118674342517).unwrap(),
    );
    let twiddle118 = Complex::new(
        T::from(-0.8881364488135446).unwrap(),
        T::from(0.45957986062148776).unwrap(),
    );
    let twiddle119 = Complex::new(
        T::from(0.11285638487348157).unwrap(),
        T::from(0.9936113105200084).unwrap(),
    );

    let row0 = fft5([x[0], x[25], x[50], x[75], x[100]]);
    let row1 = fft5([x[5], x[30], x[55], x[80], x[105]]);
    let row2 = fft5([x[10], x[35], x[60], x[85], x[110]]);
    let row3 = fft5([x[15], x[40], x[65], x[90], x[115]]);
    let row4 = fft5([x[20], x[45], x[70], x[95], x[120]]);
    let row5 = fft5([x[1], x[26], x[51], x[76], x[101]]);
    let row6 = fft5([x[6], x[31], x[56], x[81], x[106]]);
    let row7 = fft5([x[11], x[36], x[61], x[86], x[111]]);
    let row8 = fft5([x[16], x[41], x[66], x[91], x[116]]);
    let row9 = fft5([x[21], x[46], x[71], x[96], x[121]]);

    let row10 = fft5([x[2], x[2 + 25], x[2 + 50], x[2 + 75], x[2 + 100]]);
    let row11 = fft5([x[2 + 5], x[2 + 30], x[2 + 55], x[2 + 80], x[2 + 105]]);
    let row12 = fft5([x[2 + 10], x[2 + 35], x[2 + 60], x[2 + 85], x[2 + 110]]);
    let row13 = fft5([x[2 + 15], x[2 + 40], x[2 + 65], x[2 + 90], x[2 + 115]]);
    let row14 = fft5([x[2 + 20], x[2 + 45], x[2 + 70], x[2 + 95], x[2 + 120]]);
    let row15 = fft5([x[2 + 1], x[2 + 26], x[2 + 51], x[2 + 76], x[2 + 101]]);
    let row16 = fft5([x[2 + 6], x[2 + 31], x[2 + 56], x[2 + 81], x[2 + 106]]);
    let row17 = fft5([x[2 + 11], x[2 + 36], x[2 + 61], x[2 + 86], x[2 + 111]]);
    let row18 = fft5([x[2 + 16], x[2 + 41], x[2 + 66], x[2 + 91], x[2 + 116]]);
    let row19 = fft5([x[2 + 21], x[2 + 46], x[2 + 71], x[2 + 96], x[2 + 121]]);

    let row20 = fft5([x[4], x[4 + 25], x[4 + 50], x[4 + 75], x[4 + 100]]);
    let row21 = fft5([x[4 + 5], x[4 + 30], x[4 + 55], x[4 + 80], x[4 + 105]]);
    let row22 = fft5([x[4 + 10], x[4 + 35], x[4 + 60], x[4 + 85], x[4 + 110]]);
    let row23 = fft5([x[4 + 15], x[4 + 40], x[4 + 65], x[4 + 90], x[4 + 115]]);
    let row24 = fft5([x[4 + 20], x[4 + 45], x[4 + 70], x[4 + 95], x[4 + 120]]);

    let col0 = fft5([
        row0[0],
        row1[0] * twiddle0,
        row2[0] * twiddle1,
        row3[0] * twiddle2,
        row4[0] * twiddle3,
    ]);
    let col1 = fft5([
        row0[1],
        row1[1] * twiddle4,
        row2[1] * twiddle5,
        row3[1] * twiddle6,
        row4[1] * twiddle7,
    ]);
    let col2 = fft5([
        row0[2],
        row1[2] * twiddle8,
        row2[2] * twiddle9,
        row3[2] * twiddle10,
        row4[2] * twiddle11,
    ]);
    let col3 = fft5([
        row0[3],
        row1[3] * twiddle12,
        row2[3] * twiddle13,
        row3[3] * twiddle14,
        row4[3] * twiddle15,
    ]);
    let col4 = fft5([
        row0[4],
        row1[4] * twiddle16,
        row2[4] * twiddle17,
        row3[4] * twiddle18,
        row4[4] * twiddle19,
    ]);

    let col5 = fft5([
        row5[0],
        row6[0] * twiddle0,
        row7[0] * twiddle1,
        row8[0] * twiddle2,
        row9[0] * twiddle3,
    ]);
    let col6 = fft5([
        row5[1],
        row6[1] * twiddle4,
        row7[1] * twiddle5,
        row8[1] * twiddle6,
        row9[1] * twiddle7,
    ]);
    let col7 = fft5([
        row5[2],
        row6[2] * twiddle8,
        row7[2] * twiddle9,
        row8[2] * twiddle10,
        row9[2] * twiddle11,
    ]);
    let col8 = fft5([
        row5[3],
        row6[3] * twiddle12,
        row7[3] * twiddle13,
        row8[3] * twiddle14,
        row9[3] * twiddle15,
    ]);
    let col9 = fft5([
        row5[4],
        row6[4] * twiddle16,
        row7[4] * twiddle17,
        row8[4] * twiddle18,
        row9[4] * twiddle19,
    ]);

    let col10 = fft5([
        row10[0],
        row11[0] * twiddle0,
        row12[0] * twiddle1,
        row13[0] * twiddle2,
        row14[0] * twiddle3,
    ]);
    let col11 = fft5([
        row10[1],
        row11[1] * twiddle4,
        row12[1] * twiddle5,
        row13[1] * twiddle6,
        row14[1] * twiddle7,
    ]);
    let col12 = fft5([
        row10[2],
        row11[2] * twiddle8,
        row12[2] * twiddle9,
        row13[2] * twiddle10,
        row14[2] * twiddle11,
    ]);
    let col13 = fft5([
        row10[3],
        row11[3] * twiddle12,
        row12[3] * twiddle13,
        row13[3] * twiddle14,
        row14[3] * twiddle15,
    ]);
    let col14 = fft5([
        row10[4],
        row11[4] * twiddle16,
        row12[4] * twiddle17,
        row13[4] * twiddle18,
        row14[4] * twiddle19,
    ]);

    let col15 = fft5([
        row15[0],
        row16[0] * twiddle0,
        row17[0] * twiddle1,
        row18[0] * twiddle2,
        row19[0] * twiddle3,
    ]);
    let col16 = fft5([
        row15[1],
        row16[1] * twiddle4,
        row17[1] * twiddle5,
        row18[1] * twiddle6,
        row19[1] * twiddle7,
    ]);
    let col17 = fft5([
        row15[2],
        row16[2] * twiddle8,
        row17[2] * twiddle9,
        row18[2] * twiddle10,
        row19[2] * twiddle11,
    ]);
    let col18 = fft5([
        row15[3],
        row16[3] * twiddle12,
        row17[3] * twiddle13,
        row18[3] * twiddle14,
        row19[3] * twiddle15,
    ]);
    let col19 = fft5([
        row15[4],
        row16[4] * twiddle16,
        row17[4] * twiddle17,
        row18[4] * twiddle18,
        row19[4] * twiddle19,
    ]);

    let col20 = fft5([
        row20[0],
        row21[0] * twiddle0,
        row22[0] * twiddle1,
        row23[0] * twiddle2,
        row24[0] * twiddle3,
    ]);
    let col21 = fft5([
        row20[1],
        row21[1] * twiddle4,
        row22[1] * twiddle5,
        row23[1] * twiddle6,
        row24[1] * twiddle7,
    ]);
    let col22 = fft5([
        row20[2],
        row21[2] * twiddle8,
        row22[2] * twiddle9,
        row23[2] * twiddle10,
        row24[2] * twiddle11,
    ]);
    let col23 = fft5([
        row20[3],
        row21[3] * twiddle12,
        row22[3] * twiddle13,
        row23[3] * twiddle14,
        row24[3] * twiddle15,
    ]);
    let col24 = fft5([
        row20[4],
        row21[4] * twiddle16,
        row22[4] * twiddle17,
        row23[4] * twiddle18,
        row24[4] * twiddle19,
    ]);

    let row0 = fft5([
        col0[0],
        col5[0] * twiddle20,
        col10[0] * twiddle21,
        col15[0] * twiddle22,
        col20[0] * twiddle23,
    ]);
    let row1 = fft5([
        col1[0],
        col6[0] * twiddle24,
        col11[0] * twiddle25,
        col16[0] * twiddle26,
        col21[0] * twiddle27,
    ]);
    let row2 = fft5([
        col2[0],
        col7[0] * twiddle28,
        col12[0] * twiddle29,
        col17[0] * twiddle30,
        col22[0] * twiddle31,
    ]);
    let row3 = fft5([
        col3[0],
        col8[0] * twiddle32,
        col13[0] * twiddle33,
        col18[0] * twiddle34,
        col23[0] * twiddle35,
    ]);
    let row4 = fft5([
        col4[0],
        col9[0] * twiddle36,
        col14[0] * twiddle37,
        col19[0] * twiddle38,
        col24[0] * twiddle39,
    ]);

    let row5 = fft5([
        col0[1],
        col5[1] * twiddle40,
        col10[1] * twiddle41,
        col15[1] * twiddle42,
        col20[1] * twiddle43,
    ]);
    let row6 = fft5([
        col1[1],
        col6[1] * twiddle44,
        col11[1] * twiddle45,
        col16[1] * twiddle46,
        col21[1] * twiddle47,
    ]);
    let row7 = fft5([
        col2[1],
        col7[1] * twiddle48,
        col12[1] * twiddle49,
        col17[1] * twiddle50,
        col22[1] * twiddle51,
    ]);
    let row8 = fft5([
        col3[1],
        col8[1] * twiddle52,
        col13[1] * twiddle53,
        col18[1] * twiddle54,
        col23[1] * twiddle55,
    ]);
    let row9 = fft5([
        col4[1],
        col9[1] * twiddle56,
        col14[1] * twiddle57,
        col19[1] * twiddle58,
        col24[1] * twiddle59,
    ]);

    let row10 = fft5([
        col0[2],
        col5[2] * twiddle60,
        col10[2] * twiddle61,
        col15[2] * twiddle62,
        col20[2] * twiddle63,
    ]);
    let row11 = fft5([
        col1[2],
        col6[2] * twiddle64,
        col11[2] * twiddle65,
        col16[2] * twiddle66,
        col21[2] * twiddle67,
    ]);
    let row12 = fft5([
        col2[2],
        col7[2] * twiddle68,
        col12[2] * twiddle69,
        col17[2] * twiddle70,
        col22[2] * twiddle71,
    ]);
    let row13 = fft5([
        col3[2],
        col8[2] * twiddle72,
        col13[2] * twiddle73,
        col18[2] * twiddle74,
        col23[2] * twiddle75,
    ]);
    let row14 = fft5([
        col4[2],
        col9[2] * twiddle76,
        col14[2] * twiddle77,
        col19[2] * twiddle78,
        col24[2] * twiddle79,
    ]);

    let row15 = fft5([
        col0[3],
        col5[3] * twiddle80,
        col10[3] * twiddle81,
        col15[3] * twiddle82,
        col20[3] * twiddle83,
    ]);
    let row16 = fft5([
        col1[3],
        col6[3] * twiddle84,
        col11[3] * twiddle85,
        col16[3] * twiddle86,
        col21[3] * twiddle87,
    ]);
    let row17 = fft5([
        col2[3],
        col7[3] * twiddle88,
        col12[3] * twiddle89,
        col17[3] * twiddle90,
        col22[3] * twiddle91,
    ]);
    let row18 = fft5([
        col3[3],
        col8[3] * twiddle92,
        col13[3] * twiddle93,
        col18[3] * twiddle94,
        col23[3] * twiddle95,
    ]);
    let row19 = fft5([
        col4[3],
        col9[3] * twiddle96,
        col14[3] * twiddle97,
        col19[3] * twiddle98,
        col24[3] * twiddle99,
    ]);

    let row20 = fft5([
        col0[4],
        col5[4] * twiddle100,
        col10[4] * twiddle101,
        col15[4] * twiddle102,
        col20[4] * twiddle103,
    ]);
    let row21 = fft5([
        col1[4],
        col6[4] * twiddle104,
        col11[4] * twiddle105,
        col16[4] * twiddle106,
        col21[4] * twiddle107,
    ]);
    let row22 = fft5([
        col2[4],
        col7[4] * twiddle108,
        col12[4] * twiddle109,
        col17[4] * twiddle110,
        col22[4] * twiddle111,
    ]);
    let row23 = fft5([
        col3[4],
        col8[4] * twiddle112,
        col13[4] * twiddle113,
        col18[4] * twiddle114,
        col23[4] * twiddle115,
    ]);
    let row24 = fft5([
        col4[4],
        col9[4] * twiddle116,
        col14[4] * twiddle117,
        col19[4] * twiddle118,
        col24[4] * twiddle119,
    ]);

    [
        row0[0], row1[0], row2[0], row3[0], row4[0], row5[0], row6[0], row7[0], row8[0], row9[0],
        row10[0], row11[0], row12[0], row13[0], row14[0], row15[0], row16[0], row17[0], row18[0],
        row19[0], row20[0], row21[0], row22[0], row23[0], row24[0], row0[1], row1[1], row2[1],
        row3[1], row4[1], row5[1], row6[1], row7[1], row8[1], row9[1], row10[1], row11[1],
        row12[1], row13[1], row14[1], row15[1], row16[1], row17[1], row18[1], row19[1], row20[1],
        row21[1], row22[1], row23[1], row24[1], row0[2], row1[2], row2[2], row3[2], row4[2],
        row5[2], row6[2], row7[2], row8[2], row9[2], row10[2], row11[2], row12[2], row13[2],
        row14[2], row15[2], row16[2], row17[2], row18[2], row19[2], row20[2], row21[2], row22[2],
        row23[2], row24[2], row0[3], row1[3], row2[3], row3[3], row4[3], row5[3], row6[3], row7[3],
        row8[3], row9[3], row10[3], row11[3], row12[3], row13[3], row14[3], row15[3], row16[3],
        row17[3], row18[3], row19[3], row20[3], row21[3], row22[3], row23[3], row24[3], row0[4],
        row1[4], row2[4], row3[4], row4[4], row5[4], row6[4], row7[4], row8[4], row9[4], row10[4],
        row11[4], row12[4], row13[4], row14[4], row15[4], row16[4], row17[4], row18[4], row19[4],
        row20[4], row21[4], row22[4], row23[4], row24[4],
    ]
}

// #[cfg(test)]
// mod tests {
//     use num_complex::Complex;

//     use crate::fft125;

//     #[test]
//     fn test_125() {
//         let v: Vec<Complex<f64>> = (0..125).map(|i| Complex::new(i as f64, 0.0)).collect();
//         let a = fft125(&v);
//         dbg!(&a);
//     }
// }
