use monarch_butterfly::{fft, ifft};
use num_complex::Complex;

#[test]
fn rejects_wrong_input_lengths_on_each_dispatch_path() {
    fn check<const N: usize>() {
        let short = [Complex::new(0.0_f32, 0.0); N];
        assert!(std::panic::catch_unwind(|| fft::<N, _, _>(&short[..N - 1])).is_err());
        assert!(std::panic::catch_unwind(|| ifft::<N, _, _>(&short[..N - 1])).is_err());
    }
    check::<3>();
    check::<5>();
    check::<8>();
    check::<48>();
    check::<64>();
    check::<97>();
    check::<15>();
}

#[test]
fn rejects_unsupported_sizes() {
    assert!(std::panic::catch_unwind(|| fft::<0, _, _>([Complex::new(0.0_f32, 0.0); 0])).is_err());
    assert!(
        std::panic::catch_unwind(|| ifft::<201, _, _>([Complex::new(0.0_f32, 0.0); 201])).is_err()
    );
}
