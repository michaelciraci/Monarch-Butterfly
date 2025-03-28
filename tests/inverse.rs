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

macro_rules! compare_against_rustifft {
    ($test:ident, $n:expr) => {
        #[test]
        fn $test() {
            let mut p = rustfft::FftPlanner::new();
            let plan = p.plan_fft_inverse($n);
            let mut buf: Vec<_> = (0..$n)
                .map(|i| Complex::<f64>::new(i as f64, 0.0))
                .collect();
            let monarch = ifft::<$n, _, _>(&buf);
            plan.process(&mut buf);
            assert_slice_equal!(monarch, buf);
        }
    };
}

compare_against_rustifft!(inverse_1, 1);
compare_against_rustifft!(inverse_2, 2);
compare_against_rustifft!(inverse_3, 3);
compare_against_rustifft!(inverse_4, 4);
compare_against_rustifft!(inverse_5, 5);
compare_against_rustifft!(inverse_6, 6);
compare_against_rustifft!(inverse_7, 7);
compare_against_rustifft!(inverse_8, 8);
compare_against_rustifft!(inverse_9, 9);
compare_against_rustifft!(inverse_10, 10);
compare_against_rustifft!(inverse_11, 11);
compare_against_rustifft!(inverse_12, 12);
compare_against_rustifft!(inverse_13, 13);
compare_against_rustifft!(inverse_14, 14);
compare_against_rustifft!(inverse_15, 15);
compare_against_rustifft!(inverse_16, 16);
compare_against_rustifft!(inverse_17, 17);
compare_against_rustifft!(inverse_18, 18);
compare_against_rustifft!(inverse_19, 19);
compare_against_rustifft!(inverse_20, 20);
compare_against_rustifft!(inverse_21, 21);
compare_against_rustifft!(inverse_22, 22);
compare_against_rustifft!(inverse_23, 23);
compare_against_rustifft!(inverse_24, 24);
compare_against_rustifft!(inverse_25, 25);
compare_against_rustifft!(inverse_26, 26);
compare_against_rustifft!(inverse_27, 27);
compare_against_rustifft!(inverse_28, 28);
compare_against_rustifft!(inverse_29, 29);
compare_against_rustifft!(inverse_30, 30);
compare_against_rustifft!(inverse_31, 31);
compare_against_rustifft!(inverse_32, 32);
compare_against_rustifft!(inverse_33, 33);
compare_against_rustifft!(inverse_34, 34);
compare_against_rustifft!(inverse_35, 35);
compare_against_rustifft!(inverse_36, 36);
compare_against_rustifft!(inverse_37, 37);
compare_against_rustifft!(inverse_38, 38);
compare_against_rustifft!(inverse_39, 39);
compare_against_rustifft!(inverse_40, 40);
compare_against_rustifft!(inverse_41, 41);
compare_against_rustifft!(inverse_42, 42);
compare_against_rustifft!(inverse_43, 43);
compare_against_rustifft!(inverse_44, 44);
compare_against_rustifft!(inverse_45, 45);
compare_against_rustifft!(inverse_46, 46);
compare_against_rustifft!(inverse_47, 47);
compare_against_rustifft!(inverse_48, 48);
compare_against_rustifft!(inverse_49, 49);
compare_against_rustifft!(inverse_50, 50);
compare_against_rustifft!(inverse_51, 51);
compare_against_rustifft!(inverse_52, 52);
compare_against_rustifft!(inverse_53, 53);
compare_against_rustifft!(inverse_54, 54);
compare_against_rustifft!(inverse_55, 55);
compare_against_rustifft!(inverse_56, 56);
compare_against_rustifft!(inverse_57, 57);
compare_against_rustifft!(inverse_58, 58);
compare_against_rustifft!(inverse_59, 59);
compare_against_rustifft!(inverse_60, 60);
compare_against_rustifft!(inverse_61, 61);
compare_against_rustifft!(inverse_62, 62);
compare_against_rustifft!(inverse_63, 63);
compare_against_rustifft!(inverse_64, 64);
compare_against_rustifft!(inverse_65, 65);
compare_against_rustifft!(inverse_66, 66);
compare_against_rustifft!(inverse_67, 67);
compare_against_rustifft!(inverse_68, 68);
compare_against_rustifft!(inverse_69, 69);
compare_against_rustifft!(inverse_70, 70);
compare_against_rustifft!(inverse_71, 71);
compare_against_rustifft!(inverse_72, 72);
compare_against_rustifft!(inverse_73, 73);
compare_against_rustifft!(inverse_74, 74);
compare_against_rustifft!(inverse_75, 75);
compare_against_rustifft!(inverse_76, 76);
compare_against_rustifft!(inverse_77, 77);
compare_against_rustifft!(inverse_78, 78);
compare_against_rustifft!(inverse_79, 79);
compare_against_rustifft!(inverse_80, 80);
compare_against_rustifft!(inverse_81, 81);
compare_against_rustifft!(inverse_82, 82);
compare_against_rustifft!(inverse_83, 83);
compare_against_rustifft!(inverse_84, 84);
compare_against_rustifft!(inverse_85, 85);
compare_against_rustifft!(inverse_86, 86);
compare_against_rustifft!(inverse_87, 87);
compare_against_rustifft!(inverse_88, 88);
compare_against_rustifft!(inverse_89, 89);
compare_against_rustifft!(inverse_90, 90);
compare_against_rustifft!(inverse_91, 91);
compare_against_rustifft!(inverse_92, 92);
compare_against_rustifft!(inverse_93, 93);
compare_against_rustifft!(inverse_94, 94);
compare_against_rustifft!(inverse_95, 95);
compare_against_rustifft!(inverse_96, 96);
compare_against_rustifft!(inverse_97, 97);
compare_against_rustifft!(inverse_98, 98);
compare_against_rustifft!(inverse_99, 99);
compare_against_rustifft!(inverse_100, 100);
compare_against_rustifft!(inverse_101, 101);
compare_against_rustifft!(inverse_102, 102);
compare_against_rustifft!(inverse_103, 103);
compare_against_rustifft!(inverse_104, 104);
compare_against_rustifft!(inverse_105, 105);
