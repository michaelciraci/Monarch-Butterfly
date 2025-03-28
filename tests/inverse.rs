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
compare_against_rustifft!(inverse_106, 106);
compare_against_rustifft!(inverse_107, 107);
compare_against_rustifft!(inverse_108, 108);
compare_against_rustifft!(inverse_109, 109);
compare_against_rustifft!(inverse_110, 110);
compare_against_rustifft!(inverse_111, 111);
compare_against_rustifft!(inverse_112, 112);
compare_against_rustifft!(inverse_113, 113);
compare_against_rustifft!(inverse_114, 114);
compare_against_rustifft!(inverse_115, 115);
compare_against_rustifft!(inverse_116, 116);
compare_against_rustifft!(inverse_117, 117);
compare_against_rustifft!(inverse_118, 118);
compare_against_rustifft!(inverse_119, 119);
compare_against_rustifft!(inverse_120, 120);
compare_against_rustifft!(inverse_121, 121);
compare_against_rustifft!(inverse_122, 122);
compare_against_rustifft!(inverse_123, 123);
compare_against_rustifft!(inverse_124, 124);
compare_against_rustifft!(inverse_125, 125);
compare_against_rustifft!(inverse_126, 126);
compare_against_rustifft!(inverse_127, 127);
compare_against_rustifft!(inverse_128, 128);
compare_against_rustifft!(inverse_129, 129);
compare_against_rustifft!(inverse_130, 130);
compare_against_rustifft!(inverse_131, 131);
compare_against_rustifft!(inverse_132, 132);
compare_against_rustifft!(inverse_133, 133);
compare_against_rustifft!(inverse_134, 134);
compare_against_rustifft!(inverse_135, 135);
compare_against_rustifft!(inverse_136, 136);
compare_against_rustifft!(inverse_137, 137);
compare_against_rustifft!(inverse_138, 138);
compare_against_rustifft!(inverse_139, 139);
compare_against_rustifft!(inverse_140, 140);
compare_against_rustifft!(inverse_141, 141);
compare_against_rustifft!(inverse_142, 142);
compare_against_rustifft!(inverse_143, 143);
compare_against_rustifft!(inverse_144, 144);
compare_against_rustifft!(inverse_145, 145);
compare_against_rustifft!(inverse_146, 146);
compare_against_rustifft!(inverse_147, 147);
compare_against_rustifft!(inverse_148, 148);
compare_against_rustifft!(inverse_149, 149);
compare_against_rustifft!(inverse_150, 150);
compare_against_rustifft!(inverse_151, 151);
compare_against_rustifft!(inverse_152, 152);
compare_against_rustifft!(inverse_153, 153);
compare_against_rustifft!(inverse_154, 154);
compare_against_rustifft!(inverse_155, 155);
compare_against_rustifft!(inverse_156, 156);
compare_against_rustifft!(inverse_157, 157);
compare_against_rustifft!(inverse_158, 158);
compare_against_rustifft!(inverse_159, 159);
compare_against_rustifft!(inverse_160, 160);
compare_against_rustifft!(inverse_161, 161);
compare_against_rustifft!(inverse_162, 162);
compare_against_rustifft!(inverse_163, 163);
compare_against_rustifft!(inverse_164, 164);
compare_against_rustifft!(inverse_165, 165);
compare_against_rustifft!(inverse_166, 166);
compare_against_rustifft!(inverse_167, 167);
compare_against_rustifft!(inverse_168, 168);
compare_against_rustifft!(inverse_169, 169);
compare_against_rustifft!(inverse_170, 170);
compare_against_rustifft!(inverse_171, 171);
compare_against_rustifft!(inverse_172, 172);
compare_against_rustifft!(inverse_173, 173);
compare_against_rustifft!(inverse_174, 174);
compare_against_rustifft!(inverse_175, 175);
compare_against_rustifft!(inverse_176, 176);
compare_against_rustifft!(inverse_177, 177);
compare_against_rustifft!(inverse_178, 178);
compare_against_rustifft!(inverse_179, 179);
compare_against_rustifft!(inverse_180, 180);
compare_against_rustifft!(inverse_181, 181);
compare_against_rustifft!(inverse_182, 182);
compare_against_rustifft!(inverse_183, 183);
compare_against_rustifft!(inverse_184, 184);
compare_against_rustifft!(inverse_185, 185);
compare_against_rustifft!(inverse_186, 186);
compare_against_rustifft!(inverse_187, 187);
compare_against_rustifft!(inverse_188, 188);
compare_against_rustifft!(inverse_189, 189);
compare_against_rustifft!(inverse_190, 190);
compare_against_rustifft!(inverse_191, 191);
compare_against_rustifft!(inverse_192, 192);
compare_against_rustifft!(inverse_193, 193);
compare_against_rustifft!(inverse_194, 194);
compare_against_rustifft!(inverse_195, 195);
compare_against_rustifft!(inverse_196, 196);
compare_against_rustifft!(inverse_197, 197);
compare_against_rustifft!(inverse_198, 198);
compare_against_rustifft!(inverse_199, 199);
compare_against_rustifft!(inverse_200, 200);
