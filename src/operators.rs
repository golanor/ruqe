mod operators {
    use nalgebra::{Complex, Matrix2, OMatrix, SMatrix, SVector, zero};
    use std::ops::AddAssign;

    pub trait ValidOperatorType {}

    impl ValidOperatorType for Complex<f64> {}
    impl ValidOperatorType for Complex<f32> {}
    impl ValidOperatorType for f64 {}
    impl ValidOperatorType for f32 {}

    pub trait NQubitOperator {
        type OperatorType: ValidOperatorType;
    }

    pub const X: Matrix2<Complex<f32>> = Matrix2::new(
        Complex::new(0.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 0.0),
    );
    pub const Y: Matrix2<Complex<f32>> = Matrix2::new(
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 1.0),
        Complex::new(0.0, -1.0),
        Complex::new(0.0, 0.0),
    );
    pub const Z: Matrix2<Complex<f32>> = Matrix2::new(
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(-1.0, 0.0),
    );
    pub const I: Matrix2<Complex<f32>> = Matrix2::new(
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(1.0, 0.0),
    );

    #[derive(Debug, Clone)]
    pub struct Operator<T: ValidOperatorType, const N: usize>
    where
        SMatrix<Complex<T>, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>: Sized,
    {
        pub qubits: [usize; N],
        pub matrix: SMatrix<T, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>,
    }

    impl<T: ValidOperatorType, const N: usize> NQubitOperator for Operator<T, N>
    where
        SMatrix<T, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>: Sized,
    {
        type OperatorType = T;
    }

    impl<T: ValidOperatorType, const N: usize> Operator<T, N>
    where
        SMatrix<T, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>: Sized,
    {
        pub fn new(
            qubits: [usize; N],
            matrix: SMatrix<T, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>,
        ) -> Operator<T, N>
        where
            SMatrix<T, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>: Sized,
        {
            Operator { qubits, matrix }
        }
    }

    impl <const N: usize>  Operator<Complex<f32>, N>
    where
        SMatrix<Complex<f32>, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>: Sized,
    {
        pub fn from_pauli_string(
            qubits: [usize; N],
            pauli_string: &str,
        ) -> Operator<Complex<f32>, N>
        where
            SMatrix<Complex<f32>, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>: Sized
        {
            let mut matrix = zero::<SMatrix<Complex<f32>, { 2_usize.pow(N as u32) }, { 2_usize.pow(N as u32) }>>();
            let blocks = pauli_string.to_uppercase().chars().map(|x| match x {
                'X' => &X,
                'Y' => &Y,
                'Z' => &Z,
                'I' => &I,
                _ => panic!("Invalid Pauli string"),
            }).collect::<Vec<_>>();
            let mut i = 0_usize;
            for block in blocks {
                let block_size = block.nrows();
                matrix.slice_mut((i, i), (block_size, block_size)).copy_from(block);
                i += block_size;
            }
            Operator { qubits, matrix }
        }
    }

    pub struct StaticHamiltonian<O: NQubitOperator, const N: usize> {
        pub terms: [Vec<O>; N],
    }

    impl<O: NQubitOperator, const N: usize> StaticHamiltonian<O, N> {
        pub fn new(terms: [Vec<O>; N]) -> StaticHamiltonian<O, N> {
            StaticHamiltonian { terms }
        }

        // pub fn evolve(
        //     &self,
        //     state: &mut SVector<O::OperatorType, { 2_usize.pow(N as u32) }>,
        // ) -> SVector<O::OperatorType, { 2_usize.pow(N as u32) }> {
        //     let mut result = zero::<SVector<O::OperatorType, { 2_usize.pow(N as u32) }>>();
        //     for term in self.terms.iter() {
        //         for operator in term.iter() {
        //             let mut new_state = state.clone();
        //             new_state.apply_matrix(&operator.matrix);
        //             result.add_assign(&new_state);
        //         }
        //     }
        //     result
        // }
    }
}

#[cfg(test)]
mod tests {
    use crate::operators::operators::*;
    use nalgebra::Complex;

    #[test]
    fn test_pauli_commutation() {
        let commute_y = Z * X;
        assert_eq!(commute_y * Complex::new(0.0, 1.0), Y);
        let commute_z = X * Y;
        assert_eq!(commute_z * Complex::new(0.0, 1.0), Z);
        let commute_x = Y * Z;
        assert_eq!(commute_x * Complex::new(0.0, 1.0), X);
    }

    #[test]
    fn test_single_qubit_operator() {
        let _x = Operator::new([0], X);
    }

    #[test]
    fn test_two_qubit_operator() {
        let _xy = Operator::new([0, 1], X.kronecker(&Y));
    }

    #[test]
    fn test_from_pauli_string() {
        let _x = Operator::from_pauli_string([0], "X");
        let _xy = Operator::from_pauli_string([0, 1], "XY");
        assert_eq!(_xy.matrix, X.kronecker(&Y));
        assert_eq!(_x.matrix, X);
    }
}
