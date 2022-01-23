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
    pub struct Operator<T: ValidOperatorType, const N: usize, const M: usize>
    {
        pub qubits: [usize; N],
        pub matrix: SMatrix<T, M, M>,
    }

    impl<T: ValidOperatorType, const N: usize, const M: usize> NQubitOperator for Operator<T, N, M> {
        type OperatorType = T;
    }

    impl<T: ValidOperatorType, const N: usize, const M: usize> Operator<T, N, M>
    {
        pub fn new(
            qubits: [usize; N],
            matrix: SMatrix<T, M, M>,
        ) -> Operator<T, N, M>
       {
            Operator { qubits, matrix }
        }
    }

    // impl <const N: usize, const M: usize>  Operator<Complex<f32>, N, M>
    // {
    //     pub fn from_pauli_string(
    //         qubits: [usize; N],
    //         pauli_string: &str,
    //     ) -> Operator<Complex<f32>, N, M>
    //     {
    //         if pauli_string.len() != N {
    //             panic!("Invalid Pauli string length");
    //         }
    //         if N != M {
    //             panic!("Invalid Pauli string size");
    //         }
    //         let matrix = operator_from_pauli_string![pauli_string.chars()];
    //         Operator { qubits, matrix }
    //     }
    // }

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
}
