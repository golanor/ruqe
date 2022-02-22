mod operators {
    use std::fmt::Debug;
    use nalgebra::{Complex, Matrix2, SMatrix, zero};
    use num_traits::{One, Zero};

    pub trait ValidOperatorType: Clone + PartialEq + Debug + Zero + One{}

    impl ValidOperatorType for Complex<f64> {}
    impl ValidOperatorType for Complex<f32> {}
    impl ValidOperatorType for f64 {}
    impl ValidOperatorType for f32 {}

    pub trait NQubitOperator{
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
    pub struct Operator<'a, T: ValidOperatorType, const M: usize>
    {
        pub qubits: &'a[usize],
        pub matrix: SMatrix<T, M, M>,
    }

    impl<T: ValidOperatorType, const M: usize> NQubitOperator for Operator<'_, T, M> {
        type OperatorType = T;
    }

    impl<T: ValidOperatorType, const M: usize> Operator<'_, T, M>
    {
        pub fn new(
            qubits: &[usize],
            matrix: SMatrix<T, M, M>,
        ) -> Operator<T, M>
       {
            Operator { qubits, matrix }
        }

        pub fn embed_in_hilbert_space<const N: usize>(
            &self,
        ) -> SMatrix<T, N, N>
        {
            let mut result = Vec::with_capacity(N * N);
            let row_shift = N;
            for q in 0..N {
                let index = (2 as usize).pow(q as u32);
                let mat = if self.qubits.contains(&q) {
                    let qubit_index = self.qubits.iter().position(|&x| x == q).unwrap();
                    self.matrix.fixed_slice::<2, 2>(2*qubit_index, 2*qubit_index).clone_owned()
                } else {
                    SMatrix::<T, 2, 2>::identity()
                };
                result.splice(index..=(index + 1), mat.row(0));
                result.splice((index + (q + 1) * row_shift)..=(index + (q+1)*row_shift + 1), mat.row(1));
            }
            SMatrix::from_row_slice(result.as_slice())
        }
    }

    pub struct StaticHamiltonian<O: NQubitOperator, const M: usize> {
        pub matrix: SMatrix<O::OperatorType, M, M>,
        pub qubits: usize,
    }

    impl<O: NQubitOperator, const M: usize> StaticHamiltonian<O, M> {
        pub fn new(matrix: SMatrix<O::OperatorType, M, M>) -> StaticHamiltonian<O, M> {
            StaticHamiltonian { matrix, qubits: (M as u32).log2() as usize }
        }
        pub fn from_operator_sum(operators: &[O]) -> StaticHamiltonian<O, M> {
            let mut matrix = zero::<SMatrix<O::OperatorType, M, M>>();
            for operator in operators {
                matrix += &operator.matrix;
            }
            StaticHamiltonian::new(matrix)
        }
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
        let _x = Operator::new(&[0], X);
    }

    #[test]
    fn test_two_qubit_operator() {
        let _xy = Operator::new(&[0, 1], X.kronecker(&Y));
    }
}
