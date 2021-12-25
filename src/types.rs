mod types {
    use std::error::Error;
    use std::ops::Add;
    use nalgebra::{Complex, Matrix2, Matrix4, SMatrix};
    use nalgebra as na;

    pub const X: Matrix2<Complex<f32>> = Matrix2::new(Complex::new(0.0, 0.0), Complex::new(1.0, 0.0),Complex::new(1.0, 0.0), Complex::new(0.0, 0.0));
    pub const Y: Matrix2<Complex<f32>> = Matrix2::new(Complex::new(0.0, 0.0), Complex::new(0.0, 1.0),Complex::new(0.0, -1.0), Complex::new(0.0, 0.0));
    pub const Z: Matrix2<Complex<f32>> = Matrix2::new(Complex::new(1.0, 0.0), Complex::new(0.0, 0.0),Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0));

    #[derive(Debug, Clone, PartialEq)]
    pub struct Operator<R: Add> {
        pub qubits: Vec<usize>,
        pub matrix: SMatrix<Complex<R>, { 2_usize.pow(qubits) }, { 2_usize.pow(qubits) }>,
    }

    impl<R: Add> Operator<R> {
        pub fn new(qubits: Vec<usize>, matrix: SMatrix<Complex<R>, { 2_usize.pow(qubits) }, { 2_usize.pow(qubits) }>) -> Operator<R> {
            Operator {
                    qubits,
                    matrix,
                }
        }
    }

}

#[cfg(test)]
mod tests {
    use crate::types::types::*;
    use nalgebra::{Complex, Matrix2, Matrix4, SMatrix};
    use nalgebra as na;

    #[test]
    fn test_pauli_commutation(){
        let commute_y = Z * X;
        assert_eq!(commute_y * Complex::new(0.0, 1.0) , Y);
        let commute_z = X * Y;
        assert_eq!(commute_z * Complex::new(0.0, 1.0) , Z);
        let commute_x = Y * Z;
        assert_eq!(commute_x * Complex::new(0.0, 1.0) , X);
    }

    #[test]
    fn test_single_qubit_operator(){
        let x = Operator::new(vec![0], X);
    }
}