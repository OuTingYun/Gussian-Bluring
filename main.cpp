#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <bits/stdc++.h>
#include <cblas.h>

using namespace std;


class Matrix {  

public:
    vector<double> m_buffer;
    size_t n_row;
    size_t n_col;
    Matrix(size_t nrow, size_t ncol) : n_row(nrow), n_col(ncol){


        m_buffer.resize(nrow * ncol);
        for (size_t i=0;i<nrow;i++)
            for (size_t j=0;j<ncol;j++)
                m_buffer[i*nrow+j]=0;
    }



    size_t nrow() const { return n_row; }
    size_t ncol() const { return n_col; }

    double operator()(size_t row, size_t col)const{
        return m_buffer[row * n_col + col];

    }
    double &operator()(size_t row, size_t col){
			return m_buffer[row * n_col + col];
	}
    bool operator==(Matrix const &matrix) const{
        for (size_t i=0;i<n_row;i++)
            for(size_t j=0;j<n_col;j++)
                if (m_buffer[i*n_row+j]!=matrix(i,j))
                    return false;
        return true;
    }
    double* addr() { return m_buffer.data(); }

    float sum(){
        float count=0;
        for(size_t i=0;i<n_row;i++){
            for(size_t j=0;j<n_col;j++){
                count+= m_buffer[i * n_col + j];
            }
        }
        return count;
    }

};

Matrix multiply_naive(Matrix& mat1, Matrix& mat2){
	Matrix res(mat1.nrow(), mat2.ncol());
	for(size_t k=0; k < mat2.ncol(); k++)
		for(size_t i=0; i < mat1.nrow(); i++)
			for(size_t j=0; j < mat1.ncol(); j++)
				res(i, k) += mat1(i,j) * mat2(j,k);

	return res;
}

Matrix multiply_tile(Matrix const & matrix1, Matrix const & matrix2, size_t tilesize)
{
    

    size_t row_max = matrix1.nrow();
    size_t col_max = matrix2.ncol();
    size_t inner_max = matrix1.ncol();

    if(row_max != col_max){
        throw out_of_range("mat1 col is different from mat2 row");
    }

    Matrix ret(row_max, col_max);
    
    // Run for every block.
    for (size_t row = 0; row < row_max; row += tilesize){
        size_t imax = min(row_max, row + tilesize);
    
        for (size_t col = 0; col < col_max; col += tilesize) {
            size_t jmax = min(col_max, col + tilesize);
            
            for (size_t inner = 0; inner < inner_max; inner += tilesize) {
                size_t kmax = min(inner_max, inner + tilesize);

                //Runing inside the block
                for (size_t i = row; i < imax; ++i){
                    size_t base1 = i * inner_max;

                    for (size_t j = col; j < jmax; ++j){
                        size_t index = i * col_max + j;

                        for (size_t k = inner; k < kmax; ++k) {   
                            ret.m_buffer[index] += matrix1.m_buffer[base1 + k] * matrix2(k, j);
                        }
                    }
                }
            }
        }
    }

    return ret;

}
Matrix multiply_mkl(Matrix& mat1, Matrix& mat2)
{
    if(mat1.ncol() != mat2.nrow()){
        throw out_of_range("mat1 col is different from mat2 row");
    }

    Matrix ret(mat1.nrow(), mat2.ncol());

    cblas_dgemm(
        CblasRowMajor, 
        CblasNoTrans, 
        CblasNoTrans, 
        mat1.nrow(), 
        mat2.ncol(), 
        mat1.ncol(), 
        1.0, 
        mat1.addr(), 
        mat1.ncol(), 
        mat2.addr(), 
        mat2.ncol(), 
        0.0, 
        ret.addr(), 
        mat2.ncol());

    return ret;
}

Matrix convolute(Matrix mat1, Matrix mat2){
    int f_size = mat2.n_row;
    Matrix mat1_(f_size,f_size);
    Matrix mat_conv(mat1.n_row-mat2.n_row+1,mat1.n_col-mat2.n_col+1);

    for(size_t i=0; i<mat_conv.n_row; i++){
        for(size_t j=0;j<mat_conv.n_col; j++){
            for(int inner_row = 0; inner_row<f_size ;inner_row++){
                for(int inner_col = 0; inner_col<f_size ; inner_col++){
                    mat1_(inner_row,inner_col) = mat1(i+inner_row,j+inner_col)*mat2(inner_row,inner_col);   
                }
                
            }
            mat_conv(i,j) = mat1_.sum();
        }
    }
    return mat_conv;
    
}
Matrix getGaussian(int f_size, double sigma)
{
    Matrix kernel(f_size, f_size);
    double sum=0.0;
    int center = (f_size-1)/2;
    for(int i=0-center;i<=center;i++){
        for (int j=0-center;j<=center;j++){
            kernel(i+center,j+center) = exp(-(i*i+j*j)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
            sum += kernel(i+center,j+center);
        }
    }


    for (int i=0 ; i<f_size ; i++) {
        for (int j=0 ; j<f_size ; j++) {
            kernel(i,j)/= sum;
        }
    }

    return kernel;
}

Matrix getBoxBlur(int f_size){
    Matrix kernel(f_size, f_size);
    for(int i=0;i<f_size;i++){
        for(int j=0;j<f_size;j++){
            kernel(i,j)=1;
        }
    }
    for (int i=0 ; i<f_size ; i++) {
        for (int j=0 ; j<f_size ; j++) {
            kernel(i,j)/= (f_size*f_size);
        }
    }
    return kernel;
}
Matrix BoxBlur(Matrix mat_raw,int f_size){
    int row = mat_raw.n_row;
    int col = mat_raw.n_col;
    int pad = (f_size-1);
    Matrix mat1(row+pad,col+pad);
    Matrix mat2 = getBoxBlur(f_size);

    for(int i=pad/2;i<row+pad/2;i++){ 
        for(int j=pad/2;j<col+pad/2;j++){
            mat1(i,j) = mat_raw(i-1,j-1);
        }
    }
    mat_raw = convolute(mat1,mat2);
    return mat_raw;

}

Matrix Gaussian(Matrix mat_raw,int f_size,double sigma)
{

    int row = mat_raw.n_row;
    int col = mat_raw.n_col;
    int pad = (f_size-1);
    Matrix mat1(row+pad,col+pad);
    Matrix mat2 = getGaussian(f_size,sigma);

    for(int i=pad/2;i<row+pad/2;i++){ 
        for(int j=pad/2;j<col+pad/2;j++){
            mat1(i,j) = mat_raw(i-1,j-1);
        }
    }
    mat_raw = convolute(mat1,mat2);
    return mat_raw;
}

PYBIND11_MODULE(main, m) {
    pybind11::class_<Matrix>(m, "Matrix")
        .def(pybind11::init<size_t, size_t>())
        .def("__getitem__", [](const Matrix &mat, array<int, 2> i){ return mat(i[0], i[1]);}) // lambda function.
        .def("__setitem__",[](Matrix &mat, pair<size_t, size_t> idx, double val) { return mat(idx.first, idx.second) = val; })
        .def("__eq__", [](const Matrix &mat, const Matrix &other) { return mat == other; })
        .def_property_readonly("nrow", &Matrix::nrow)
        .def_property_readonly("ncol", &Matrix::ncol);
    m.def("Gaussian",&Gaussian);
    m.def("getGaussian",&getGaussian);
    m.def("BoxBlur",&BoxBlur);
    m.def("getBoxBlur",&getBoxBlur);
    m.def("multiply_naive",&multiply_naive);
    m.def("multiply_tile",&multiply_tile);
    m.def("multiply_mkl",&multiply_mkl);

}
