/*
 * conv.h Matlab-like convolution
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_CONV_H
#define ES_FLOW_CONV_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>
#include <math.h>

//#include "cpplot.h"
//using namespace cpplot;

namespace utilities {


/** @brief Find the next good size for an fft, to pad with minimal number of zeros.
 *
 * @param N Length of a signal
 * @return M Length to pad the signal to, for optimal FFT performance
 */
template<typename T>
T fft_next_good_size(const T n) {
    T result = n;
    if (n <= 2) {
        result = 2;
        return result;
    }
    while (true) {
        T m = result;
        while ((m % 2) == 0) m = m / 2;
        while ((m % 3) == 0) m = m / 3;
        while ((m % 5) == 0) m = m / 5;
        if (m <= 1) {
            return (result);
        }
        result = result + 1;
    }
}


/** @brief Convolve an input vector with a kernel, returning an output of the same length as the input.
 *
 * Uses zero-padded fft based output
 *
 *
 * @param out
 * @param input
 * @param kernel
 */
Eigen::VectorXd conv(const Eigen::VectorXd &input, const Eigen::VectorXd &kernel) {

    // TODO template this function signature to also accept arrays

    // Map the input signature data to tensors (shared memory)
    auto input_len = input.rows();
    auto kernel_len = kernel.rows();

    // Compute cumulative length of input and kernel;
    auto N = input_len + kernel_len - 1;
    auto M = fft_next_good_size(N);

    // Create padded FFT inputs
    Eigen::FFT<double> fft;
    Eigen::VectorXd input_padded(M);
    Eigen::VectorXd kernel_padded(M);
    input_padded.setZero();
    kernel_padded.setZero();
    input_padded.topRows(input_len) = input;
    kernel_padded.topRows(kernel_len) = kernel;

    // Take the forward ffts
    Eigen::VectorXcd input_transformed(M);
    Eigen::VectorXcd kernel_transformed(M);
    fft.fwd(input_transformed, input_padded);
    fft.fwd(kernel_transformed, kernel_padded);

    // Convolve by element-wise multiplication
    Eigen::VectorXcd inter = input_transformed.array() * kernel_transformed.array();

    // Inverse FFT
    Eigen::VectorXd out(M);
    out.setZero();
    fft.inv(out, inter);

    // Crop to the size of the input vector
    out = out.topRows(input_len);
    return out;
}


/** @brief Determines the convolution matrix C for impulse response (kernel) vector `k`.
 *
 *  convolution_matrix(k, n) * x is equivalent to conv(k, x)
 *
 *  C is a toeplitz matrix where the upper diagonal is entirely zero, the lower diagonal is zero for cases where
 *
 *  TODO sparse alternative. Large signals produce epically big, mostly empty, matrices.
 *
 *  Example:
 *  \code
 *      Eigen::VectorXd k(5)
 *      Eigen::VectorXd x(7)
 *      k << 1, 2, 3, 2, 1;
 *      x << 1, 2, 1, 2, 1, 2, 1;
 *      Eigen::MatrixXd c = convolution_matrix(k, 7);
 *      std::cout << "Convolved result c * x: " << c * x << std::endl
 *  \endcode
 *
 * @param[in]  k Column vector input
 * @param[in]  n Desired dimension of the output matrix. i.e. if convolving with vector x of length 8, n should be 8.
 * @param[out] c Convolution matrix
 */
Eigen::MatrixXd convolution_matrix(const Eigen::VectorXd &k, const Eigen::Index n) {

    Eigen::Index mk = k.rows();
    Eigen::VectorXd toeplitz_col = Eigen::VectorXd::Zero(mk + n - 1);
    toeplitz_col.topRows(mk) = k;

    Eigen::Index mt = toeplitz_col.rows();
    Eigen::MatrixXd toeplitz = Eigen::MatrixXd::Zero(mt, n);

    // Fill the lower diagonal of a toeplitz matrix.
    // If it's fat rectangular we don't need to fill the right-most columns.
    Eigen::Index nt_limit = Eigen::Vector2i(n, mt).minCoeff();
    for (Eigen::Index j = 0; j < nt_limit; j++) {
        Eigen::Index extent = Eigen::Vector2i(mt-j, k.rows()).minCoeff();
        toeplitz.block(j,j, extent, 1) = toeplitz_col.topRows(extent);
    }
    return toeplitz;
}


/** @brief stably solve Fredholm Integral of the first kind by deconvolution with diagonal loading
 *
 * Useful site: http://eeweb.poly.edu/iselesni/lecture_notes/least_squares/LeastSquares_SPdemos/deconvolution/html/deconv_demo.html
 *
 * @param[in] input Input signal
 * @param[in] kernel Kernel (impulse response) to deconvolve out of the input signal
 * @param[in] alpha  Stabilisation parameter alpha. Default 0.1.
 * @return deconvolved signal, same length as input signal.
 */
Eigen::VectorXd diagonal_loading_deconv(const Eigen::VectorXd &input, const Eigen::VectorXd &kernel, const double alpha=0.1) {
    Eigen::Index n = input.rows();
    Eigen::MatrixXd eye(n, n);
    eye.setIdentity();
    Eigen::MatrixXd c = convolution_matrix(kernel, n).topRows(n);
    Eigen::MatrixXd ct = c.transpose();
    Eigen::BDCSVD<Eigen::MatrixXd> svd(ct * c + alpha * eye, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd x = svd.solve(ct * input);
    return x;
}


/** @brief stably solve Fredholm Integral of the first kind by fft based deconvolution with a gaussian low pass filter
 *
 * @param input
 * @param kernel
 * @param stab  Stabilisation parameter as a proportion of the max magnitude of the kernel. Fixes the lowpass total
 * cutoff point. For example, default 0.01 creates a low pass gaussian filter, between 0 and the point where the FFT
 * of the kernel reaches 1% of its max (usually DC) magnitude. Set <= 0 to circumvent low pass filtering.
 * @param flag
 * @return
 */
Eigen::VectorXd lowpass_fft_deconv(const Eigen::VectorXd &input, const Eigen::VectorXd &kernel, const std::string &flag, double stab=0.01) {

    // Map the input signature data to tensors (shared memory)
    auto input_len = input.rows();
    auto kernel_len = kernel.rows();

    // Compute cumulative length of input and kernel;
    auto N = input_len + kernel_len - 1;
    auto M = fft_next_good_size(N);

    // Create padded FFT inputs
    Eigen::FFT<double> fft;
    Eigen::VectorXd input_padded(M);
    Eigen::VectorXd kernel_padded(M);
    input_padded.setZero();
    kernel_padded.setZero();
    input_padded.topRows(input_len) = input;
    kernel_padded.topRows(kernel_len) = kernel;

    // Take the forward ffts
    Eigen::VectorXcd input_transformed(M);
    Eigen::VectorXcd kernel_transformed(M);
    fft.fwd(input_transformed, input_padded);
    fft.fwd(kernel_transformed, kernel_padded);

    // Magnitude of the transformed kernel and input
    Eigen::ArrayXd input_mag = pow(input_transformed.array().real(), 2.0) + pow(input_transformed.array().imag(), 2.0);
    input_mag = pow(input_mag, 0.5);
    Eigen::ArrayXd kernel_mag = pow(kernel_transformed.array().real(), 2.0) + pow(kernel_transformed.array().imag(), 2.0);
    kernel_mag = pow(kernel_mag, 0.5);

    // Double check plot
//    Figure fig = Figure();
//    ScatterPlot p3 = ScatterPlot();
//    p3.x = Eigen::ArrayXd::LinSpaced(M, 1, M);
//    p3.y = input_mag;
//    p3.name = "input mag";
//    fig.add(p3);
//    ScatterPlot p2 = ScatterPlot();
//    p2.x = Eigen::ArrayXd::LinSpaced(M, 1, M);
//    p2.y = kernel_mag;
//    p2.name = "kernel mag";
//    fig.add(p2);

    // Deconvolve by element-wise division, stabilising divide-by-0 errors based on the 1% of magnitude of the kernel
    Eigen::VectorXcd inter(M);
    inter.setZero();

    if (stab > 0) {
        double kernel_cutoff_magnitude = kernel_mag.maxCoeff() * stab;
        Eigen::Index ctr = 0;
        Eigen::Index location = -1;
        for (ctr = 0; ctr < M; ctr++) {
            if (kernel_mag(ctr) > kernel_cutoff_magnitude) {
                inter(ctr) = input_transformed(ctr) / kernel_transformed(ctr);
            } else {
                if (location == -1) {
                    location = ctr;
                };
//                break;
                inter(ctr) = 0;
            }
        }
        std::cout << "LOCATION " << location << std::endl;

        // The above works on its own, but produces high frequency ringing, because of the sharp cutoff low pass filter.
        // The ctr is now set at the point where we want the lowpass filter to fade out entirely, so determine a
        // gaussian form between position 0 and here

//        inter = input_transformed.array() / kernel_transformed.array();

        // TODO this can be done in-place to avoid duplicating memory but we want to plot for the time being
        Eigen::ArrayXd map_frequency = Eigen::ArrayXd::LinSpaced(location, -5.0, 5);
        Eigen::ArrayXd low_pass(M);
        low_pass.setZero();
        for (auto i = 0; i < map_frequency.rows(); i++) {
            low_pass(i) = 0.5 * (1.0 - std::erf(map_frequency(i)));
        }

        // We're working with a dual-sided FFT so mirror the filter
        low_pass = low_pass + low_pass.reverse();

        // Apply the low pass filter
        inter = inter.array() * low_pass;

        // Debug plot
//        ScatterPlot p5 = ScatterPlot();
//        p5.x = Eigen::ArrayXd::LinSpaced(low_pass.rows(), 1, low_pass.rows());
//        p5.y = low_pass;
//        p5.name = "lowpass magnitude";
//        fig.add(p5);

    } else {
        inter = input_transformed.array() / kernel_transformed.array();
    }

//    for (auto k = 0; k<M; k++) {
//        if ((k > 4) && (k < M-5)) {
//            inter(k) = 0.0;
//        }
//    }

//    std::cout << "Smoothing the inter" <<std::endl;
//    Eigen::ArrayXd inter_sm_imag(inter.size());
//    inter_sm_imag.segment(2, M-5) = (inter.imag().segment(0, M-5) + inter.imag().segment(1, M-5) + inter.imag().segment(2, M-5) + inter.imag().segment(3, M-5)  + inter.imag().segment(4, M-5)) / 5.0;
//    inter_sm_imag(1) = (inter.imag()(0) + inter.imag()(1) + inter.imag()(2))/3.0;
//    inter_sm_imag(M-2) = (inter.imag()(M-3) + inter.imag()(M-2) + inter.imag()(M-1))/3.0;
//    inter.imag() = inter_sm_imag;
//    Eigen::ArrayXd inter_sm_real(inter.size());
//    inter_sm_real.segment(2, M-5) = (inter.real().segment(0, M-5) + inter.real().segment(1, M-5) + inter.real().segment(2, M-5) + inter.real().segment(3, M-5)  + inter.real().segment(4, M-5)) / 5.0;
//    inter_sm_real(1) = (inter.real()(0) + inter.real()(1) + inter.real()(2))/3.0;
//    inter_sm_real(M-2) = (inter.real()(M-3) + inter.real()(M-2) + inter.real()(M-1))/3.0;
//    inter.real() = inter_sm_real;


//    Eigen::ArrayXd inter_mag = pow(inter.array().real(), 2.0) + pow(inter.array().imag(), 2.0);
//    inter_mag = pow(inter_mag, 0.5);
    Eigen::ArrayXd inter_mag = inter.array().imag();
//    ScatterPlot p4 = ScatterPlot();
//    p4.x = Eigen::ArrayXd::LinSpaced(M, 1, M);
//    p4.y = inter_mag;
//    p4.name = "inter mag";
//    fig.add(p4);
//
//    fig.write("check_that_fft_deconv_behaves_" + flag + ".json");
    // Inverse FFT
    Eigen::VectorXd out(M);
    out.setZero();
    fft.inv(out, inter);

    // Smooth the output
    std::cout << "Smoothing the out" <<std::endl;
    Eigen::ArrayXd inter_sm_imag(out.size());
    Eigen::ArrayXd out_sm(out.size());
    out_sm.segment(2, M-5) = (out.segment(0, M-5) + out.segment(1, M-5) + out.segment(2, M-5) + out.segment(3, M-5)  + out.segment(4, M-5)) / 5.0;
    out_sm(1) = (out(0) + out(1) + out(2))/3.0;
    out_sm(M-2) = (out(M-3) + out(M-2) + out(M-1))/3.0;
    out = out_sm;

    // Crop to the size of the input vector
    out = out.topRows(input_len);

    return out;
}

}  /* namespace utilities */



#endif //ES_FLOW_CONV_H
