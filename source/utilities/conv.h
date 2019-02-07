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


namespace utilities {


/** @brief Find the next good size for an fft, to pad with minimal number of zeros
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
    while (TRUE) {
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

/** @brief Convolve an input vector with a kernel, returning an output of the same length as the input
 *
 * Uses zero-padded fft based output
 *
 * @param out
 * @param input
 * @param kernel
 */
Eigen::VectorXd conv(const Eigen::VectorXd &input, const Eigen::VectorXd &kernel) {

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

    // Inverse FFT to deconvolve
    Eigen::VectorXd out(M);
    out.setZero();
    fft.inv(out, inter);

    // Crop to the size of the input vector
    out = out.topRows(input_len);
    return out;
}

}  /* namespace utilities */

#endif //ES_FLOW_CONV_H
