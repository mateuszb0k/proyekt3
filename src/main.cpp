#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include <vector>
#include <cmath>
#include <matplot/matplot.h>
#include <set>
#include <sndfile.hh>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace py = pybind11;
namespace m = matplot;
typedef std::complex<double> Complex;
void sinus(int freq, int num_samples = 1000) {
    std::vector<double> x(num_samples);
    std::vector<double> y(num_samples);
    double step = 2 * 3.14159265359 * freq / num_samples;
    for (int i = 0; i < num_samples; ++i) {
        x[i] = i * step;
        y[i] = sin(i * step);
    }
    m::plot(x, y);
    m::show();
}

void cosinus(int freq, int num_samples = 1000) {
    std::vector<double> x(num_samples);
    std::vector<double> y(num_samples);
    double step = 2 * 3.14159265359 * freq / num_samples;
    for (int i = 0; i < num_samples; ++i) {
        x[i] = i * step;
        y[i] = cos(i * step);
    }
    m::plot(x, y);
    m::show();
}

void pila(int freq, int num_samples = 1000) {
    std::vector<double> x(num_samples);
    std::vector<double> y(num_samples);
    double step = 2 * 3.14159265359 * freq / num_samples;
    for (int i = 0; i < num_samples; ++i) {
        x[i] = i * step;
        y[i] = std::fmod(i * step, 2 * 3.14159265359);
    }
    m::plot(x, y);
    m::show();
}

void prostokatny(int freq, int num_samples = 1000) { //default 1000 najoptymalniejsze rezultaty?
    std::vector<double> x(num_samples);
    std::vector<double> y(num_samples);
    double step = 2 * 3.14159265359 * freq / num_samples;
    for (int i = 0; i < num_samples; ++i) {
        x[i] = i * step;
        y[i] = sin(i * step) > 0 ? 1 : -1;
    }
    m::plot(x, y);
    m::show();
}
void visualize_wav(const std::string& filename, int downsample_factor = 10) { //optymalny czas/jakosc downsample factor ok 40!
    SndfileHandle file { filename };
    if (file.error()) {
        std::cerr << "Error: blad pliku\n";
        return;
    }
    if (file.channels() > 2) {
        std::cerr << "Error: tylko pliki mono i stereo\n";
        return;
    }
    size_t n_frames = file.frames();
    std::vector<double> data(n_frames * file.channels());
    file.read(data.data(), n_frames);
    std::vector<double> time(n_frames);
    std::vector<double> downsampled_data;
    std::vector<double> downsampled_time;

    for (size_t i = 0; i < n_frames; i += downsample_factor) {
        downsampled_time.push_back(static_cast<double>(i) / file.samplerate());
        downsampled_data.push_back(data[i]);
    }
    matplot::plot(downsampled_time, downsampled_data);
    matplot::show();
}
std::vector<Complex> DFT(const std::vector<double>& input) {
    int N = input.size();
    std::vector<Complex> output(N);

    // Dla każdej wartości wyjściowej
    for (int k = 0; k < N; ++k) {
        // Sumujemy iloczyn wartości wejściowych i eksponencjalnej zespolonej
        for (int n = 0; n < N; ++n) {
            double theta = 2.0 * M_PI * k * n / N;
            output[k] += input[n] * Complex(cos(theta), -sin(theta));
        }
    }

    return output;
}

// Funkcja do wizualizacji DFT
void visualize_DFT(const std::vector<double>& input) {
    // Obliczamy DFT sygnału wejściowego
    std::vector<Complex> output = DFT(input);

    // Tworzymy dwa wektory: częstotliwości i amplitudy
    std::vector<double> frequencies(output.size());
    std::vector<double> magnitudes(output.size());

    // Wypełniamy wektory wartościami
    for (size_t i = 0; i < output.size(); ++i) {
        frequencies[i] = i;
        magnitudes[i] = std::abs(output[i]);
    }

    // Rysujemy wykres amplitud w funkcji częstotliwości
    matplot::plot(frequencies, magnitudes);
    matplot::xlabel("Częstotliwość");
    matplot::ylabel("Amplituda");
    matplot::show();
}


PYBIND11_MODULE(cmake_example, m) {//nagłówki funkcji dla cmake 

    m.def("sinus", &sinus);
    m.def("cosinus", &cosinus);
    m.def("pila", &pila);
    m.def("prostokatny", &prostokatny);
    m.def("visualize_wav", &visualize_wav);
    m.def("visualize_DFT", &visualize_DFT);

}