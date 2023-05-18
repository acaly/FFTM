# FFTM
Fast Fourier Transform in Managed.

This library implements the FFT algorithm in pure C# code. It internally uses iterative SRFFT algorithm with SIMD acceleration.

### Benchmark

Below is the result of a simple benchmark of 64-bit complex 1D FFT with N = 1024, in comparison with FFTW.

| Method |    N |     Mean |   Error |  StdDev | Code Size |
|------- |----- |---------:|--------:|--------:|----------:|
|   FFTM | 1024 | 960.1 ns | 9.16 ns | 5.45 ns |   2,016 B |
|   FFTW | 1024 | 842.7 ns | 5.16 ns | 3.41 ns |      38 B |
