# FFTM
Fast Fourier Transform in Managed.

This library implements the FFT algorithm in pure C# code. It internally uses iterative SRFFT algorithm with SIMD acceleration.

### Benchmark

Below is the result of a simple benchmark of 64-bit complex 1D FFT with N = 1024, in comparison with FFTW.

|       Method |       Mean |    Error |  StdDev | Code Size |
|------------- |-----------:|---------:|--------:|----------:|
|     F64N1024 | 1,181.7 ns | 11.19 ns | 7.40 ns |   2,802 B |
| FFTWF64N1024 |   843.3 ns |  6.63 ns | 4.38 ns |      38 B |

