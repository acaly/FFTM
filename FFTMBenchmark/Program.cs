using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Diagnosers;
using BenchmarkDotNet.Running;
using LibFFTM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace FFTMBenchmark
{
    using DoubleTransformer = AbstractTransformer<double, double>;

    internal class Program
    {
        static void Main()
        {
            BenchmarkRunner.Run<Benchmark>();
        }
    }

    [DisassemblyDiagnoser(printInstructionAddresses: true, syntax: DisassemblySyntax.Masm, maxDepth: 3)]
    [SimpleJob(iterationCount: 10, warmupCount: 10)]
    public class Benchmark
    {
        private static readonly DoubleTransformer _transformerF64N1024 = CreateTransformerF64(1024);
        private static readonly FFTWTransformer _transformerFFTWF64N1024 = new(1024);

        private static DoubleTransformer CreateTransformerF64(int n)
        {
            return TransformerFactory.ComplexF64.CreateProcedure(n).CreateTransformer();
        }

        public static void WriteTestDataF64(DoubleTransformer t, double k, double mag)
        {
            var inputBuffer = MemoryMarshal.Cast<double, Complex>(t.Input);
            for (int i = 0; i < inputBuffer.Length; ++i)
            {
                var c = Math.Cos(i * k);
                inputBuffer[i] = new Complex(c * mag, 0);
            }
        }

        public Benchmark()
        {
            WriteTestDataF64(_transformerF64N1024, 0.2, 1);
            WriteTestDataF64(_transformerFFTWF64N1024, 0.2, 1);
        }

        #pragma warning disable CA1822 // Mark members as static

        [Benchmark]
        public void F64N1024()
        {
            _transformerF64N1024.Transform();
        }

        [Benchmark]
        public void FFTWF64N1024()
        {
            _transformerFFTWF64N1024.Transform();
        }

        #pragma warning restore CA1822 // Mark members as static
    }
}
