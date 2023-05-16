using LibFFTM;
using LibFFTM.AvxImpl;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using Xunit;

namespace FFTMTest.AvxImpl
{
    using DoubleTransformer = AbstractTransformer<double, double>;

    public class TransformerFactoryD1F64Test
    {
        private static DoubleTransformer CreateTransformer(int n)
        {
            return new TransformerFactoryD1F64().CreateProcedure(n).CreateTransformer();
        }

        private static void Fill(DoubleTransformer t)
        {
            var span = t.Input;
            for (int i = 0; i < span.Length / 2; ++i)
            {
                span[i * 2 + 0] = i;
                span[i * 2 + 1] = 0;
            }
        }

        private static void CheckResult(int n, DoubleTransformer t)
        {
            var refTransformer = new ReferenceTransformer(n);
            t.Input.CopyTo(refTransformer.Input);
            refTransformer.Transform();

            for (int i = 0; i < n; ++i)
            {
                var v = new Complex(t.Output[i * 2], t.Output[i * 2 + 1]);
                var refv = new Complex(refTransformer.Output[i * 2], refTransformer.Output[i * 2 + 1]);
                if ((v - refv).Magnitude > 0.01 * (v + refv).Magnitude)
                {
                    Assert.Fail($"Result mismatch at i, expected = {refv}, actual = {v}");
                }
            }
        }

        [Fact]
        public void Test16()
        {
            var transformer = CreateTransformer(16);
            Fill(transformer);
            transformer.Transform();
            CheckResult(16, transformer);
        }

        [Fact]
        public void Test1024()
        {
            var transformer = CreateTransformer(1024);
            Fill(transformer);
            transformer.Transform();
            CheckResult(1024, transformer);
        }
    }
}
