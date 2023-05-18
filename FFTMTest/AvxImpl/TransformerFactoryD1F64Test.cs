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

        [Theory]
        [InlineData(16)]
        [InlineData(32)]
        [InlineData(64)]
        [InlineData(128)]
        [InlineData(256)]
        [InlineData(512)]
        [InlineData(1024)]
        public void TestAny(int n)
        {
            var transformer = CreateTransformer(n);
            Fill(transformer);
            transformer.Transform();
            CheckResult(n, transformer);
        }
    }
}
