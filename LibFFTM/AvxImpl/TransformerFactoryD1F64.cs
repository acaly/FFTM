using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace LibFFTM.AvxImpl
{
    internal sealed class TransformerFactoryD1F64 : TransformerFactoryD1<double, double>
    {
        private sealed class TwoPowerProcedure : AbstractTransformerProcedure<double, double>
        {
            private readonly int _n;
            private readonly int _bits;
            private readonly double[] _table;

            public TwoPowerProcedure(int n)
            {
                _n = n;
                _bits = 32 - BitOperations.LeadingZeroCount((uint)n / 16);
                _table = CreateTable(n);
            }

            private static double[] CreateTable(int n)
            {
                var (im, re) = Math.SinCos(Math.PI * 2 / n);
                Complex w = new(re, -im);
                Complex wi = Complex.One;
                var ret = new double[n / 2];

                int[] permArray = new int[n / 8];
                MakeSwapArray(permArray.AsSpan(), 0, 1, 0);

                for (int j = 0; j < n / 8; ++j)
                {
                    var s = permArray[j];
                    ret[s * 4 + 0] = wi.Real;
                    ret[s * 4 + 1] = wi.Imaginary;
                    wi *= w;
                    ret[s * 4 + 2] = wi.Real;
                    ret[s * 4 + 3] = wi.Imaginary;
                    wi *= w;
                }
                return ret;
            }

            private static void MakeSwapArray(Span<int> span, int offset, int interval, int index)
            {
                if (interval == span.Length)
                {
                    span[offset] = index;
                    return;
                }
                MakeSwapArray(span, offset, interval * 2, index);
                MakeSwapArray(span, offset + interval, interval * 2, index + span.Length / 2 / interval);
            }

            public override AbstractTransformer<double, double> CreateTransformer()
            {
                switch (_bits & 3)
                {
                case 1:
                    return new TransformerD1F64M1(_n, _table);
                case 2:
                    return new TransformerD1F64M2(_n, _table);
                case 3:
                    return new TransformerD1F64M3(_n, _table);
                case 0:
                    return new TransformerD1F64M4(_n, _table);
                default:
                    return null!;
                }
            }
        }

        public override AbstractTransformerProcedure<double, double> CreateProcedure(int n)
        {
            if (n < 16)
            {
                //TODO
                throw new NotImplementedException();
            }
            if ((n & (n - 1)) != 0)
            {
                //TODO
                throw new NotImplementedException();
            }
            return new TwoPowerProcedure(n);
        }
    }
}
