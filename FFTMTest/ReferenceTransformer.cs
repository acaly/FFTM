using LibFFTM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace FFTMTest
{
    internal sealed class ReferenceTransformer : AbstractTransformer<double, double>
    {
        private readonly Complex[] _input;
        private readonly Complex[] _output;
        private readonly int[] _swap;

        public override Span<double> Input => MemoryMarshal.Cast<Complex, double>(_input);
        public override Span<double> Output => MemoryMarshal.Cast<Complex, double>(_output);

        public ReferenceTransformer(int n)
        {
            _input = new Complex[n];
            _output = new Complex[n];
            _swap = new int[n];
            MakeSwapArray(_swap, 0, 1, 0);
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

        public override void Transform()
        {
            Span<Complex> input = _input;
            Span<Complex> output = _output;
            var n = input.Length;
            for (int i = 0; i < n; ++i)
            {
                output[_swap[i]] = input[i];
            }
            for (var m2 = 1; m2 < n; m2 *= 2)
            {
                var (s, c) = Math.SinCos(2 * Math.PI / (m2 * 2));
                var wm = new Complex(c, s);
                for (var k = 0; k < n; k += m2 * 2)
                {
                    var w = new Complex(1, 0);
                    for (var j = 0; j < m2; ++j)
                    {
                        ref var refl = ref output[k + j];
                        ref var refr = ref output[k + j + m2];
                        var t = w * refr;
                        (refl, refr) = (refl + t, refl - t);
                        w *= wm;
                    }
                }
            }
        }
    }
}
