using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics.X86;
using System.Runtime.Intrinsics;
using System.Numerics;
using System.Runtime.InteropServices;

namespace LibFFTM.AvxImpl
{
    internal abstract unsafe class AbstractTransformerD1F64 : AbstractTransformer<double, double>
    {
        private static T* AllocAligned<T>(int count, int alignment, out T[] array) where T : unmanaged
        {
            var additional = alignment / Unsafe.SizeOf<T>();
            array = GC.AllocateArray<T>(count + additional, pinned: true);
            fixed (T* p = array)
            {
                var unalignedPtr = (nint)p;
                return (T*)((unalignedPtr + alignment - 1) & ~(alignment - 1));
            }
        }

        private static readonly double _w8v = Math.Sqrt(0.5);

        protected readonly int _n;
        private readonly double[] _buffers;
        protected readonly double* _input, _output, _table;

        public override Span<double> Input => new(_input, _n * 2);
        public override Span<double> Output => new(_output, _n * 2);

        protected readonly (int, int)[] _outputTable;
        private readonly double[] _tableAllBuffer;
        protected readonly double* _tableAll;

        public AbstractTransformerD1F64(int n, ReadOnlySpan<double> table)
        {
            _n = n;
            _input = AllocAligned(n * 5, 32, out _buffers);
            _output = &_input[n * 2];
            _table = &_input[n * 4];
            table.CopyTo(new(_table, n));

            //TODO move to procedure
            var outputTableFull = new int[n / 4];
            MakeSwapArray(outputTableFull, 0, 1, 0);
            _outputTable = Enumerable.Range(0, n / 4)
                .Select(i => (A: i, B: outputTableFull[i]))
                .Where(t => t.A < t.B)
                .ToArray();

            var tableAll = InitTableAll(n);
            _tableAll = AllocAligned(tableAll.Length * 2, 32, out _tableAllBuffer);
            tableAll.AsSpan().CopyTo(new Span<Complex>((Complex*)_tableAll, tableAll.Length));
        }

        private static Complex[] InitTableAll(int n)
        {
            List<Complex> ret = new();
            int[] permArray = new int[n];
            for (int i = 4; i <= n; i *= 4)
            {
                var (s, c) = Math.SinCos(Math.PI / 4 / i);
                Complex w = new(c, s);
                Complex ww = w * w;
                Complex wi0 = Complex.One;
                Complex wi1 = new(_w8v, _w8v);
                Complex wi2 = Complex.One;

                var index = ret.Count;
                for (int j = 0; j < i * 3; ++j)
                {
                    ret.Add(default);
                }

                MakeSwapArray(permArray.AsSpan()[..i], 0, 1, 0);

                for (int j = 0; j < i; ++j)
                {
                    var perm = permArray.AsSpan()[j];
                    ret[index + perm * 3 + 0] = wi0;
                    ret[index + perm * 3 + 1] = wi1;
                    ret[index + perm * 3 + 2] = wi2;
                    wi0 *= w;
                    wi1 *= w;
                    wi2 *= ww;
                }
            }
            return ret.ToArray();
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

        //The initial radix-4 SRFFT pass (A-AC/C-AB/B-CB).
        //Input layout is different from the normal pass. No outer loop.
        protected static void Iteration_R4V2_I(nint n4, nint len4, nint count, double* input, double* output)
        {
            Debug.Assert(len4 == 1);

            var countB = count / 6;
            var countA = countB + 1;
            var halfCount = count >> 1;
            var countC = halfCount - countA - countB;

            var n8 = n4 >> 1;
            n4 <<= 2;
            n8 <<= 2;

            double* inputPtr = input;
            double* outputPtr = output;
            Vector256<double> wi1 = Vector256.Create(_w8v, _w8v, _w8v, _w8v);
            Vector256<double> signMulI = Vector256.Create(-0.0, 0.0, -0.0, 0.0);

            //A-AC.
            for (nint j = 0; j < countA; ++j)
            {
                var p = Vector256.LoadAligned(&inputPtr[0]);
                var q = Vector256.LoadAligned(&inputPtr[n4]);
                var u = Vector256.LoadAligned(&inputPtr[n8]);
                var v = Vector256.LoadAligned(&inputPtr[n4 + n8]);
                inputPtr += 4;

                var (a0, a1) = TypeA(p, q);
                var (c0, c1) = TypeCI(u, v, signMulI);
                var (o0, o2) = TypeA(a0, c0);
                var (o1, o3) = TypeA(a1, c1);

                o0.StoreAligned(&outputPtr[0]);
                o1.StoreAligned(&outputPtr[n4]);
                o2.StoreAligned(&outputPtr[n8]);
                o3.StoreAligned(&outputPtr[n4 + n8]);
                outputPtr += 4;
            }

            //C-AB.
            var wi10 = Avx.Permute(wi1, 0b0000);
            var wi11 = Avx.Permute(wi1, 0b1111);
            for (nint j = 0; j < countC; ++j)
            {
                var p = Vector256.LoadAligned(&inputPtr[0]);
                var q = Vector256.LoadAligned(&inputPtr[n4]);
                var u = Vector256.LoadAligned(&inputPtr[n8]);
                var v = Vector256.LoadAligned(&inputPtr[n4 + n8]);
                inputPtr += 4;

                var (a0, a1) = TypeA(p, q);
                var (b0, b1) = TypeB(u, v);
                var (o0, o2) = TypeCI(a0, b0, signMulI);
                var (o1, o3) = TypeC(a1, b1, signMulI, wi10, wi11);

                o0.StoreAligned(&outputPtr[0]);
                o1.StoreAligned(&outputPtr[n4]);
                o2.StoreAligned(&outputPtr[n8]);
                o3.StoreAligned(&outputPtr[n4 + n8]);
                outputPtr += 4;
            }

            //B-CB.
            for (nint j = 0; j < countB; ++j)
            {
                var p = Vector256.LoadAligned(&inputPtr[0]);
                var q = Vector256.LoadAligned(&inputPtr[n4]);
                var u = Vector256.LoadAligned(&inputPtr[n8]);
                var v = Vector256.LoadAligned(&inputPtr[n4 + n8]);
                inputPtr += 4;

                var (c0, c1) = TypeCI(p, q, signMulI);
                var (b0, b1) = TypeB(u, v);
                var (o0, o2) = TypeB(c0, b0);
                var (o1, o3) = TypeB(c1, b1);

                o0.StoreAligned(&outputPtr[0]);
                o1.StoreAligned(&outputPtr[n4]);
                o2.StoreAligned(&outputPtr[n8]);
                o3.StoreAligned(&outputPtr[n4 + n8]);
                outputPtr += 4;
            }
        }

        //The normal radix-4 SRFFT pass (A-AC/C-AB/B-CB).
        protected static void Iteration_R4V2(nint n4, nint len4, nint count, double* input, ref double* table)
        {
            var countB = count / 6;
            var countA = countB + 1;
            var halfCount = count >> 1;
            var countC = halfCount - countA - countB;

            var n8 = n4 >> 1;
            count <<= 1;
            n4 <<= 2;
            n8 <<= 2;

            double* inputPtr = input;
            Vector256<double> signMulI = Vector256.Create(-0.0, 0.0, -0.0, 0.0);

            for (nint i = 0; i < len4; ++i)
            {
                Vector128<double> wi0h = Vector128.Load(table);
                Vector128<double> wi1h = Vector128.Load(table + 2);
                Vector128<double> wi2h = Vector128.Load(table + 4);
                Vector256<double> wi0 = Vector256.Create(wi0h, wi0h);
                Vector256<double> wi1 = Vector256.Create(wi1h, wi1h);
                Vector256<double> wi2 = Vector256.Create(wi2h, wi2h);
                table += 6;

                //A-AC.
                var wi20 = Avx.Permute(wi2, 0b0000);
                var wi21 = Avx.Permute(wi2, 0b1111);
                for (nint j = 0; j < countA; ++j)
                {
                    var p = Vector256.LoadAligned(&inputPtr[0]);
                    var q = Vector256.LoadAligned(&inputPtr[count * 2]);
                    var u = Vector256.LoadAligned(&inputPtr[count]);
                    var v = Vector256.LoadAligned(&inputPtr[count * 3]);

                    var (a0, a1) = TypeA(p, q);
                    var (c0, c1) = TypeC(u, v, signMulI, wi20, wi21);
                    var (o0, o2) = TypeA(a0, c0);
                    var (o1, o3) = TypeA(a1, c1);

                    o0.StoreAligned(&inputPtr[0]);
                    o1.StoreAligned(&inputPtr[count * 2]);
                    o2.StoreAligned(&inputPtr[count]);
                    o3.StoreAligned(&inputPtr[count * 3]);
                    inputPtr += 4;
                }

                //C-AB.
                var wi00 = Avx.Permute(wi0, 0b0000);
                var wi01 = Avx.Permute(wi0, 0b1111);
                var wi10 = Avx.Permute(wi1, 0b0000);
                var wi11 = Avx.Permute(wi1, 0b1111);
                for (nint j = 0; j < countC; ++j)
                {
                    var p = Vector256.LoadAligned(&inputPtr[0]);
                    var q = Vector256.LoadAligned(&inputPtr[count * 2]);
                    var u = Vector256.LoadAligned(&inputPtr[count]);
                    var v = Vector256.LoadAligned(&inputPtr[count * 3]);

                    var (a0, a1) = TypeA(p, q);
                    var (b0, b1) = TypeB(u, v);
                    var (o0, o2) = TypeC(a0, b0, signMulI, wi00, wi01);
                    var (o1, o3) = TypeC(a1, b1, signMulI, wi10, wi11);

                    o0.StoreAligned(&inputPtr[0]);
                    o1.StoreAligned(&inputPtr[count * 2]);
                    o2.StoreAligned(&inputPtr[count]);
                    o3.StoreAligned(&inputPtr[count * 3]);
                    inputPtr += 4;
                }

                //B-CB.
                for (nint j = 0; j < countB; ++j)
                {
                    var p = Vector256.LoadAligned(&inputPtr[0]);
                    var q = Vector256.LoadAligned(&inputPtr[count * 2]);
                    var u = Vector256.LoadAligned(&inputPtr[count]);
                    var v = Vector256.LoadAligned(&inputPtr[count * 3]);

                    var (c0, c1) = TypeC(p, q, signMulI, wi20, wi21);
                    var (b0, b1) = TypeB(u, v);
                    var (o0, o2) = TypeB(c0, b0);
                    var (o1, o3) = TypeB(c1, b1);

                    o0.StoreAligned(&inputPtr[0]);
                    o1.StoreAligned(&inputPtr[count * 2]);
                    o2.StoreAligned(&inputPtr[count]);
                    o3.StoreAligned(&inputPtr[count * 3]);
                    inputPtr += 4;
                }

                inputPtr += count * 3;
            }
        }

        //The last radix-2 SRFFT pass (A) followed by the final vectorization 2x pass.
        protected static void Iteration_R2V2_F(nint n4, nint len4, nint count, double* input, double* table)
        {
            Debug.Assert(count == 1);
            n4 <<= 1;
            var len8 = len4 >> 1;
            var pwi = table;
            for (nint j = 0; j < len8; ++j)
            {
                var p = Vector256.LoadAligned(input);
                var q = Vector256.LoadAligned(input + 4);
                var u = Vector256.LoadAligned(input + n4 * 2);
                var v = Vector256.LoadAligned(input + n4 * 2 + 4);

                var (x0, x1) = (p + q, p - q);
                var (y0, y1) = (u + v, u - v);

                var wi = Vector256.LoadAligned(pwi);
                var w0 = Avx.Permute(wi, 0b0000);
                var w1 = Avx.Permute(wi, 0b1111);
                pwi += 4;

                var xy00 = Avx.Permute2x128(x0, y0, 0b00100000);
                var xy01 = Avx.Permute2x128(x0, y0, 0b00110001);
                var xy10 = Avx.Permute2x128(x1, y1, 0b00100000);
                var xy11 = Avx.Permute2x128(x1, y1, 0b00110001);

                var xy01m = MulT(xy01, w0, w1);
                var xy11m = Mul(xy11, w1, w0);

                (xy00 + xy01m).StoreAligned(input);
                (xy00 - xy01m).StoreAligned(input + n4 * 2);
                (xy10 + xy11m).StoreAligned(input + 4);
                (xy10 - xy11m).StoreAligned(input + n4 * 2 + 4);
                input += 8;
            }
        }

        //The final vectorization 2x pass.
        protected static void Iteration_F(nint n4, nint len4, nint count, double* input, double* table)
        {
            Debug.Assert(count == 0);
            var n2 = n4 << 1;
            n2 <<= 1;
            var len16 = len4 >> 2;
            var pwi = table;
            for (nint j = 0; j < len16; ++j)
            {
                var p = Vector256.LoadAligned(input);
                var q = Vector256.LoadAligned(input + n2);
                var u = Vector256.LoadAligned(input + 4);
                var v = Vector256.LoadAligned(input + n2 + 4);

                var wi = Vector256.LoadAligned(pwi);
                pwi += 4;
                var w0 = Avx.Permute(wi, 0b0000);
                var w1 = Avx.Permute(wi, 0b1111);

                var pq0 = Avx.Permute2x128(p, q, 0b00100000);
                var pq1 = Avx.Permute2x128(p, q, 0b00110001);
                var uv0 = Avx.Permute2x128(u, v, 0b00100000);
                var uv1 = Avx.Permute2x128(u, v, 0b00110001);
                var pq1m = MulT(pq1, w0, w1);
                var uv1m = Mul(uv1, w1, w0);

                (pq0 + pq1m).StoreAligned(input);
                (pq0 - pq1m).StoreAligned(input + n2);
                (uv0 + uv1m).StoreAligned(input + 4);
                (uv0 - uv1m).StoreAligned(input + n2 + 4);
                input += 8;
            }
        }

        protected static void Iteration_Swap(nint n4, double* input, (int, int)[] outputTable)
        {
            var table = outputTable;
            var input1 = input;
            var input2 = input + (n4 << 2);
            for (nint i = 0; i < table.Length; ++i)
            {
                var (a, b) = table[i];
                var va = Vector256.LoadAligned(&input1[a * 4]);
                var vb = Vector256.LoadAligned(&input1[b * 4]);
                va.StoreAligned(&input1[b * 4]);
                vb.StoreAligned(&input1[a * 4]);
                var vc = Vector256.LoadAligned(&input2[a * 4]);
                var vd = Vector256.LoadAligned(&input2[b * 4]);
                vc.StoreAligned(&input2[b * 4]);
                vd.StoreAligned(&input2[a * 4]);
            }
        }

        //A type-A radix-2 butterfly unit.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static (Vector256<double>, Vector256<double>) TypeA(Vector256<double> a, Vector256<double> c)
        {
            return (a + c, a - c);
        }

        //A type-B radix-2 butterfly unit.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static (Vector256<double>, Vector256<double>) TypeB(Vector256<double> c, Vector256<double> b)
        {
            return (b + c, b - c);
        }

        //A type-C radix-2 butterfly unit with wi = 1.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static (Vector256<double>, Vector256<double>) TypeCI(Vector256<double> a, Vector256<double> b, Vector256<double> sign)
        {
            return (a + b, Avx.Permute(a - b, 0b0101) ^ sign);
        }

        //A type-C radix-2 butterfly unit.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static (Vector256<double>, Vector256<double>) TypeC(Vector256<double> a, Vector256<double> b,
            Vector256<double> sign, Vector256<double> w0, Vector256<double> w1)
        {
            var aw = Mul(a, w0, w1);
            var bw = MulT(b, w0, w1);
            return (aw + bw, Avx.Permute(aw - bw, 0b0101) ^ sign);
        }

        //Complex multiplication, vectorized 2x.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Vector256<double> Mul(Vector256<double> a, Vector256<double> w)
        {
            return Mul(a, Avx.Permute(w, 0b0000), Avx.Permute(w, 0b1111));
        }

        //Complex multiplication, vectorized 2x, with one multiplier split into real and imaginary parts.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Vector256<double> Mul(Vector256<double> a, Vector256<double> w0, Vector256<double> w1)
        {
            return Fma.MultiplyAddSubtract(a, w0, Avx.Permute(a, 0b0101) * w1);
        }

        //Complex multiplication, vectorized 2x, with one multiplier split into real and imaginary parts
        //and the imaginary part negated.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Vector256<double> MulT(Vector256<double> a, Vector256<double> w0, Vector256<double> w1)
        {
            return Fma.MultiplySubtractAdd(a, w0, Avx.Permute(a, 0b0101) * w1);
        }
    }
}
