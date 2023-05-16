using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics.X86;
using System.Runtime.Intrinsics;

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
        protected readonly double* _input, _output, _temp, _table;

        public override Span<double> Input => new(_input, _n * 2);
        public override Span<double> Output => new(_output, _n * 2);

        public AbstractTransformerD1F64(int n, ReadOnlySpan<double> table)
        {
            _n = n;
            _input = AllocAligned(n * 7, 32, out _buffers);
            _output = &_input[n * 2];
            _temp = &_input[n * 4];
            _table = &_input[n * 6];
            table.CopyTo(new(_table, n));
        }

        //The initial radix-4 SRFFT pass (A-AC/C-AB/B-CB).
        //Input layout is different from the normal pass. No outer loop.
        protected static void Iteration_R4V2_I(nint n4, nint len4, nint count, double* input, double* output)
        {
            Debug.Assert(len4 == 1);

            var countB = count / 6;
            var countA = countB + 1;
            var halfCount = count >> 2;
            var countC1 = halfCount - countA;
            var countC2 = halfCount - countB;

            var n8 = n4 >> 1;
            n4 <<= 2;
            n8 <<= 2;

            double* inputPtr = input;
            double* outputPtrStart = output;
            Vector256<double> wi1 = Vector256.Create(_w8v, _w8v, _w8v, _w8v);
            Vector256<double> signMulI = Vector256.Create(-0.0, 0.0, -0.0, 0.0);

            double* outputPtr = outputPtrStart;

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
                o2.StoreAligned(&outputPtr[n4]);
                o1.StoreAligned(&outputPtr[n8]);
                o3.StoreAligned(&outputPtr[n4 + n8]);
                outputPtr += 8;
            }

            //C-AB.
            for (nint j = 0; j < countC1; ++j)
            {
                var p = Vector256.LoadAligned(&inputPtr[0]);
                var q = Vector256.LoadAligned(&inputPtr[n4]);
                var u = Vector256.LoadAligned(&inputPtr[n8]);
                var v = Vector256.LoadAligned(&inputPtr[n4 + n8]);
                inputPtr += 4;

                var (a0, a1) = TypeA(p, q);
                var (b0, b1) = TypeB(u, v);
                var (o0, o2) = TypeCI(a0, b0, signMulI);
                var (o1, o3) = TypeC(a1, b1, signMulI, wi1);

                o0.StoreAligned(&outputPtr[0]);
                o2.StoreAligned(&outputPtr[n4]);
                o1.StoreAligned(&outputPtr[n8]);
                o3.StoreAligned(&outputPtr[n4 + n8]);
                outputPtr += 8;
            }

            outputPtr = outputPtrStart + 4;

            //C-AB.
            for (nint j = 0; j < countC2; ++j)
            {
                var p = Vector256.LoadAligned(&inputPtr[0]);
                var q = Vector256.LoadAligned(&inputPtr[n4]);
                var u = Vector256.LoadAligned(&inputPtr[n8]);
                var v = Vector256.LoadAligned(&inputPtr[n4 + n8]);
                inputPtr += 4;

                var (a0, a1) = TypeA(p, q);
                var (b0, b1) = TypeB(u, v);
                var (o0, o2) = TypeCI(a0, b0, signMulI);
                var (o1, o3) = TypeC(a1, b1, signMulI, wi1);

                o0.StoreAligned(&outputPtr[0]);
                o2.StoreAligned(&outputPtr[n4]);
                o1.StoreAligned(&outputPtr[n8]);
                o3.StoreAligned(&outputPtr[n4 + n8]);
                outputPtr += 8;
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
                o2.StoreAligned(&outputPtr[n4]);
                o1.StoreAligned(&outputPtr[n8]);
                o3.StoreAligned(&outputPtr[n4 + n8]);
                outputPtr += 8;
            }
        }

        //The normal radix-4 SRFFT pass (A-AC/C-AB/B-CB).
        protected static void Iteration_R4V2(nint n4, nint len4, nint count, double* input, double* output)
        {
            var countB = count / 6;
            var countA = countB + 1;
            var halfCount = count >> 2;
            var countC1 = halfCount - countA;
            var countC2 = halfCount - countB;

            var n8 = n4 >> 1;
            count <<= 2;
            n4 <<= 2;
            n8 <<= 2;

            double* inputPtr = input;
            double* outputPtrStart = output;
            var (s, c) = Math.SinCos(Math.PI / 4 / len4);
            Vector256<double> w = Vector256.Create(c, s, c, s);
            Vector256<double> ww = Mul(w, w);
            Vector256<double> wi0 = Vector256.Create(1.0, 0.0, 1.0, 0.0);
            Vector256<double> wi1 = Vector256.Create(_w8v, _w8v, _w8v, _w8v);
            Vector256<double> wi2 = Vector256.Create(1.0, 0.0, 1.0, 0.0);
            Vector256<double> signMulI = Vector256.Create(-0.0, 0.0, -0.0, 0.0);

            for (nint i = 0; i < len4; ++i)
            {
                double* outputPtr = outputPtrStart;

                //A-AC.
                for (nint j = 0; j < countA; ++j)
                {
                    var p = Vector256.LoadAligned(&inputPtr[0]);
                    var q = Vector256.LoadAligned(&inputPtr[4]);
                    var u = Vector256.LoadAligned(&inputPtr[count]);
                    var v = Vector256.LoadAligned(&inputPtr[count + 4]);
                    inputPtr += 8;

                    var (a0, a1) = TypeA(p, q);
                    var (c0, c1) = TypeC(u, v, signMulI, wi2);
                    var (o0, o2) = TypeA(a0, c0);
                    var (o1, o3) = TypeA(a1, c1);

                    o0.StoreAligned(&outputPtr[0]);
                    o2.StoreAligned(&outputPtr[n4]);
                    o1.StoreAligned(&outputPtr[n8]);
                    o3.StoreAligned(&outputPtr[n4 + n8]);
                    outputPtr += 8;
                }

                //C-AB.
                for (nint j = 0; j < countC1; ++j)
                {
                    var p = Vector256.LoadAligned(&inputPtr[0]);
                    var q = Vector256.LoadAligned(&inputPtr[4]);
                    var u = Vector256.LoadAligned(&inputPtr[count]);
                    var v = Vector256.LoadAligned(&inputPtr[count + 4]);
                    inputPtr += 8;

                    var (a0, a1) = TypeA(p, q);
                    var (b0, b1) = TypeB(u, v);
                    var (o0, o2) = TypeC(a0, b0, signMulI, wi0);
                    var (o1, o3) = TypeC(a1, b1, signMulI, wi1);

                    o0.StoreAligned(&outputPtr[0]);
                    o2.StoreAligned(&outputPtr[n4]);
                    o1.StoreAligned(&outputPtr[n8]);
                    o3.StoreAligned(&outputPtr[n4 + n8]);
                    outputPtr += 8;
                }

                outputPtr = outputPtrStart + 4;

                //C-AB.
                for (nint j = 0; j < countC2; ++j)
                {
                    var p = Vector256.LoadAligned(&inputPtr[0]);
                    var q = Vector256.LoadAligned(&inputPtr[4]);
                    var u = Vector256.LoadAligned(&inputPtr[count]);
                    var v = Vector256.LoadAligned(&inputPtr[count + 4]);
                    inputPtr += 8;

                    var (a0, a1) = TypeA(p, q);
                    var (b0, b1) = TypeB(u, v);
                    var (o0, o2) = TypeC(a0, b0, signMulI, wi0);
                    var (o1, o3) = TypeC(a1, b1, signMulI, wi1);

                    o0.StoreAligned(&outputPtr[0]);
                    o2.StoreAligned(&outputPtr[n4]);
                    o1.StoreAligned(&outputPtr[n8]);
                    o3.StoreAligned(&outputPtr[n4 + n8]);
                    outputPtr += 8;
                }

                //B-CB.
                for (nint j = 0; j < countB; ++j)
                {
                    var p = Vector256.LoadAligned(&inputPtr[0]);
                    var q = Vector256.LoadAligned(&inputPtr[4]);
                    var u = Vector256.LoadAligned(&inputPtr[count]);
                    var v = Vector256.LoadAligned(&inputPtr[count + 4]);
                    inputPtr += 8;

                    var (c0, c1) = TypeC(p, q, signMulI, wi2);
                    var (b0, b1) = TypeB(u, v);
                    var (o0, o2) = TypeB(c0, b0);
                    var (o1, o3) = TypeB(c1, b1);

                    o0.StoreAligned(&outputPtr[0]);
                    o2.StoreAligned(&outputPtr[n4]);
                    o1.StoreAligned(&outputPtr[n8]);
                    o3.StoreAligned(&outputPtr[n4 + n8]);
                    outputPtr += 8;
                }

                outputPtrStart += count >> 1;
                inputPtr += count;
                wi0 = Mul(wi0, w);
                wi1 = Mul(wi1, w);
                wi2 = Mul(wi2, ww);
            }
        }

        //The last radix-2 SRFFT pass (A) followed by the final vectorization 2x pass.
        protected static void Iteration_R2V2_F(nint n4, nint len4, nint count, double* input, double* output, double* table)
        {
            Debug.Assert(count == 1);
            n4 <<= 1;
            var len8 = len4 >> 1;
            var pwi = table;
            for (nint j = 0; j < len8; ++j)
            {
                var p = Vector256.LoadAligned(input);
                var q = Vector256.LoadAligned(input + 4);
                var u = Vector256.LoadAligned(input + 8);
                var v = Vector256.LoadAligned(input + 12);
                input += 16;

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

                (xy00 + xy01m).StoreAligned(&output[0]);
                (xy10 + xy11m).StoreAligned(&output[n4]);
                (xy00 - xy01m).StoreAligned(&output[n4 * 2]);
                (xy10 - xy11m).StoreAligned(&output[n4 * 3]);
                output += 4;
            }
        }

        //The final vectorization 2x pass.
        protected static void Iteration_F(nint n4, nint len4, nint count, double* input, double* output, double* table)
        {
            Debug.Assert(count == 0);
            var n2 = n4 << 1;
            n2 <<= 1;
            n4 <<= 1;
            var len16 = len4 >> 2;
            var pwi = table;
            for (nint j = 0; j < len16; ++j)
            {
                var p = Vector256.LoadAligned(input);
                var q = Vector256.LoadAligned(input + 4);
                var u = Vector256.LoadAligned(input + n2);
                var v = Vector256.LoadAligned(input + n2 + 4);
                input += 8;

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

                (pq0 + pq1m).StoreAligned(&output[0]);
                (pq0 - pq1m).StoreAligned(&output[n2]);
                (uv0 + uv1m).StoreAligned(&output[n4]);
                (uv0 - uv1m).StoreAligned(&output[n2 + n4]);
                output += 4;
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
            Vector256<double> sign, Vector256<double> w)
        {
            var w0 = Avx.Permute(w, 0b0000);
            var w1 = Avx.Permute(w, 0b1111);
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
