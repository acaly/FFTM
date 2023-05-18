using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace LibFFTM.AvxImpl
{
    internal unsafe sealed class TransformerD1F64M1 : AbstractTransformerD1F64
    {
        public TransformerD1F64M1(int n, ReadOnlySpan<double> table) : base(n, table)
        {
        }

        public override void Transform()
        {
            Debug.Assert(_n >= 16);
            nint n4 = _n >> 2;
            var inputPtr = _input;
            var outputPtr = _output;
            var tempPtr = _temp;

            nint len4 = 1;
            nint count = n4;

            Iteration_R4V2_I(n4, len4, count, inputPtr, outputPtr);
            len4 <<= 2;
            count >>= 2;

            double* table = _tableAll;
            while (count >= 16)
            {
                Iteration_R4V2(n4, len4, count, outputPtr, ref table);
                Iteration_R4V2(n4, len4 << 2, count >> 2, outputPtr, ref table);
                len4 <<= 4;
                count >>= 4;
            }

            Debug.Assert(count == 1);
            Iteration_R2V2_F(n4, len4, count, outputPtr, _table);
            Iteration_Swap(n4, outputPtr, _outputTable);
        }
    }

    internal unsafe sealed class TransformerD1F64M2 : AbstractTransformerD1F64
    {
        public TransformerD1F64M2(int n, ReadOnlySpan<double> table) : base(n, table)
        {
        }

        public override void Transform()
        {
            Debug.Assert(_n >= 16);
            nint n4 = _n >> 2;
            var inputPtr = _input;
            var outputPtr = _output;
            var tempPtr = _temp;

            nint len4 = 1;
            nint count = n4;

            Iteration_R4V2_I(n4, len4, count, inputPtr, outputPtr);
            len4 <<= 2;
            count >>= 2;

            double* table = _tableAll;
            while (count >= 16)
            {
                Iteration_R4V2(n4, len4, count, outputPtr, ref table);
                Iteration_R4V2(n4, len4 << 2, count >> 2, outputPtr, ref table);
                len4 <<= 4;
                count >>= 4;
            }

            Debug.Assert(count == 2);
            Iteration_R4V2(n4, len4, count, outputPtr, ref table);
            Iteration_F(n4, len4 << 2, count >> 2, outputPtr, _table);
            Iteration_Swap(n4, outputPtr, _outputTable);
        }
    }

    internal unsafe sealed class TransformerD1F64M3 : AbstractTransformerD1F64
    {
        public TransformerD1F64M3(int n, ReadOnlySpan<double> table) : base(n, table)
        {
        }

        public override void Transform()
        {
            Debug.Assert(_n >= 16);
            nint n4 = _n >> 2;
            var inputPtr = _input;
            var outputPtr = _output;
            var tempPtr = _temp;

            nint len4 = 1;
            nint count = n4;

            Iteration_R4V2_I(n4, len4, count, inputPtr, outputPtr);
            len4 <<= 2;
            count >>= 2;

            double* table = _tableAll;
            while (count >= 16)
            {
                Iteration_R4V2(n4, len4, count, outputPtr, ref table);
                Iteration_R4V2(n4, len4 << 2, count >> 2, outputPtr, ref table);
                len4 <<= 4;
                count >>= 4;
            }

            Debug.Assert(count == 4);
            Iteration_R4V2(n4, len4, count, outputPtr, ref table);
            Iteration_R2V2_F(n4, len4 << 2, count >> 2, outputPtr, _table);
            Iteration_Swap(n4, outputPtr, _outputTable);
        }
    }

    internal unsafe sealed class TransformerD1F64M4 : AbstractTransformerD1F64
    {
        public TransformerD1F64M4(int n, ReadOnlySpan<double> table) : base(n, table)
        {
        }

        public override void Transform()
        {
            Debug.Assert(_n >= 16);
            nint n4 = _n >> 2;
            var inputPtr = _input;
            var outputPtr = _output;
            var tempPtr = _temp;

            nint len4 = 1;
            nint count = n4;

            Iteration_R4V2_I(n4, len4, count, inputPtr, outputPtr);
            len4 <<= 2;
            count >>= 2;

            double* table = _tableAll;
            while (count >= 16)
            {
                Iteration_R4V2(n4, len4, count, outputPtr, ref table);
                Iteration_R4V2(n4, len4 << 2, count >> 2, outputPtr, ref table);
                len4 <<= 4;
                count >>= 4;
            }

            Debug.Assert(count == 8);
            Iteration_R4V2(n4, len4, count, outputPtr, ref table);
            Iteration_R4V2(n4, len4 << 2, count >> 2, outputPtr, ref table);
            Iteration_F(n4, len4 << 4, count >> 4, outputPtr, _table);
            Iteration_Swap(n4, outputPtr, _outputTable);
        }
    }
}
