using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.Intrinsics.X86;
using System.Text;
using System.Threading.Tasks;

namespace LibFFTM
{
    public abstract class TransformerFactoryD1<TIn, TOut>
        where TIn : struct
        where TOut : struct
    {
        public abstract AbstractTransformerProcedure<TIn, TOut> CreateProcedure(int n);
    }

    public static class TransformerFactory
    {
        public static readonly TransformerFactoryD1<double, double> ComplexF64 = CreateComplexF64();

        private static TransformerFactoryD1<double, double> CreateComplexF64()
        {
            if (Avx.IsSupported)
            {
                return new AvxImpl.TransformerFactoryD1F64();
            }
            throw new PlatformNotSupportedException();
        }
    }
}
