using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace LibFFTM
{
    public abstract class AbstractTransformer
    {
        public abstract void Transform();
    }

    public abstract class AbstractTransformer<TIn, TOut> : AbstractTransformer
        where TIn : struct
        where TOut : struct
    {
        public abstract Span<TIn> Input { get; }
        public abstract Span<TOut> Output { get; }
    }
}
