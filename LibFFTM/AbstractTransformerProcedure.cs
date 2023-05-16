using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace LibFFTM
{
    public abstract class AbstractTransformerProcedure<TIn, TOut>
        where TIn : struct
        where TOut : struct
    {
        public abstract AbstractTransformer<TIn, TOut> CreateTransformer();
    }
}
