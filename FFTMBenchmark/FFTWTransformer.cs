using LibFFTM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace FFTMBenchmark
{
    internal unsafe partial class FFTWTransformer : AbstractTransformer<double, double>
    {
        private readonly int _n;
        private readonly void* _dataIn = null;
        private readonly void* _dataOut = null;
        private readonly IntPtr _plan;

        public override Span<double> Input => new(_dataIn, _n * 2);
        public override Span<double> Output => new(_dataOut, _n * 2);

        public FFTWTransformer(int n)
        {
            _dataIn = Malloc((nuint)n * 16);
            _dataOut = Malloc((nuint)n * 16);
            _n = n;
            _plan = PlanDft1D(n, _dataIn, _dataOut, -1 /* FFTW_FORWARD */, 0 /* FFTW_MEASURE */);
        }

        public override void Transform()
        {
            Execute(_plan);
        }

        public delegate void* MallocDelegate(UIntPtr size);
        public delegate void FreeDelegate(void* ptr);
        public delegate IntPtr PlanDft1DDelegate(int n0, void* pIn, void* pOut, int direction, uint flags);
        public delegate void DestroyPlanDelegate(IntPtr plan);
        public delegate void ExecuteDelegate(IntPtr plan);

        public static readonly MallocDelegate Malloc;
        public static readonly FreeDelegate Free;
        public static readonly PlanDft1DDelegate PlanDft1D;
        public static readonly DestroyPlanDelegate DestroyPlan;
        public static readonly ExecuteDelegate Execute;

        [LibraryImport("kernel32")]
        private static partial IntPtr LoadLibraryW([MarshalAs(UnmanagedType.LPWStr)] string dllName);
        [LibraryImport("kernel32")]
        private static partial IntPtr GetProcAddress(IntPtr m, [MarshalAs(UnmanagedType.LPStr)] string name);

        static FFTWTransformer()
        {
            var basePath = Path.GetDirectoryName(typeof(FFTWTransformer).Assembly.Location)!;
            var dllName = IntPtr.Size == 4 ? "fftw\\libfftw3d_x86.dll" : "fftw\\libfftw3d_x64.dll";
            var m = LoadLibraryW(Path.Combine(basePath, dllName));
            InitDelegate(out Malloc, m, "fftw_malloc");
            InitDelegate(out Free, m, "fftw_free");
            InitDelegate(out PlanDft1D, m, "fftw_plan_dft_1d");
            InitDelegate(out DestroyPlan, m, "fftw_destroy_plan");
            InitDelegate(out Execute, m, "fftw_execute");
        }

        private static void InitDelegate<T>(out T fieldRef, IntPtr m, string procName)
        {
            var addr = GetProcAddress(m, procName);
            fieldRef = Marshal.GetDelegateForFunctionPointer<T>(addr);
        }
    }
}
