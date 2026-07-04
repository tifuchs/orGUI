#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <intrin.h>
#include <stdint.h>
#include <stdio.h>

static void cpuid(int out[4], int leaf, int subleaf)
{
    __cpuidex(out, leaf, subleaf);
}

static bool has_leaf(int leaf)
{
    int regs[4] = {0, 0, 0, 0};
    cpuid(regs, 0, 0);
    return regs[0] >= leaf;
}

static uint64_t xgetbv0()
{
    return _xgetbv(0);
}

int main()
{
    int regs[4] = {0, 0, 0, 0};
    cpuid(regs, 1, 0);

    const bool sse42 = (regs[2] & (1 << 20)) != 0;
    const bool avx_cpu = (regs[2] & (1 << 28)) != 0;
    const bool xsave = (regs[2] & (1 << 26)) != 0;
    const bool osxsave = (regs[2] & (1 << 27)) != 0;
    const bool avx_os = xsave && osxsave && ((xgetbv0() & 0x6) == 0x6);

    bool avx2 = false;
    bool avx512 = false;
    if (has_leaf(7)) {
        cpuid(regs, 7, 0);
        avx2 = (regs[1] & (1 << 5)) != 0;

        const bool avx512f = (regs[1] & (1 << 16)) != 0;
        const bool avx512dq = (regs[1] & (1 << 17)) != 0;
        const bool avx512cd = (regs[1] & (1 << 28)) != 0;
        const bool avx512bw = (regs[1] & (1 << 30)) != 0;
        const bool avx512vl = (regs[1] & (1 << 31)) != 0;
        const bool avx512_os = avx_os && ((xgetbv0() & 0xe6) == 0xe6);
        avx512 = avx512_os && avx512f && avx512dq && avx512cd && avx512bw && avx512vl;
    }

    if (avx512) {
        printf("AVX512 ");
    }
    if (avx_cpu && avx_os && avx2) {
        printf("AVX2 ");
    }
    if (avx_cpu && avx_os) {
        printf("AVX ");
    }
    if (sse42) {
        printf("SSE4.2 ");
    }
    puts("");
    return 0;
}
#else
int main()
{
    return 0;
}
#endif
