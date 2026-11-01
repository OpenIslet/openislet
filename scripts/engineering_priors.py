"""Shared conditional engineering utility accounting; never a measured yield."""
from __future__ import annotations

import math
from decimal import Decimal, localcontext, MAX_EMAX


def engineering_complexity_discount(
    K: int, B: int, lambda_K: float = 0.15, gamma_B: float = 0.965
) -> float:
    """Return the declared locus-complexity utility discount, without rounding.

    K counts distinct DSB/indel loci and B distinct base/prime-edit loci. A value
    of one means no *assigned utility penalty*, not perfect manufacturing yield.
    """
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in (K, B)):
        raise ValueError("K and B must be nonnegative integers")
    def finite_real(value):
        try:
            return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)
        except OverflowError:
            return False

    if not finite_real(lambda_K) or lambda_K < 0:
        raise ValueError("lambda_K must be finite and nonnegative")
    if not finite_real(gamma_B) or not 0 < gamma_B <= 1:
        raise ValueError("gamma_B must be finite in (0, 1]")
    excess = max(0, K - 4)
    if lambda_K == 0 or excess == 0:
        locus = 1.0
    else:
        try:
            exponent = lambda_K * excess
        except OverflowError:
            # Python first converts an arbitrary-size count to float. Its
            # product with a subnormal penalty can nevertheless be modest.
            with localcontext() as context:
                context.prec = 40
                context.Emax = MAX_EMAX
                exponent = float(Decimal.from_float(float(lambda_K)) * excess)
        locus = math.exp(-exponent)
    if gamma_B == 1 or B == 0:
        base = 1.0
    else:
        try:
            base = gamma_B ** B
        except OverflowError:
            # Any integer too large for binary64 with gamma < 1 has a
            # discount far below the smallest representable positive float.
            base = 0.0
    return locus * base
