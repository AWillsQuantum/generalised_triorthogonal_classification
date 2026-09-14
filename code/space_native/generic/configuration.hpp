#pragma once

#ifndef UTSP_LENGTH
#error Define UTSP_LENGTH as the support length.
#endif
#ifndef UTSP_DIMENSION
#error Define UTSP_DIMENSION as the affine dimension.
#endif

static_assert(16 <= UTSP_LENGTH && UTSP_LENGTH <= 54 && UTSP_LENGTH % 2 == 0);
static_assert(8 <= UTSP_DIMENSION && UTSP_DIMENSION <= 15);
