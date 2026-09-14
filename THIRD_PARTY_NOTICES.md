# Third-party components

The unmodified `code/native/vendor/bliss-0.77.zip` contains Bliss 0.77,
including its copyright notices and `COPYING` and `COPYING.LESSER` files.
The graph library sources identify their licence as GNU LGPL version 3.
The build verifies the archive's SHA-256 before extracting it. The complete
source is included; no prebuilt Bliss library is distributed.

The native protocol build also links to the separately installed GNU
Multiple Precision Arithmetic Library (GMP). Python interfaces use the
separately installed NumPy and `zstandard` packages. These dependencies
are not bundled here. Their respective notices accompany their own
distributions.
