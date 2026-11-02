#!/usr/bin/env bash
# Build NIST NBIS 5.0.0 (mindtct, bozorth3, nfiq, cjpegl) from the official source, on Linux.
# Revision step 1.3 (C3, reference matcher). Reproducible recipe; records the source hash.
#
#   bash Code/reference_matcher/build_nbis.sh <build_dir> <install_dir>
#
# Two changes to NIST's generated rules.mak are needed on a modern toolchain (gcc >= 10) and
# on non-x86 CPUs: the x86-only "-m64" flag is dropped and "-fcommon" is added (NBIS relies on
# tentative definitions that gcc >= 10 rejects by default). No source file is modified.
set -euo pipefail
BUILD=${1:?build dir}; INSTALL=${2:?install dir}
URL="https://nigos.nist.gov/nist/nbis/nbis_v5_0_0.zip"
SHA="0adf8ab0f6b0e4208de50ca00ba21d3d77112ecd66288757ddfed21f6bee92c3"
mkdir -p "$BUILD" "$INSTALL"; cd "$BUILD"
[ -f nbis_v5_0_0.zip ] || curl -sSL --max-time 300 -o nbis_v5_0_0.zip "$URL"
echo "$SHA  nbis_v5_0_0.zip" | sha256sum -c -
[ -d Rel_5.0.0 ] || unzip -q nbis_v5_0_0.zip
cd Rel_5.0.0
./setup.sh "$INSTALL" --without-X11 --STDLIBS --64 > setup.log 2>&1
sed -i 's/^ARCH_FLAG\(\s*\):= -m64 -fPIC/ARCH_FLAG\1:= -fPIC -fcommon/' rules.mak
grep -q '^ARCH_FLAG.*-fPIC -fcommon' rules.mak
make config > make_config.log 2>&1
make it > make_it.log 2>&1
make install LIBNBIS=no > make_install.log 2>&1
for t in mindtct bozorth3 nfiq cjpegl; do [ -x "$INSTALL/bin/$t" ] || { echo "missing $t"; exit 1; }; done
echo "NBIS 5.0.0 built: $INSTALL/bin (source sha256 $SHA; rules.mak ARCH_FLAG := -fPIC -fcommon)"
