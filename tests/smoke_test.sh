#!/usr/bin/env bash
# ============================================================
# smoke_test.sh  —  Structural dry-run / smoke test for
#                   TTA_HBV-OBI_NGS Nextflow pipeline
# ============================================================
# This script:
#   1. Validates the Python helper scripts with unit checks
#   2. Validates the test samplesheet
#   3. Runs Nextflow in -stub mode (no actual tools needed)
#
# Usage:
#   bash tests/smoke_test.sh [--stub-only] [--python-only]
#
# Requirements:
#   Python >=3.9, pandas, biopython (for Python checks)
#   Nextflow >=23.04 (for stub run)
# ============================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="${REPO_ROOT}/tests"
PASS=0
FAIL=0

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok()   { echo -e "${GREEN}  [PASS]${NC} $*"; PASS=$((PASS+1)); }
fail() { echo -e "${RED}  [FAIL]${NC} $*"; FAIL=$((FAIL+1)); }
info() { echo -e "${YELLOW}  [INFO]${NC} $*"; }

echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║    TTA_HBV-OBI_NGS  —  Smoke / Structural Test      ║"
echo "╚══════════════════════════════════════════════════════╝"
echo ""

# ── 1. Check required files exist ────────────────────────

info "Checking required pipeline files..."

REQUIRED_FILES=(
    "${REPO_ROOT}/main.nf"
    "${REPO_ROOT}/nextflow.config"
    "${REPO_ROOT}/environment.yml"
    "${REPO_ROOT}/workflow/modules/fastqc.nf"
    "${REPO_ROOT}/workflow/modules/fastp.nf"
    "${REPO_ROOT}/workflow/modules/multiqc.nf"
    "${REPO_ROOT}/workflow/modules/bwa_mem2.nf"
    "${REPO_ROOT}/workflow/modules/samtools.nf"
    "${REPO_ROOT}/workflow/modules/ivar.nf"
    "${REPO_ROOT}/workflow/modules/mafft.nf"
    "${REPO_ROOT}/workflow/modules/iqtree.nf"
    "${REPO_ROOT}/workflow/modules/annotate.nf"
    "${REPO_ROOT}/workflow/modules/genotype.nf"
    "${REPO_ROOT}/workflow/modules/report.nf"
    "${REPO_ROOT}/workflow/bin/validate_samplesheet.py"
    "${REPO_ROOT}/workflow/bin/annotate_variants.py"
    "${REPO_ROOT}/workflow/bin/genotype_assign.py"
    "${REPO_ROOT}/workflow/bin/generate_report.py"
    "${REPO_ROOT}/resources/knowledge_base/obi_mutations.tsv"
    "${REPO_ROOT}/schema/samplesheet_schema.json"
    "${REPO_ROOT}/data/samplesheet_template.csv"
    "${REPO_ROOT}/tests/samplesheet_test.csv"
    "${REPO_ROOT}/tests/data/hbv_ref_small.fasta"
)

for f in "${REQUIRED_FILES[@]}"; do
    if [[ -f "$f" ]]; then
        ok "File exists: $(basename "$f")"
    else
        fail "Missing file: $f"
    fi
done

# ── 2. Python syntax checks ───────────────────────────────

info "Checking Python script syntax..."

for py in "${REPO_ROOT}/workflow/bin/"*.py; do
    if python3 -m py_compile "$py" 2>&1; then
        ok "Syntax OK: $(basename "$py")"
    else
        fail "Syntax error in: $py"
    fi
done

# ── 3. Validate test samplesheet ────────────────────────

info "Validating test samplesheet..."

SAMPLESHEET="${TESTS_DIR}/samplesheet_test.csv"
if python3 "${REPO_ROOT}/workflow/bin/validate_samplesheet.py" \
        --samplesheet "${SAMPLESHEET}" 2>&1; then
    ok "Test samplesheet validation passed"
else
    fail "Test samplesheet validation failed"
fi

# ── 4. OBI knowledge base sanity check ──────────────────

info "Checking OBI knowledge base..."

OBI_DB="${REPO_ROOT}/resources/knowledge_base/obi_mutations.tsv"
N_ROWS=$(tail -n +2 "${OBI_DB}" | grep -c "^[^#]" || true)
if [[ "${N_ROWS}" -ge 10 ]]; then
    ok "OBI knowledge base has ${N_ROWS} entries"
else
    fail "OBI knowledge base has only ${N_ROWS} entries (expected >=10)"
fi

# ── 5. annotate_variants.py unit test ───────────────────

info "Running annotate_variants.py unit test..."

TMPDIR_TEST="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_TEST"' EXIT

# Create minimal iVar-style variants TSV
cat > "${TMPDIR_TEST}/test.ivar.tsv" <<'EOF'
REGION	POS	REF	ALT	REF_DP	REF_RV	REF_QUAL	ALT_DP	ALT_RV	ALT_QUAL	ALT_FREQ	TOTAL_DP	PVAL	PASS
HBV_genotype_A2	578	C	T	100	50	35	10	5	35	0.09	110	0.001	TRUE
HBV_genotype_A2	653	G	A	90	45	35	15	7	35	0.14	105	0.0001	TRUE
EOF

if python3 "${REPO_ROOT}/workflow/bin/annotate_variants.py" \
        --sample TEST_SAMPLE \
        --variants "${TMPDIR_TEST}/test.ivar.tsv" \
        --obi_db "${OBI_DB}" \
        --out_annotated "${TMPDIR_TEST}/annotated.tsv" \
        --out_flags "${TMPDIR_TEST}/flags.tsv" 2>&1; then
    ok "annotate_variants.py ran successfully"
    if [[ -f "${TMPDIR_TEST}/annotated.tsv" ]]; then
        N_LINES=$(wc -l < "${TMPDIR_TEST}/annotated.tsv")
        ok "Annotation output has ${N_LINES} line(s)"
    fi
else
    fail "annotate_variants.py failed"
fi

# ── 6. generate_report.py unit test ─────────────────────

info "Running generate_report.py unit test..."

# Create minimal inputs
cat > "${TMPDIR_TEST}/s1.genotype.tsv" <<'EOF'
sample	genotype	subgenotype	method	confidence	notes
TEST_SAMPLE	A	A2	pairwise_identity	high	stub
EOF

if python3 "${REPO_ROOT}/workflow/bin/generate_report.py" \
        --annotated "${TMPDIR_TEST}/annotated.tsv" \
        --genotypes "${TMPDIR_TEST}/s1.genotype.tsv" \
        --depths    "" \
        --multiqc   "" \
        --outdir    "${TMPDIR_TEST}" \
        --html 2>&1; then
    ok "generate_report.py ran successfully"
    if [[ -f "${TMPDIR_TEST}/cohort_summary.tsv" ]]; then
        ok "cohort_summary.tsv generated"
    fi
else
    fail "generate_report.py failed"
fi

# ── 7. Nextflow stub-run (if Nextflow available) ─────────

if command -v nextflow &>/dev/null; then
    info "Running Nextflow stub-run..."
    cd "${REPO_ROOT}"
    if nextflow run main.nf \
            --samplesheet "${TESTS_DIR}/samplesheet_test.csv" \
            --hbv_ref "${TESTS_DIR}/data/hbv_ref_small.fasta" \
            --genotype_refs "${TESTS_DIR}/data/hbv_ref_small.fasta" \
            --outdir "${TMPDIR_TEST}/nf_stub_out" \
            -stub-run \
            -profile test \
            2>&1 | tail -30; then
        ok "Nextflow stub-run completed successfully"
    else
        fail "Nextflow stub-run failed"
    fi
else
    info "Nextflow not found in PATH — skipping stub-run (install Nextflow to enable)"
fi

# ── Summary ───────────────────────────────────────────────

echo ""
echo "══════════════════════════════════════════════════════"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "══════════════════════════════════════════════════════"

if [[ "${FAIL}" -gt 0 ]]; then
    echo -e "${RED}  SMOKE TEST FAILED${NC}"
    exit 1
else
    echo -e "${GREEN}  SMOKE TEST PASSED${NC}"
    exit 0
fi
