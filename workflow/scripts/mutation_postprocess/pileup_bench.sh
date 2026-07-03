#!/bin/bash
# Instrumented copy of pileup.sh: reports per-PHASE wall-time and peak RSS.
# Same arguments and same commands as pileup.sh; only adds measurement.
# Run on a node with the FULL dataset (small test data will not reveal the
# real bottleneck). Outputs:
#   <OUTPUT_DIR>/pileup_bench.phases.log  : epoch <TAB> phase <TAB> START|END
#   <OUTPUT_DIR>/pileup_bench.rss.tsv     : epoch <TAB> total_RSS_KB (whole process group)
# and prints a per-phase summary table at the end.

set -uo pipefail   # NOTE: no errexit, so the summary always runs

INPUT_BAM=$1
PILEUP_TASKS=$2
REFERENCE_FA=$3
SAMPLE=$4
OUTPUT_DIR=$5
THREADS=$6
MEM_MB=${7:-8000}   # 割当メモリ(MB)。sort のバッファ上限に使う

mkdir -p "${OUTPUT_DIR}/pileup/workspace"
PHLOG="${OUTPUT_DIR}/pileup_bench.phases.log"
RSSLOG="${OUTPUT_DIR}/pileup_bench.rss.tsv"
: > "${PHLOG}"
: > "${RSSLOG}"

# --- background RSS sampler: sum RSS (KB) of every process in our process group
MYPGID=$(ps -o pgid= -p $$ | tr -d ' ')
sampler() {
    while :; do
        ts=$(date +%s)
        tot=$(ps -A -o pgid=,rss= | awk -v g="${MYPGID}" '$1==g{s+=$2} END{print s+0}')
        printf '%s\t%s\n' "${ts}" "${tot}" >> "${RSSLOG}"
        sleep 2
    done
}
sampler &
SAMPLER_PID=$!
trap 'kill ${SAMPLER_PID} 2>/dev/null' EXIT

phase() { printf '%s\t%s\t%s\n' "$(date +%s)" "$1" "$2" >> "${PHLOG}"; }

export INPUT_BAM REFERENCE_FA
run_pileup() {
    local bed_file=$1
    local pileup_file=$2
    samtools mpileup -l "${bed_file}" -f "${REFERENCE_FA}" --output-QNAME "${INPUT_BAM}" | \
        awk '{print $1 "\t" $2 -1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > "${pileup_file}"
}
export -f run_pileup

# ---------------- PHASE 1: parallel mpileup ----------------
phase mpileup_parallel START
xargs -n 2 -P "${THREADS}" bash -c 'run_pileup "$@"' _ < "${PILEUP_TASKS}"
phase mpileup_parallel END

# ---------------- PHASE 2: merge + sort --------------------
phase merge_sort START
cat ${OUTPUT_DIR}/pileup/workspace/*.pileup | sort -S $(( MEM_MB / 2 < 8192 ? MEM_MB / 2 : 8192 ))M --parallel="${THREADS}" -k 1,1 -k 2,2n > ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
phase merge_sort END

# ---------------- PHASE 3: bgzip ---------------------------
phase bgzip START
bgzip -@ ${THREADS} -f ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
phase bgzip END

# ---------------- PHASE 4: tabix ---------------------------
phase tabix START
tabix -p bed ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz
phase tabix END

kill ${SAMPLER_PID} 2>/dev/null

# ---------------- summary ----------------------------------
python3 - "${PHLOG}" "${RSSLOG}" <<'PY'
import sys
phlog, rsslog = sys.argv[1], sys.argv[2]
spans=[]; starts={}
for ln in open(phlog):
    ts,name,ev=ln.rstrip("\n").split("\t")
    ts=int(ts)
    if ev=="START": starts[name]=ts
    else: spans.append((name, starts[name], ts))
rss=[]
for ln in open(rsslog):
    a,b=ln.rstrip("\n").split("\t"); rss.append((int(a),int(b)))
print()
print("PHASE                 duration       peakRSS(GB)   meanRSS(GB)  samples")
print("-"*72)
total=0
for name,s,e in spans:
    total+=e-s
    win=[v for (t,v) in rss if s<=t<e]
    pk=max(win)/1024/1024 if win else 0
    mn=(sum(win)/len(win))/1024/1024 if win else 0
    dur=e-s
    print(f"{name:20s} {dur//3600}h{(dur%3600)//60:02d}m{dur%60:02d}s   {pk:10.1f}    {mn:10.1f}   {len(win)}")
gpk=max((v for _,v in rss),default=0)/1024/1024
print("-"*72)
print(f"{'TOTAL':20s} {total//3600}h{(total%3600)//60:02d}m{total%60:02d}s   {gpk:10.1f} (global peak)")
PY

echo 0
