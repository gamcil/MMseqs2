#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -e "$1" ]
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

notExists "$1.dbtype" && echo "$1.dbtype not found!" && exit 1;
notExists "$2" && echo "$2 not found!" && exit 1;

hasCommand sort || fail "sort not found in PATH"

TARGET="$1"
INPUT="$2"
OUTDB="$3"
TMP_PATH="$4"
LOOKUP="${TARGET}.lookup"
notExists "$LOOKUP" && echo "$LOOKUP not found!" && exit 1;

TAB=$(printf '\t')
SORT_MEMORY_PAR=""
if [ -n "$SORT_MEMORY" ]; then
    SORT_MEMORY_PAR="-S ${SORT_MEMORY}"
fi
SORT_PARALLEL_PAR=""
if sort --help 2>&1 | grep -q -- "--parallel"; then
    SORT_PARALLEL_PAR="--parallel=${SORT_THREADS}"
fi

FEATURES_BY_ID="${TMP_PATH}/features.by_id"
LOOKUP_SORTED="${TMP_PATH}/lookup.sorted"
FEATURES_RESOLVED="${TMP_PATH}/features.resolved"
FEATURES_RESOLVED_SORTED="${TMP_PATH}/features.resolved.sorted"
CONTEXT_ROWS="${TMP_PATH}/context.rows"
CONTEXT_ROWS_SORTED="${TMP_PATH}/context.rows.sorted"

# Sort by first field accession for matching to sequence db
# shellcheck disable=SC2086
LC_ALL=C sort -t "$TAB" ${SORT_MEMORY_PAR} ${SORT_PARALLEL_PAR} -T "$TMP_PATH" \
    -k1,1 \
    "$INPUT" > "$FEATURES_BY_ID" \
    || fail "feature table id sort failed"

if [ "$CONTEXT_ID_MODE" = "0" ]; then
    LOOKUP_SORT_KEYS="-k1,1"
else
    LOOKUP_SORT_KEYS="-k2,2 -k1,1"
fi

# Sort lookup by the same identifier mode as before
# shellcheck disable=SC2086
LC_ALL=C sort -t "$TAB" ${SORT_MEMORY_PAR} ${SORT_PARALLEL_PAR} -T "$TMP_PATH" \
    ${LOOKUP_SORT_KEYS} \
    "$LOOKUP" > "$LOOKUP_SORTED" \
    || fail "lookup sort failed"

# Resolve target identifiers to mmseqs db keys
# shellcheck disable=SC2086
"$MMSEQS" createcontextresolve "$LOOKUP_SORTED" "$FEATURES_BY_ID" "$FEATURES_RESOLVED" ${CREATECONTEXTRESOLVE_PAR} \
    || fail "createcontextresolve died"

# Sort features into genomic order by scaffold and coordinates
# shellcheck disable=SC2086
LC_ALL=C sort -t "$TAB" ${SORT_MEMORY_PAR} ${SORT_PARALLEL_PAR} -T "$TMP_PATH" \
    -k4,4 -k5,5n -k6,6n -k3,3 \
    "$FEATURES_RESOLVED" > "$FEATURES_RESOLVED_SORTED" \
    || fail "resolved feature coordinate sort failed"

# Collapse redundant features and write complete context rows
# shellcheck disable=SC2086
"$MMSEQS" createcontextcontexts "$FEATURES_RESOLVED_SORTED" "$CONTEXT_ROWS" ${CREATECONTEXTCONTEXTS_PAR} \
    || fail "createcontextcontexts died"

# Group context rows by target key
# shellcheck disable=SC2086
LC_ALL=C sort -t "$TAB" ${SORT_MEMORY_PAR} ${SORT_PARALLEL_PAR} -T "$TMP_PATH" \
    -k1,1n \
    "$CONTEXT_ROWS" > "$CONTEXT_ROWS_SORTED" \
    || fail "context sort failed"

# Pack grouped context rows into a standalone mmseqs context db
# shellcheck disable=SC2086
"$MMSEQS" createcontextdbcore "$CONTEXT_ROWS_SORTED" "$OUTDB" ${CREATECONTEXTDBCORE_PAR} \
    || fail "createcontextdbcore died"

if [ -n "$REMOVE_TMP" ]; then
    rm -f "$FEATURES_BY_ID" "$LOOKUP_SORTED" "$FEATURES_RESOLVED" "$FEATURES_RESOLVED_SORTED" "$CONTEXT_ROWS" "$CONTEXT_ROWS_SORTED" "${TMP_PATH}/createcontextdb.sh"
fi
