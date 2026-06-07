/**
 * @file gbz_ffi.h
 * @brief C interface to gbz-base for in-process GBZ/GAF-base interval queries.
 *
 * Replaces subprocess-based query invocations with direct library calls.
 * The GBZ and GAF databases are opened once and reused across all chunk queries.
 */

#ifndef PGPHASE_GBZ_FFI_H
#define PGPHASE_GBZ_FFI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Lifecycle ─────────────────────────────────────────────────────────── */

/**
 * Open a GBZ-base database. Returns an opaque handle, or NULL on error.
 * On error, *err_out (if non-null) receives a malloc'd string; free with
 * pgphase_gbz_free_string().
 */
void* pgphase_gbz_open(const char* path, char** err_out);

/** Close a GBZ-base handle. */
void pgphase_gbz_close(void* handle);

/**
 * Open a GAF-base database. Returns an opaque handle, or NULL on error.
 */
void* pgphase_gaf_open(const char* path, char** err_out);

/** Close a GAF-base handle. */
void pgphase_gaf_close(void* handle);

/** Free a string returned by the FFI (error messages). */
void pgphase_gbz_free_string(char* s);

/**
 * Validate that a GBZ-base and GAF-base are compatible (built from the
 * same graph). Returns 0 if compatible, -1 on error.
 */
int pgphase_gbz_gaf_validate(void* gbz_handle, void* gaf_handle, char** err_out);

/* ── Interval query ────────────────────────────────────────────────────── */

/**
 * Callback invoked once per GAF alignment line.
 * @param line     Pointer to GAF line bytes (NOT null-terminated).
 * @param line_len Length in bytes (excludes trailing newline).
 * @param user_data Opaque pointer passed through from the query call.
 */
typedef void (*pgphase_gaf_line_callback)(
    const unsigned char* line,
    size_t line_len,
    void* user_data
);

/**
 * Query overlapping reads for a genomic interval and invoke @p callback
 * once per GAF alignment line.
 *
 * @param gbz_handle  Handle from pgphase_gbz_open().
 * @param gaf_handle  Handle from pgphase_gaf_open().
 * @param sample      Reference sample name (e.g. "CHM13"), or NULL.
 * @param contig      Contig name (e.g. "chr20").
 * @param beg         0-based start coordinate.
 * @param end         0-based end coordinate (exclusive).
 * @param callback    Called once per GAF line.
 * @param user_data   Passed through to callback.
 * @param err_out     On error, receives a malloc'd error string.
 * @return 0 on success, -1 on error.
 *
 * Thread safety: each concurrent call needs its own gbz_handle (one
 * SQLite connection per thread). gaf_handle may also need to be
 * per-thread depending on the SQLite threading mode.
 */
int pgphase_gbz_query_interval(
    void* gbz_handle,
    void* gaf_handle,
    const char* sample,
    const char* contig,
    uint64_t beg,
    uint64_t end,
    pgphase_gaf_line_callback callback,
    void* user_data,
    char** err_out
);

/* ── Structured interval query (no GAF text serialization) ─────────────── */

/**
 * Callback invoked once per alignment with structured fields.
 *
 * GBWT node encoding: handle = 2 * node_id + orientation
 *   (Forward=0, Reverse=1).
 * Decode: node_id = handle / 2, is_reverse = (handle & 1) != 0.
 *
 * @param name       Read name bytes (NOT null-terminated).
 * @param name_len   Length of name in bytes.
 * @param nodes      Array of GBWT node handles (target path).
 * @param node_count Number of elements in nodes.
 * @param mapq       Mapping quality (255 = missing).
 * @param user_data  Opaque pointer passed through from the query call.
 */
typedef void (*pgphase_alignment_callback)(
    const unsigned char* name,
    size_t name_len,
    const uint64_t* nodes,
    size_t node_count,
    int mapq,
    void* user_data
);

/**
 * Query overlapping reads for a genomic interval and invoke @p callback
 * once per alignment with structured data (no GAF serialization).
 *
 * Same parameters as pgphase_gbz_query_interval except the callback type.
 * Returns 0 on success, -1 on error.
 */
int pgphase_gbz_query_interval_structured(
    void* gbz_handle,
    void* gaf_handle,
    const char* sample,
    const char* contig,
    uint64_t beg,
    uint64_t end,
    pgphase_alignment_callback callback,
    void* user_data,
    char** err_out
);

/* ── Path range query ──────────────────────────────────────────────────── */

/**
 * Query the genome-coordinate range covered by a reference path.
 *
 * For subgraphed GBZ files the path may start at a non-zero genome offset
 * (the "fragment" field).  This function returns that offset as *out_start
 * and a conservative end coordinate as *out_end (based on the last indexed
 * position; the true path may extend up to ~1000 bp beyond *out_end).
 *
 * @param gbz_handle  Handle from pgphase_gbz_open().
 * @param sample      Reference sample name (e.g. "CHM13"), or NULL.
 * @param contig      Contig name (e.g. "chr6").
 * @param out_start   Receives the genome-coordinate start of the path.
 * @param out_end     Receives a conservative genome-coordinate end.
 * @param err_out     On error, receives a malloc'd error string.
 * @return 0 on success, -1 on error (path not found / not indexed).
 */
int pgphase_gbz_path_range(
    void* gbz_handle,
    const char* sample,
    const char* contig,
    uint64_t* out_start,
    uint64_t* out_end,
    char** err_out
);

#ifdef __cplusplus
}
#endif

#endif /* PGPHASE_GBZ_FFI_H */
