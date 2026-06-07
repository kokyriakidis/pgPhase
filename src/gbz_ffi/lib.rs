//! C FFI for gbz-base interval queries.
//!
//! Provides opaque handles to GBZBase/GAFBase databases and a function to
//! query overlapping reads for a genomic interval, returning GAF lines via
//! callback instead of spawning a subprocess and writing temp files.

use gbz_base::{GBZBase, GraphInterface, GAFBase};
use gbz_base::{Subgraph, SubgraphQuery, AlignmentOutput, ReadSet};
use gbz_base::GraphReference;
use gbz_base::utils;

use gbz::FullPathName;

use std::ffi::{CStr, c_char, c_int, c_void};
use std::ptr;

// ---------------------------------------------------------------------------
// Opaque handle: holds the opened databases for the lifetime of the session.
// ---------------------------------------------------------------------------

struct GbzSession {
    db: GBZBase,
}

struct GafSession {
    db: GAFBase,
}

// ---------------------------------------------------------------------------
// Lifecycle: open / close
// ---------------------------------------------------------------------------

/// Open a GBZ-base database. Returns an opaque handle, or null on error.
/// On error, if `err_out` is non-null, it receives a malloc'd error string
/// that the caller must free with `pgphase_gbz_free_string`.
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gbz_open(
    path: *const c_char,
    err_out: *mut *mut c_char,
) -> *mut c_void {
    let path = match unsafe { CStr::from_ptr(path) }.to_str() {
        Ok(s) => s,
        Err(e) => return set_error(err_out, &format!("invalid UTF-8 in path: {}", e)),
    };
    match GBZBase::open(path) {
        Ok(db) => {
            let session = Box::new(GbzSession { db });
            Box::into_raw(session) as *mut c_void
        }
        Err(e) => set_error(err_out, &e),
    }
}

/// Close a GBZ-base handle opened with `pgphase_gbz_open`.
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gbz_close(handle: *mut c_void) {
    if !handle.is_null() {
        unsafe { drop(Box::from_raw(handle as *mut GbzSession)); }
    }
}

/// Open a GAF-base database. Returns an opaque handle, or null on error.
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gaf_open(
    path: *const c_char,
    err_out: *mut *mut c_char,
) -> *mut c_void {
    let path = match unsafe { CStr::from_ptr(path) }.to_str() {
        Ok(s) => s,
        Err(e) => return set_error(err_out, &format!("invalid UTF-8 in path: {}", e)),
    };
    match GAFBase::open(path) {
        Ok(db) => {
            let session = Box::new(GafSession { db });
            Box::into_raw(session) as *mut c_void
        }
        Err(e) => set_error(err_out, &e),
    }
}

/// Close a GAF-base handle opened with `pgphase_gaf_open`.
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gaf_close(handle: *mut c_void) {
    if !handle.is_null() {
        unsafe { drop(Box::from_raw(handle as *mut GafSession)); }
    }
}

/// Free a string returned by the FFI (error messages).
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gbz_free_string(s: *mut c_char) {
    if !s.is_null() {
        unsafe { drop(std::ffi::CString::from_raw(s)); }
    }
}

/// Validate that a GBZ-base and GAF-base are compatible (built from the
/// same graph). Returns 0 if compatible, -1 on error.
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gbz_gaf_validate(
    gbz_handle: *mut c_void,
    gaf_handle: *mut c_void,
    err_out: *mut *mut c_char,
) -> c_int {
    if gbz_handle.is_null() || gaf_handle.is_null() {
        return set_error_int(err_out, "null GBZ or GAF handle");
    }
    let gbz_session = unsafe { &*(gbz_handle as *const GbzSession) };
    let gaf_session = unsafe { &*(gaf_handle as *const GafSession) };

    let mut graph_iface = match GraphInterface::new(&gbz_session.db) {
        Ok(g) => g,
        Err(e) => return set_error_int(err_out, &format!("GraphInterface::new: {}", e)),
    };
    let reference = match graph_iface.graph_name() {
        Ok(n) => n,
        Err(e) => return set_error_int(err_out, &format!("GBZ graph_name: {}", e)),
    };
    let alignments = match gaf_session.db.graph_name() {
        Ok(n) => n,
        Err(e) => return set_error_int(err_out, &format!("GAF graph_name: {}", e)),
    };
    if let Err(e) = utils::require_valid_reference(&alignments, &reference) {
        return set_error_int(err_out, &e);
    }
    0
}

// ---------------------------------------------------------------------------
// Interval query: extract overlapping reads as GAF lines
// ---------------------------------------------------------------------------

/// Callback type: called once per GAF line.
/// `line` is a pointer to the GAF line bytes (not null-terminated).
/// `line_len` is the length in bytes (excludes trailing newline).
/// `user_data` is the opaque pointer passed to the query function.
pub type GafLineCallback = extern "C" fn(
    line: *const u8,
    line_len: usize,
    user_data: *mut c_void,
);

/// Query overlapping reads for a genomic interval and invoke `callback`
/// once per GAF alignment line.
///
/// Returns 0 on success, -1 on error. On error, if `err_out` is non-null,
/// it receives a malloc'd error string.
///
/// Thread safety: each call must use its own `gbz_handle` or be externally
/// synchronized. The GBZBase uses SQLite which is not safe for concurrent
/// reads from the same connection.
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gbz_query_interval(
    gbz_handle: *mut c_void,
    gaf_handle: *mut c_void,
    sample: *const c_char,   // may be null
    contig: *const c_char,
    beg: u64,
    end: u64,
    callback: GafLineCallback,
    user_data: *mut c_void,
    err_out: *mut *mut c_char,
) -> c_int {
    // Validate handles.
    if gbz_handle.is_null() || gaf_handle.is_null() {
        return set_error_int(err_out, "null GBZ or GAF handle");
    }
    let gbz_session = unsafe { &*(gbz_handle as *const GbzSession) };
    let gaf_session = unsafe { &*(gaf_handle as *const GafSession) };

    // Parse string arguments.
    let contig_str = match unsafe { CStr::from_ptr(contig) }.to_str() {
        Ok(s) => s,
        Err(e) => return set_error_int(err_out, &format!("invalid contig: {}", e)),
    };
    let sample_str = if sample.is_null() {
        None
    } else {
        match unsafe { CStr::from_ptr(sample) }.to_str() {
            Ok(s) => Some(s),
            Err(e) => return set_error_int(err_out, &format!("invalid sample: {}", e)),
        }
    };

    // Build the path name for the interval query.
    let sample_name = sample_str.unwrap_or(gbz::GENERIC_SAMPLE);
    let path_name = FullPathName::reference(sample_name, contig_str);
    let query = SubgraphQuery::path_interval(&path_name, (beg as usize)..(end as usize));

    // Extract subgraph from the database.
    let mut graph_iface = match GraphInterface::new(&gbz_session.db) {
        Ok(g) => g,
        Err(e) => return set_error_int(err_out, &format!("GraphInterface::new: {}", e)),
    };
    let mut subgraph = Subgraph::new();
    if let Err(e) = subgraph.from_db(&mut graph_iface, &query) {
        return set_error_int(err_out, &format!("subgraph extraction: {}", e));
    }

    // Extract overlapping reads.
    let graph_ref = GraphReference::Db(&mut graph_iface);
    let read_set = match ReadSet::new(graph_ref, &subgraph, &gaf_session.db, AlignmentOutput::Overlapping) {
        Ok(rs) => rs,
        Err(e) => return set_error_int(err_out, &format!("ReadSet extraction: {}", e)),
    };

    // Serialize to GAF in memory, then split by newlines and invoke callback.
    let mut buf: Vec<u8> = Vec::new();
    if let Err(e) = read_set.to_gaf(&mut buf) {
        return set_error_int(err_out, &format!("to_gaf: {}", e));
    }

    // Each GAF line ends with '\n'. Split and invoke callback per line.
    for line in buf.split(|&b| b == b'\n') {
        if !line.is_empty() {
            callback(line.as_ptr(), line.len(), user_data);
        }
    }

    0
}

// ---------------------------------------------------------------------------
// Structured interval query: pass alignment fields directly, no GAF text
// ---------------------------------------------------------------------------

/// Callback type for structured alignment data.
/// Called once per alignment with the read name, GBWT node handles, and mapq.
///
/// GBWT node encoding: `handle = 2 * node_id + orientation` where
/// Forward=0, Reverse=1. Decode with `node_id = handle / 2`,
/// `is_reverse = (handle & 1) != 0`.
///
/// `name` is NOT null-terminated; use `name_len`.
pub type AlignmentCallback = extern "C" fn(
    name: *const u8,
    name_len: usize,
    nodes: *const u64,
    node_count: usize,
    mapq: i32,
    user_data: *mut c_void,
);

/// Query overlapping reads for a genomic interval and invoke `callback`
/// once per alignment with structured data (no GAF serialization).
///
/// Returns 0 on success, -1 on error.
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gbz_query_interval_structured(
    gbz_handle: *mut c_void,
    gaf_handle: *mut c_void,
    sample: *const c_char,
    contig: *const c_char,
    beg: u64,
    end: u64,
    callback: AlignmentCallback,
    user_data: *mut c_void,
    err_out: *mut *mut c_char,
) -> c_int {
    if gbz_handle.is_null() || gaf_handle.is_null() {
        return set_error_int(err_out, "null GBZ or GAF handle");
    }
    let gbz_session = unsafe { &*(gbz_handle as *const GbzSession) };
    let gaf_session = unsafe { &*(gaf_handle as *const GafSession) };

    let contig_str = match unsafe { CStr::from_ptr(contig) }.to_str() {
        Ok(s) => s,
        Err(e) => return set_error_int(err_out, &format!("invalid contig: {}", e)),
    };
    let sample_str = if sample.is_null() {
        None
    } else {
        match unsafe { CStr::from_ptr(sample) }.to_str() {
            Ok(s) => Some(s),
            Err(e) => return set_error_int(err_out, &format!("invalid sample: {}", e)),
        }
    };

    let sample_name = sample_str.unwrap_or(gbz::GENERIC_SAMPLE);
    let path_name = FullPathName::reference(sample_name, contig_str);
    let query = SubgraphQuery::path_interval(&path_name, (beg as usize)..(end as usize));

    let mut graph_iface = match GraphInterface::new(&gbz_session.db) {
        Ok(g) => g,
        Err(e) => return set_error_int(err_out, &format!("GraphInterface::new: {}", e)),
    };
    let mut subgraph = Subgraph::new();
    if let Err(e) = subgraph.from_db(&mut graph_iface, &query) {
        return set_error_int(err_out, &format!("subgraph extraction: {}", e));
    }

    let graph_ref = GraphReference::Db(&mut graph_iface);
    let read_set = match ReadSet::new(graph_ref, &subgraph, &gaf_session.db, AlignmentOutput::Overlapping) {
        Ok(rs) => rs,
        Err(e) => return set_error_int(err_out, &format!("ReadSet extraction: {}", e)),
    };

    // Emit each alignment as structured data.
    // On 64-bit platforms usize == u64, so we can pass the GBWT handle
    // slice directly without allocating a conversion Vec per alignment.
    #[cfg(not(target_pointer_width = "64"))]
    compile_error!("pgphase FFI requires a 64-bit target (usize == u64)");

    for alignment in read_set.iter() {
        let path = match alignment.target_path() {
            Some(p) => p,
            None => continue, // skip alignments without an extracted path
        };
        let mapq = alignment.mapq.unwrap_or(255) as i32;

        // Safety: on 64-bit targets usize and u64 have identical size,
        // alignment, and representation. The slice is valid for the
        // duration of the callback.
        let nodes_ptr = path.as_ptr() as *const u64;

        callback(
            alignment.name.as_ptr(),
            alignment.name.len(),
            nodes_ptr,
            path.len(),
            mapq,
            user_data,
        );
    }

    0
}

// ---------------------------------------------------------------------------
// Path range query: expose the genome-coordinate extent of a reference path
// ---------------------------------------------------------------------------

/// Query the genome-coordinate range covered by a reference path.
///
/// For subgraphed GBZ files the path may start at a non-zero genome offset
/// (the "fragment" field).  This function returns that offset as `*out_start`
/// and a conservative end coordinate as `*out_end`.  The end is computed from
/// the last indexed position in the ReferenceIndex, so the true path may
/// extend up to ~INDEX_INTERVAL (1000 bp) beyond `*out_end`.
///
/// Returns 0 on success, -1 on error (path not found / not indexed).
#[unsafe(no_mangle)]
pub extern "C" fn pgphase_gbz_path_range(
    gbz_handle: *mut c_void,
    sample: *const c_char,
    contig: *const c_char,
    out_start: *mut u64,
    out_end: *mut u64,
    err_out: *mut *mut c_char,
) -> c_int {
    if gbz_handle.is_null() {
        return set_error_int(err_out, "null GBZ handle");
    }
    let gbz_session = unsafe { &*(gbz_handle as *const GbzSession) };

    let contig_str = match unsafe { CStr::from_ptr(contig) }.to_str() {
        Ok(s) => s,
        Err(e) => return set_error_int(err_out, &format!("invalid contig: {}", e)),
    };
    let sample_str = if sample.is_null() {
        None
    } else {
        match unsafe { CStr::from_ptr(sample) }.to_str() {
            Ok(s) => Some(s),
            Err(e) => return set_error_int(err_out, &format!("invalid sample: {}", e)),
        }
    };

    let sample_name = sample_str.unwrap_or(gbz::GENERIC_SAMPLE);
    // Use a large fragment so find_path's "fragment <= ?" query matches
    // any path regardless of its starting offset.  Stay within i64::MAX
    // to avoid overflow when rusqlite binds usize as SQLite INTEGER.
    let mut path_name = FullPathName::reference(sample_name, contig_str);
    path_name.fragment = i64::MAX as usize;

    let mut graph_iface = match GraphInterface::new(&gbz_session.db) {
        Ok(g) => g,
        Err(e) => return set_error_int(err_out, &format!("GraphInterface::new: {}", e)),
    };

    // Look up the path record to get the fragment (genome-coordinate offset).
    let path = match graph_iface.find_path(&path_name) {
        Ok(Some(p)) => p,
        Ok(None) => return set_error_int(err_out, &format!(
            "no path found for {}#{}", sample_name, contig_str)),
        Err(e) => return set_error_int(err_out, &format!("find_path: {}", e)),
    };
    if !path.is_indexed {
        return set_error_int(err_out, &format!(
            "path {}#{} is not indexed for random access", sample_name, contig_str));
    }

    let fragment_start = path.name.fragment as u64;

    // Find the last indexed position by querying with a very large offset.
    // indexed_position returns the last entry with path_offset <= query.
    // Stay within i64::MAX to avoid rusqlite integer overflow.
    let max_indexed_offset = match graph_iface.indexed_position(path.handle, i64::MAX as usize) {
        Ok(Some((offset, _pos))) => offset as u64,
        Ok(None) => 0,
        Err(e) => return set_error_int(err_out, &format!("indexed_position: {}", e)),
    };

    unsafe {
        *out_start = fragment_start;
        *out_end = fragment_start + max_indexed_offset;
    }

    0
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn set_error(err_out: *mut *mut c_char, msg: &str) -> *mut c_void {
    if !err_out.is_null() {
        if let Ok(c) = std::ffi::CString::new(msg) {
            unsafe { *err_out = c.into_raw(); }
        }
    }
    ptr::null_mut()
}

fn set_error_int(err_out: *mut *mut c_char, msg: &str) -> c_int {
    set_error(err_out, msg);
    -1
}
