use pyo3::prelude::*;
pub mod subsample_bam;

/// Formats the sum of two numbers as string.
#[pyfunction]
#[pyo3(name = "subsample_bam")]
fn subsample_bam_py(
    bam_file: String,
    barcodes_file: String,
    bam_tag: String,
    to_replace: Option<String>,
    replacement: Option<String>,
    outfile: String,
    n_threads: usize,
)  -> PyResult<String> 
{   

    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    let out = subsample_bam::subsample_bam(
        bam_file,
        barcodes_file,
        bam_tag,
        to_replace,
        replacement,
        outfile,
        n_threads,
    );

    Ok(out.unwrap().as_path().display().to_string())

}

#[pymodule]
fn rust_bam_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(subsample_bam_py, m)?)?;

    Ok(())
}
