use failure::Error;
use log::{debug, error, info};
use rayon::prelude::*;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use rust_htslib::bam::{self, Read};
use simplelog::{Config, LevelFilter, SimpleLogger};
use std::cmp;
use std::collections::HashSet;
use std::fs;
use std::io::prelude::*;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process;
use tempfile::tempdir;

pub struct SliceArgs<'a> {
    cell_barcodes: &'a HashSet<Vec<u8>>,
    i: usize,
    bam_file: &'a Path,
    tmp_dir: &'a Path,
    bam_tag: String,
    virtual_start: Option<i64>,
    virtual_stop: Option<i64>,
    to_replace: Option<String>,
    replacement: Option<String>,
}

pub fn load_barcodes(filename: impl AsRef<Path>) -> Result<HashSet<Vec<u8>>, Error> {
    let r = fs::File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut bc_set = HashSet::new();

    for l in reader.lines() {
        let seq = l?.into_bytes();
        bc_set.insert(seq);
    }
    let num_bcs = bc_set.len();
    if num_bcs == 0 {
        error!("Loaded 0 barcodes. Is your barcode file gzipped or empty?");
        process::exit(1);
    }
    debug!("Loaded {} barcodes", num_bcs);
    Ok(bc_set)
}

pub fn get_record_tag<'a>(rec: &'a Record, bam_tag: &str) -> Option<Vec<u8>> {
    let tag = rec.aux(bam_tag.as_bytes());
    match tag {
        Ok(t) => match t {
            Aux::String(t) => Some(t.as_bytes().to_vec()),
            _ => None,
        },
        Err(t) => None,
    }
}

pub fn load_writer(bam: &bam::Reader, out_bam_path: &Path) -> Result<bam::Writer, Error> {
    use rust_htslib::bam::Read; // collides with fs::Read
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr, bam::Format::Bam)?;
    Ok(out_handle)
}

pub fn bgzf_noffsets<P: AsRef<Path>>(
    bam_path: &P,
    num_chunks: &u64,
) -> Result<Vec<(Option<i64>, Option<i64>)>, Error> {
    use std::io::Read;

    fn vec_diff(input: &Vec<u64>) -> Vec<u64> {
        let vals = input.iter();
        let next_vals = input.iter().skip(1);

        vals.zip(next_vals).map(|(cur, next)| next - cur).collect()
    }

    // if we only have one thread, this is easy
    if *num_chunks == 1 as u64 {
        let final_offsets = vec![(None, None)];
        return Ok(final_offsets);
    }

    let bam_bytes = fs::metadata(bam_path)?.len();
    let mut initial_offsets = Vec::new();
    let step_size = bam_bytes / num_chunks;
    for n in 1..*num_chunks {
        initial_offsets.push((step_size * n) as u64);
    }

    let num_bytes = if initial_offsets.len() > 1 {
        let diff = vec_diff(&initial_offsets);
        let m = diff.iter().max().unwrap();
        cmp::min(1 << 16, *m)
    } else {
        1 << 16
    };

    // linear search to the right of each possible offset until
    // a valid virtual offset is found
    let mut adjusted_offsets = Vec::new();
    let mut fp = fs::File::open(bam_path)?;
    for offset in initial_offsets {
        fp.seek(io::SeekFrom::Start(offset))?;
        let mut buffer = [0; 2 << 16];
        fp.read(&mut buffer)?;
        for i in 0..num_bytes {
            if is_valid_bgzf_block(&buffer[i as usize..]) {
                adjusted_offsets.push(offset + i);
                break;
            }
        }
    }
    // bit-shift and produce start/stop intervals
    let mut final_offsets = Vec::new();

    // handle special case where we only found one offset
    if adjusted_offsets.len() == 1 {
        final_offsets.push((None, None));
        return Ok(final_offsets);
    }

    final_offsets.push((None, Some(((adjusted_offsets[1]) as i64) << 16)));
    for n in 2..num_chunks - 1 {
        let n = n as usize;
        final_offsets.push((
            Some((adjusted_offsets[n - 1] as i64) << 16),
            Some((adjusted_offsets[n] as i64) << 16),
        ));
    }
    final_offsets.push((
        Some(((adjusted_offsets[adjusted_offsets.len() - 1]) as i64) << 16),
        None,
    ));
    Ok(final_offsets)
}

pub fn is_valid_bgzf_block(block: &[u8]) -> bool {
    // look for the bgzip magic characters \x1f\x8b\x08\x04
    // TODO: is this sufficient?
    if block.len() < 18 {
        return false;
    }
    if (block[0] != 31) | (block[1] != 139) | (block[2] != 8) | (block[3] != 4) {
        return false;
    }
    true
}

pub fn read_bam_slice(args: &SliceArgs) -> Result<PathBuf, rust_htslib::tpool::Error> {
    let mut bam = bam::Reader::from_path(args.bam_file).unwrap();
    let out_bam_file = args.tmp_dir.join(format!("{}.bam", args.i));

    let mut out_bam = load_writer(&bam, &out_bam_file).unwrap();

    for r in bam.iter_chunk(args.virtual_start, args.virtual_stop) {
        let mut rec = r?;
        let tag = get_record_tag(&rec, &args.bam_tag);

        if args.to_replace.is_some() && tag.is_some() {
            substitute_text_in_tag(
                &mut rec,
                &args.bam_tag,
                &args.to_replace.clone().unwrap(),
                &args.replacement.clone().unwrap(),
            )
            .expect("Missing tag");
        }

        if tag.is_some() && args.cell_barcodes.contains(&tag.unwrap()) {
            out_bam.write(&rec).expect("Cannot write to temp bam file")
        }
    }

    Ok(out_bam_file.to_path_buf())
}

pub fn merge_bams<P: AsRef<Path>>(tmp_bams: Vec<&PathBuf>, out_bam_file: P) {
    use rust_htslib::bam::Read; // collides with fs::Read
    let bam = bam::Reader::from_path(tmp_bams[0].clone()).unwrap();
    let mut out_bam = load_writer(&bam, out_bam_file.as_ref()).unwrap();
    for b in tmp_bams.iter() {
        let mut rdr = bam::Reader::from_path(b).unwrap();
        for _rec in rdr.records() {
            let rec = _rec.unwrap();
            out_bam.write(&rec).unwrap();
        }
    }
}

fn substitute_text_in_tag(
    rec: &mut bam::Record,
    bam_tag: &str,
    to_replace: &str,
    replacement: &str,
) -> Result<(), rust_htslib::tpool::Error> {
    let bam_tag_bytes = bam_tag.as_bytes();
    let bc = get_record_tag(&rec, &bam_tag);

    match bc {
        Some(b) => match rec.remove_aux(&bam_tag_bytes) {
            Ok(res) => {
                let new_tag = std::str::from_utf8(&b)
                    .expect("Not UTF-8 formatted")
                    .replace(to_replace, replacement);
                rec.push_aux(&bam_tag_bytes, Aux::String(&new_tag)).unwrap();
                Ok(())
            }

            Err(res) => Err(rust_htslib::tpool::Error::BamAuxTagNotFound),
        },

        None => Err(rust_htslib::tpool::Error::BamAuxTagNotFound),
    }
}

pub fn subsample_bam<P: AsRef<Path>>(
    bam_file: P,
    barcodes_file: P,
    bam_tag: String,
    to_replace: Option<String>,
    replacement: Option<String>,
    out_bam_file: P,
    cores: usize,
) -> Result<PathBuf, Error> {

    let _ = SimpleLogger::init(LevelFilter::Info, Config::default());

    let cell_barcodes = load_barcodes(&barcodes_file).unwrap();
    let tmp_dir = tempdir().unwrap();
    let virtual_offsets = bgzf_noffsets(&bam_file, &(cores as u64)).unwrap();

    let mut chunks = Vec::new();

    for (i, (virtual_start, virtual_stop)) in virtual_offsets.iter().enumerate() {
        let c = SliceArgs {
            cell_barcodes: &cell_barcodes,
            i: i,
            bam_file: &bam_file.as_ref(),
            tmp_dir: tmp_dir.path(),
            bam_tag: bam_tag.clone(),
            virtual_start: *virtual_start,
            virtual_stop: *virtual_stop,
            to_replace: to_replace.clone(),
            replacement: replacement.clone(),
        };
        chunks.push(c);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cores as usize)
        .build()
        .unwrap();

    let results: Vec<_> = pool.install(|| {
        chunks
            .par_iter()
            .map(|chunk| read_bam_slice(chunk))
            .collect()
    });

    let tmp_bams: Vec<_> = results.iter().map(|r| r.as_ref().unwrap()).collect();
    merge_bams(tmp_bams, &out_bam_file);

    Ok(PathBuf::from(&out_bam_file.as_ref()))
}
