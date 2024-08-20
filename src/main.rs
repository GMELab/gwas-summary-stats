use std::{
    collections::{HashMap, HashSet},
    io::Write,
    mem::MaybeUninit,
    path::Path,
    sync::Mutex,
};

use clap::Parser;
use rayon::prelude::*;
use tracing::{debug, error, info, warn};

const GOOGLE_SHEETS_API_KEY: &str = "AIzaSyA91UNqny43WENob6M3VpLKS0ayr-H-Lcw";
const COLS_MUST_BE_PRESENT: [&str; 20] = [
    "rsid",
    "chr",
    "pos",
    "ref",
    "alt",
    "effect_size",
    "effect_is_OR",
    "standard_error",
    "EAF",
    "pvalue",
    "pvalue_het",
    "N_total_column",
    "N_case_column",
    "N_ctrl_column",
    "column_delim",
    "hg_version",
    "file_path",
    "N_total",
    "N_case",
    "N_ctrl",
];
const COLS_MUST_NOT_BE_NA: [&str; 4] = ["chr", "pos", "ref", "alt"];
const ASSIGN_COL_NAMES: [&str; 13] = [
    "rsid",
    "chr",
    "pos",
    "ref",
    "alt",
    "effect_size",
    "standard_error",
    "EAF",
    "pvalue",
    "pvalue_het",
    "N_total_column",
    "N_case_column",
    "N_ctrl_column",
];

#[derive(Clone, Debug, clap::Parser)]
#[command(version)]
pub struct Args {
    #[arg(short, long)]
    google_sheets_id: String,
    #[arg(short, long)]
    trait_name:       String,
    #[arg(short = 'i', long)]
    raw_input_dir:    String,
    #[arg(short, long)]
    liftover:         String,
    #[arg(long)]
    liftover_dir:     String,
    #[arg(short = 'r', long)]
    grs_dir:          String,
    #[arg(short, long)]
    dbsnp_file:       String,
    #[arg(short, long)]
    samtools:         String,
    #[arg(short, long)]
    fasta_ref:        String,
    #[arg(short, long)]
    output_file:      String,
    #[arg(short, long)]
    samtools_threads: Option<usize>,
}

pub struct Ctx {
    args:  Args,
    sheet: Data,
}

#[derive(Clone)]
pub struct Data {
    // raw:    String,
    header: Vec<String>,
    data:   Vec<Vec<String>>,
}

impl Data {
    #[track_caller]
    pub fn idx(&self, key: &str) -> usize {
        self.idx_opt(key).unwrap()
    }

    pub fn idx_opt(&self, key: &str) -> Option<usize> {
        self.header.iter().position(|x| x == key)
    }

    pub fn col(&self, key: &str) -> impl Iterator<Item = &'_ str> {
        let idx = self.idx(key);
        self.data.iter().map(move |x| x[idx].as_str())
    }

    pub fn matching_rows(
        &self,
        key: &str,
        f: impl Fn(&str) -> bool,
    ) -> impl Iterator<Item = &'_ [String]> {
        let idx = self.idx(key);
        debug!(key, idx, "Matching rows");
        self.data
            .iter()
            .filter(move |x| f(x[idx].as_str()))
            .map(|x| x.as_slice())
    }

    pub fn get_from_row<'a>(&self, row: &'a [String], key: &str) -> &'a String {
        &row[self.idx(key)]
    }

    pub fn col_mut(&mut self, key: &str) -> impl Iterator<Item = &'_ mut String> {
        debug!(key, "Mutating column");
        let idx = self.idx(key);
        debug!(key, idx, "Mutating column");
        self.data.iter_mut().map(move |x| &mut x[idx])
    }

    pub fn write(&self, name: impl AsRef<Path>) {
        let file = std::fs::File::create(name).unwrap();
        let mut writer = flate2::write::GzEncoder::new(&file, flate2::Compression::default());
        debug!(len = self.data.len(), "Writing rows",);
        writeln!(writer, "{}", self.header.join("\t")).unwrap();
        for r in &self.data {
            writeln!(writer, "{}", r.join("\t")).unwrap();
        }
        writer.finish().unwrap();
    }

    #[track_caller]
    pub fn reorder(&mut self, new_order: &[&str]) {
        let new_order_idxs = new_order
            .iter()
            .map(|x| self.idx_opt(x))
            .collect::<Vec<_>>();
        let new_len = new_order.len();
        let data = std::mem::take(&mut self.data);
        self.data = data
            .into_par_iter()
            .map(|mut r| {
                let mut new_r = Vec::with_capacity(new_len);
                for idx in &new_order_idxs {
                    match idx {
                        Some(idx) => new_r.push(std::mem::take(&mut r[*idx])),
                        None => new_r.push("NA".to_string()),
                    }
                }
                new_r
            })
            .collect::<Vec<_>>();
        self.header = new_order.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    }

    pub fn read(delim: char, mut file: impl std::io::Read, has_header: bool) -> Self {
        let mut raw = String::new();
        file.read_to_string(&mut raw).unwrap();
        let (header, content) = if has_header {
            let (header, content) = raw.split_once('\n').unwrap();
            let header = header
                .split(delim)
                // .map(|x| unsafe { String::from_raw_parts(x.as_ptr().cast_mut(), x.len(), x.len()) })
                .map(|x| x.to_string())
                .collect::<Vec<_>>();
            (header, content)
        } else {
            (vec![], raw.as_str())
        };
        let data = content
            .par_lines()
            .map(|x| {
                x.split(delim)
                    // .map(|x| unsafe {
                    //     String::from_raw_parts(x.as_ptr().cast_mut(), x.len(), x.len())
                    // })
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        // Data { raw, header, data }
        Data { header, data }
    }
}

fn read_raw_data(delim: &str, file: impl std::io::Read) -> Data {
    let delim = if delim == "\t" || delim == "tab" {
        '\t'
    } else if delim == "," || delim == "comma" {
        ','
    } else if delim == "space" {
        ' '
    } else {
        error!("Invalid column delimiter {}", delim);
        panic!();
    };
    Data::read(delim, file, true)
}

fn reserve_to(r: &mut Vec<String>, len: usize) -> usize {
    match len.checked_sub(r.capacity()) {
        Some(res) => {
            r.reserve_exact(res);
            res
        },
        None => 0,
    }
}

#[tracing::instrument(skip(ctx))]
fn preformat(ctx: &Ctx) -> Data {
    let rows = ctx
        .sheet
        .matching_rows("trait_name", |x| x == ctx.args.trait_name)
        .collect::<Vec<_>>();
    if rows.is_empty() {
        error!(
            "No rows found in the GWAS formatting legend for trait_name={}",
            ctx.args.trait_name
        );
        panic!();
    }
    if rows.len() > 1 {
        error!(
            "Multiple rows found in the GWAS formatting legend for trait_name={}",
            ctx.args.trait_name
        );
        panic!();
    }
    let row = rows[0];
    for col in COLS_MUST_BE_PRESENT.iter() {
        let val = ctx.sheet.get_from_row(row, col);
        if val.is_empty() {
            error!(
                "Column {} is missing in the GWAS formatting legend for trait_name={}",
                col, ctx.args.trait_name
            );
            panic!();
        }
    }
    for col in COLS_MUST_NOT_BE_NA.iter() {
        let val = ctx.sheet.get_from_row(row, col);
        if val == "NA" || val == "NaN" {
            error!(
                "Column {} is NA in the GWAS formatting legend for trait_name={}",
                col, ctx.args.trait_name
            );
            panic!();
        }
    }
    let raw_input_dir = std::path::Path::new(&ctx.args.raw_input_dir);
    if !raw_input_dir.exists() {
        error!(
            "Raw input directory {} does not exist",
            ctx.args.raw_input_dir
        );
        panic!();
    }
    if !raw_input_dir.is_dir() {
        error!(
            "Raw input directory {} is not a directory",
            ctx.args.raw_input_dir
        );
        panic!();
    }
    let mut file_path = ctx.sheet.get_from_row(row, "file_path").as_str();
    if file_path.starts_with('/') {
        file_path = file_path.strip_prefix('/').unwrap();
    }
    let raw_input_file = raw_input_dir.join(file_path);
    if !raw_input_file.exists() {
        error!(
            "Raw input file {} does not exist",
            raw_input_file.to_string_lossy()
        );
        panic!();
    }
    if !raw_input_file.is_file() {
        error!(
            "Raw input file {} is not a file",
            raw_input_file.to_string_lossy()
        );
        panic!();
    }
    info!(raw_input_file = %raw_input_file.to_string_lossy(), "Reading raw input file");
    let gz = raw_input_file.to_string_lossy().ends_with(".gz");
    let delim = ctx.sheet.get_from_row(row, "column_delim");
    let file = std::fs::File::open(&raw_input_file).unwrap();
    let mut raw_data = if gz {
        let gz = flate2::read::GzDecoder::new(file);
        read_raw_data(delim, gz)
    } else {
        read_raw_data(delim, file)
    };
    debug!(header = ?raw_data.header, "Header");
    for col in ASSIGN_COL_NAMES.iter() {
        let val = ctx.sheet.get_from_row(row, col);
        if val != "NA" {
            for r in raw_data.header.iter_mut() {
                if r == val {
                    *r = col.to_string();
                }
            }
        }
    }
    debug!(header = ?raw_data.header, "Header");
    for chr in raw_data.col_mut("chr") {
        // a) Remove "chr" prefix
        if let Some(c) = chr.strip_prefix("chr") {
            *chr = c.to_string();
        }
        // b) Convert 23-25 to X, Y, M
        if *chr == "23" {
            *chr = "X".to_string();
        } else if *chr == "24" {
            *chr = "Y".to_string();
        } else if *chr == "25" {
            *chr = "M".to_string();
        }
    }
    // c) Change alleles to uppercase
    for r in raw_data.col_mut("ref") {
        *r = r.to_ascii_uppercase();
    }
    for a in raw_data.col_mut("alt") {
        *a = a.to_ascii_uppercase();
    }
    debug!(len = raw_data.data.len(), "Raw data before d and e");
    let data = std::mem::take(&mut raw_data.data);
    raw_data.data = data
        .into_par_iter()
        .filter(|x| {
            let r = raw_data.get_from_row(x.as_slice(), "ref");
            let a = raw_data.get_from_row(x.as_slice(), "alt");
            let effect_size = raw_data.get_from_row(x.as_slice(), "effect_size");
            // debug!(?x, r, a, effect_size, "Checking ref, alt, and effect size");
            // d) Remove SNPs with ambiguous ref or alt
            r != "I"
                && r != "D"
                && r != "IND"
                && r != "DEL"
                && a != "I"
                && a != "D"
                && a != "IND"
                && a != "DEL"
            // e) Remove variants with nonsensical effect estimates
                && effect_size != "Nan"
                && effect_size != "NaN"
                && effect_size != "NA"
                && effect_size != "Inf"
                && effect_size != "-Inf"
                && effect_size != "inf"
                && effect_size != "-inf"
        })
        .collect::<Vec<_>>();
    debug!(len = raw_data.data.len(), "Raw data after d and e");
    // f) Convert OR to beta
    let effect_is_or = ctx.sheet.get_from_row(row, "effect_is_OR");
    let effect_sizes = raw_data
        .col("effect_size")
        .map(|x| x.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    if effect_is_or == "N" && effect_sizes.iter().all(|x| *x > 0.0) {
        warn!(
            "All effect sizes are positive yet effect_is_OR has been set to N. Please double \
             check that effect estimates from the raw data file are indeed regression \
             coefficients and not odds ratios"
        );
    }
    if effect_is_or == "Y" && effect_sizes.iter().any(|x| *x < 0.0) {
        warn!(
            "Some effect sizes are negative yet effect_is_OR has been set to Y. Please double \
             check that effect estimates from the raw data file are indeed odds or hazard ratios \
             and not regression coefficients"
        );
    }
    if effect_is_or == "Y" {
        let data = std::mem::take(&mut raw_data.data);
        let effect_size = raw_data.idx("effect_size");
        raw_data.data = data
            .into_par_iter()
            .zip(effect_sizes)
            .filter_map(|(mut r, e)| {
                let l = e.ln();
                if l.is_nan() || l.is_infinite() {
                    None
                } else {
                    r[effect_size] = l.to_string();
                    Some(r)
                }
            })
            .collect::<Vec<_>>();
    }
    debug!(len = raw_data.data.len(), "Raw data after f");
    // g) Tabulate columns for sample sizes
    for var in ["total", "case", "ctrl"] {
        let var_col_name = ctx.sheet.get_from_row(row, &format!("N_{}_column", var));
        let var_value = ctx.sheet.get_from_row(row, &format!("N_{}", var));
        if var_col_name != "NA" {
            // rename column if values are present
            for r in raw_data.header.iter_mut() {
                if *r == format!("N_{}_column", var) {
                    *r = format!("N_{}", var);
                }
            }
        } else if var_value != "NA" {
            // update column
            for r in raw_data.col_mut(&format!("N_{}", var)) {
                r.clone_from(var_value);
            }
        }
    }
    let na = "NA".to_string();
    // if no sample sizes indicated and gwas legend input is NA then set all three
    // columns to NA
    debug!("g: Adding header");
    for var in ["total", "case", "ctrl"] {
        if !raw_data.header.contains(&format!("N_{}", var)) {
            raw_data.header.push(format!("N_{}", var));
        }
    }
    debug!("g: Added header");
    let header_len = raw_data.header.len();
    raw_data.data.par_iter_mut().for_each(|r| {
        let res = reserve_to(r, header_len);
        for _ in 0..res {
            r.push(na.clone());
        }
    });
    debug!("g: Added NAs");
    // compile case control or total sample sizes if inoformation is available
    let n_case = raw_data.idx("N_case");
    let n_ctrl = raw_data.idx("N_ctrl");
    let n_total = raw_data.idx("N_total");
    raw_data.data.par_iter_mut().for_each(|r| {
        if r[n_case] != "NA" && r[n_ctrl] != "NA" {
            r[n_total] =
                (r[n_case].parse::<f64>().unwrap() + r[n_ctrl].parse::<f64>().unwrap()).to_string();
        }
        if r[n_ctrl] != "NA" && r[n_total] != "NA" && r[n_case] == "NA" {
            r[n_case] = (r[n_total].parse::<f64>().unwrap() - r[n_ctrl].parse::<f64>().unwrap())
                .to_string();
        }
        if r[n_case] != "NA" && r[n_total] != "NA" && r[n_ctrl] == "NA" {
            r[n_ctrl] = (r[n_total].parse::<f64>().unwrap() - r[n_case].parse::<f64>().unwrap())
                .to_string();
        }
    });
    debug!(len = raw_data.data.len(), "Raw data after g");
    raw_data.reorder(&[
        "chr",
        "pos",
        "ref",
        "alt",
        "EAF",
        "effect_size",
        "standard_error",
        "pvalue",
        "pvalue_het",
        "N_total",
        "N_case",
        "N_ctrl",
    ]);
    let pos = raw_data.idx("pos");
    let chr = raw_data.idx("chr");
    let hg_version = ctx.sheet.get_from_row(row, "hg_version");
    raw_data.header[pos] = format!("pos_{}", hg_version);
    raw_data.header[chr] = format!("chr_{}", hg_version);
    debug!(header = ?raw_data.header, "Header");
    assert_eq!(raw_data.header.len(), raw_data.data[0].len());
    raw_data
}

#[tracing::instrument(skip(ctx, raw_data))]
fn liftover(ctx: &Ctx, raw_data: &Data) {
    let current_dir = std::env::current_dir().unwrap();
    let liftover_dir = std::path::Path::new(&ctx.args.liftover_dir);
    let mut bed = std::fs::File::create(current_dir.join("input.bed")).unwrap();
    let pos_hg17 = raw_data.header.contains(&"pos_hg17".to_string());
    let pos_hg18 = raw_data.header.contains(&"pos_hg18".to_string());
    let pos_hg19 = raw_data.header.contains(&"pos_hg19".to_string());
    let pos_hg38 = raw_data.header.contains(&"pos_hg38".to_string());
    debug!(
        pos_hg17,
        pos_hg18, pos_hg19, pos_hg38, "Checking position columns"
    );
    if pos_hg17 || pos_hg18 || pos_hg19 || pos_hg38 {
        let chr_idx = raw_data.idx(if pos_hg17 {
            "chr_hg17"
        } else if pos_hg18 {
            "chr_hg18"
        } else if pos_hg19 {
            "chr_hg19"
        } else {
            "chr_hg38"
        });
        let pos_idx = raw_data.idx(if pos_hg17 {
            "pos_hg17"
        } else if pos_hg18 {
            "pos_hg18"
        } else if pos_hg19 {
            "pos_hg19"
        } else {
            "pos_hg38"
        });
        for (i, r) in raw_data.data.iter().enumerate() {
            writeln!(
                bed,
                "chr{}\t{}\t{}\t{}",
                r[chr_idx],
                r[pos_idx].parse::<i64>().unwrap() - 1,
                r[pos_idx],
                i + 2
            )
            .unwrap();
        }
        drop(bed);
        if pos_hg17 || pos_hg18 {
            std::process::Command::new(&ctx.args.liftover)
                .arg(current_dir.join("input.bed"))
                .arg(liftover_dir.join(if pos_hg17 {
                    "hg17ToHg19.over.chain.gz"
                } else {
                    "hg18ToHg19.over.chain.gz"
                }))
                .arg(current_dir.join("input2.bed"))
                .arg(current_dir.join("1unlifted.bed"))
                .status()
                .unwrap();
            let mut hg19 = std::fs::File::create(current_dir.join("hg19.bed")).unwrap();
            for line in std::fs::read_to_string(current_dir.join("input2.bed"))
                .unwrap()
                .lines()
            {
                writeln!(hg19, "{}", line.strip_prefix("chr").unwrap_or(line)).unwrap();
            }
        } else {
            std::fs::rename(
                current_dir.join("input.bed"),
                current_dir.join("input2.bed"),
            )
            .unwrap();
        }
        std::process::Command::new(&ctx.args.liftover)
            .arg(current_dir.join("input2.bed"))
            .arg(liftover_dir.join(if pos_hg38 {
                "hg38ToHg19.over.chain.gz"
            } else {
                "hg19ToHg38.over.chain.gz"
            }))
            .arg(current_dir.join("final.bed"))
            .arg(current_dir.join("unlifted.bed"))
            .status()
            .unwrap();
        let hg38_input = if pos_hg38 { "input2.bed" } else { "final.bed" };
        debug!(hg38_input, "Reading hg38 bed file");
        let mut hg38 = std::fs::File::create(current_dir.join("hg38.bed")).unwrap();
        for line in std::fs::read_to_string(current_dir.join(hg38_input))
            .unwrap()
            .lines()
        {
            writeln!(hg38, "{}", line.strip_prefix("chr").unwrap_or(line)).unwrap();
        }
        std::fs::remove_file(current_dir.join(hg38_input)).unwrap();
        if pos_hg19 || pos_hg38 {
            let hg19_input = if pos_hg38 { "final.bed" } else { "input2.bed" };
            debug!(hg19_input, "Reading hg19 bed file");
            let mut hg19 = std::fs::File::create(current_dir.join("hg19.bed")).unwrap();
            for line in std::fs::read_to_string(current_dir.join(hg19_input))
                .unwrap()
                .lines()
            {
                writeln!(hg19, "{}", line.strip_prefix("chr").unwrap_or(line)).unwrap();
            }
            std::fs::remove_file(current_dir.join(hg19_input)).unwrap();
        }
    } else {
        error!("No position columns found in the raw data file");
        panic!();
    }
}

#[tracing::instrument(skip(ctx, raw_data))]
fn dbsnp_matching(ctx: &Ctx, mut raw_data: Data) -> (Data, Data) {
    debug!("Reading hg19 and hg38 bed files");
    let hg19 = {
        if raw_data.header.contains(&"chr_hg19".to_string()) {
            None
        } else {
            raw_data.header.push("chr_hg19".to_string());
            raw_data.header.push("pos_hg19".to_string());
            let file =
                std::fs::File::open(std::env::current_dir().unwrap().join("hg19.bed")).unwrap();
            Some(
                Data::read('\t', file, false)
                    .data
                    .into_iter()
                    .map(|x| (x.get(3).unwrap().parse::<usize>().unwrap() - 2, x))
                    .collect::<HashMap<usize, _>>(),
            )
        }
    };
    let hg38 = {
        if raw_data.header.contains(&"chr_hg38".to_string()) {
            None
        } else {
            raw_data.header.push("chr_hg38".to_string());
            raw_data.header.push("pos_hg38".to_string());
            let file =
                std::fs::File::open(std::env::current_dir().unwrap().join("hg38.bed")).unwrap();
            Some(
                Data::read('\t', file, false)
                    .data
                    .into_iter()
                    .map(|x| (x.get(3).unwrap().parse::<usize>().unwrap() - 2, x))
                    .collect::<HashMap<usize, _>>(),
            )
        }
    };
    debug!(
        raw_data = raw_data.data.len(),
        "Read hg19 and hg38 bed files"
    );
    let header_len = raw_data.header.len();
    raw_data
        .data
        .par_iter_mut()
        .enumerate()
        .for_each(move |(i, r)| {
            reserve_to(r, header_len);
            if let Some(ref hg19) = hg19 {
                let hg19 = hg19.get(&i);
                if let Some(hg19) = hg19 {
                    r.push(hg19.first().unwrap().to_string());
                    r.push(hg19.get(2).unwrap().to_string());
                } else {
                    r.push("NA".to_string());
                    r.push("NA".to_string());
                }
            }
            if let Some(ref hg38) = hg38 {
                let hg38 = hg38.get(&i);
                if let Some(hg38) = hg38 {
                    r.push(hg38.first().unwrap().to_string());
                    r.push(hg38.get(2).unwrap().to_string());
                } else {
                    r.push("NA".to_string());
                    r.push("NA".to_string());
                }
            }
        });

    debug!("Reordering columns");
    raw_data.reorder(&[
        "chr_hg19",
        "pos_hg19",
        "ref",
        "alt",
        "effect_size",
        "standard_error",
        "EAF",
        "pvalue",
        "pvalue_het",
        "N_total",
        "N_case",
        "N_ctrl",
        "chr_hg38",
        "pos_hg38",
    ]);
    // raw_data.write("dbsnp.e.txt.gz");
    debug!(len = raw_data.data.len(), "Raw data after bed matching");

    debug!("Reading dbSNP file");
    let dbsnp = flate2::read::GzDecoder::new(std::fs::File::open(&ctx.args.dbsnp_file).unwrap());
    let dbsnp = Data::read('\t', dbsnp, true);
    debug!("Merging dbSNP data");
    let dbsnp_idxs = [
        dbsnp.idx("chr"),
        dbsnp.idx("pos_hg19"),
        dbsnp.idx("ref"),
        dbsnp.idx("alt"),
        dbsnp.idx("pos_hg38"),
    ];
    debug!("Creating dbsnp map");
    let dbsnp_map: HashMap<(&str, &str, &str, &str, &str), &Vec<String>> =
        HashMap::from_par_iter(dbsnp.data.par_iter().map(|x| {
            (
                (
                    x[dbsnp_idxs[0]].as_str(),
                    x[dbsnp_idxs[1]].as_str(),
                    x[dbsnp_idxs[2]].as_str(),
                    x[dbsnp_idxs[3]].as_str(),
                    x[dbsnp_idxs[4]].as_str(),
                ),
                x,
            )
        }));
    debug!("Getting raw data indexes");
    let raw_data_idxs = [
        raw_data.idx("chr_hg19"),
        raw_data.idx("pos_hg19"),
        raw_data.idx("ref"),
        raw_data.idx("alt"),
        raw_data.idx("pos_hg38"),
    ];
    let raw_data_merged_flipped_idxs = [
        raw_data.idx("chr_hg19"),
        raw_data.idx("pos_hg19"),
        raw_data.idx("alt"),
        raw_data.idx("ref"),
        raw_data.idx("pos_hg38"),
    ];
    let mut raw_data_merged = raw_data.clone();
    let raw_data_merged_data = std::mem::take(&mut raw_data_merged.data);
    for i in 0..dbsnp.header.len() {
        if !dbsnp_idxs.contains(&i) {
            debug!(i, header = dbsnp.header[i], "Adding missing column");
            raw_data_merged.header.push(dbsnp.header[i].clone());
        }
    }
    raw_data_merged.header.push("unique_id".to_string());
    let unique_id_idx = raw_data_merged.idx("unique_id");
    let mut raw_data_flipped = raw_data_merged.clone();
    debug!(header = ?raw_data_merged.header, "Header");
    debug!(idxs = ?raw_data_idxs, "Raw data indexes");
    let header_len = raw_data_merged.header.len();
    raw_data_merged.data = raw_data_merged_data
        .into_par_iter()
        .filter_map(|mut r| {
            reserve_to(&mut r, header_len);
            let key = (
                r[raw_data_idxs[0]].as_str(),
                r[raw_data_idxs[1]].as_str(),
                r[raw_data_idxs[2]].as_str(),
                r[raw_data_idxs[3]].as_str(),
                r[raw_data_idxs[4]].as_str(),
            );
            let dbsnp_data = *dbsnp_map.get(&key)?;
            (0..dbsnp.header.len()).for_each(|i| {
                if !dbsnp_idxs.contains(&i) {
                    r.push(dbsnp_data[i].clone());
                }
            });
            r.push(format!(
                "{}_{}_{}_{}",
                r[raw_data_idxs[0]], r[raw_data_idxs[1]], r[raw_data_idxs[2]], r[raw_data_idxs[3]],
            ));
            Some(r)
        })
        .collect::<Vec<_>>();
    debug!("Flipping alleles");
    let mut raw_data_flipped_data = std::mem::take(&mut raw_data_flipped.data);
    let header_len = raw_data_flipped.header.len();
    raw_data_flipped_data = raw_data_flipped_data
        .into_par_iter()
        .filter_map(|mut r| {
            reserve_to(&mut r, header_len);
            let key = (
                r[raw_data_merged_flipped_idxs[0]].as_str(),
                r[raw_data_merged_flipped_idxs[1]].as_str(),
                r[raw_data_merged_flipped_idxs[2]].as_str(),
                r[raw_data_merged_flipped_idxs[3]].as_str(),
                r[raw_data_merged_flipped_idxs[4]].as_str(),
            );
            let dbsnp_data = *dbsnp_map.get(&key)?;
            (0..dbsnp.header.len()).for_each(|i| {
                if !dbsnp_idxs.contains(&i) {
                    r.push(dbsnp_data[i].clone());
                }
            });
            r.push(format!(
                "{}_{}_{}_{}",
                r[raw_data_idxs[0]], r[raw_data_idxs[1]], r[raw_data_idxs[2]], r[raw_data_idxs[3]],
            ));
            Some(r)
        })
        .collect::<Vec<_>>();
    debug!("Merging flipped alleles");
    let unique_ids: HashSet<&str> = HashSet::from_iter(
        raw_data_merged
            .data
            .iter()
            .map(|x| x[unique_id_idx].as_str()),
    );
    raw_data_flipped.data = raw_data_flipped_data
        .into_par_iter()
        .filter(|x| !unique_ids.contains(x[unique_id_idx].as_str()))
        .collect::<Vec<_>>();
    let alt = raw_data_flipped.idx("alt");
    let ref_ = raw_data_flipped.idx("ref");
    let effect_size = raw_data_flipped.idx("effect_size");
    let eaf = raw_data_flipped.idx("EAF");
    raw_data_flipped.data.par_iter_mut().for_each(|r| {
        let (one, two) = r.split_at_mut(alt.max(ref_));
        let min = alt.min(ref_);
        let max = alt.max(ref_);
        std::mem::swap(&mut one[min], &mut two[max]);
        let es = r[effect_size].parse::<f64>().unwrap();
        r[effect_size] = (-es).to_string();
        let e = r[eaf].parse::<f64>().unwrap();
        r[eaf] = (1.0 - e).to_string();
        let unique_id = r.len() - 1;
        r[unique_id] = format!(
            "{}_{}_{}_{}",
            r[raw_data_idxs[0]], r[raw_data_idxs[1]], r[raw_data_idxs[2]], r[raw_data_idxs[3]]
        );
    });
    raw_data_merged.data.extend(raw_data_flipped.data);
    let mut seen = HashSet::new();
    raw_data_merged
        .data
        .retain(|x| seen.insert(x[unique_id_idx].as_str().to_string()));
    debug!("Merging missing data");
    let new_order = [
        "rsid",
        "unique_id",
        "chr_hg19",
        "pos_hg19",
        "ref",
        "alt",
        "effect_size",
        "standard_error",
        "EAF",
        "pvalue",
        "pvalue_het",
        "N_total",
        "N_case",
        "N_ctrl",
        "chr_hg38",
        "pos_hg38",
        "gnomAD_AF_EUR",
        "gnomAD_AF_AMR",
        "gnomAD_AF_AFR",
        "gnomAD_AF_EAS",
        "gnomAD_AF_SAS",
    ];
    debug!("Constructing raw unique ids");
    let raw_unique_ids: HashSet<(&str, &str, &str, &str)> = HashSet::from_par_iter(
        raw_data_merged
            .data
            .par_iter()
            .map(|r| {
                (
                    r[raw_data_idxs[0]].as_str(),
                    r[raw_data_idxs[1]].as_str(),
                    r[raw_data_idxs[2]].as_str(),
                    r[raw_data_idxs[3]].as_str(),
                )
            })
            .chain(raw_data_merged.data.par_iter().map(|r| {
                (
                    r[raw_data_idxs[0]].as_str(),
                    r[raw_data_idxs[1]].as_str(),
                    r[raw_data_idxs[3]].as_str(),
                    r[raw_data_idxs[2]].as_str(),
                )
            })),
    );
    let pos_hg19 = raw_data.idx("pos_hg19");
    let pos_hg38 = raw_data.idx("pos_hg38");
    debug!("Constructing missing data");
    let header = raw_data.header.clone();
    let raw_data_missing = raw_data
        .data
        .into_par_iter()
        .filter(|r| {
            !raw_unique_ids.contains(&(
                r[raw_data_idxs[0]].as_str(),
                r[raw_data_idxs[1]].as_str(),
                r[raw_data_idxs[2]].as_str(),
                r[raw_data_idxs[3]].as_str(),
            )) && r[pos_hg19] != "NA"
                && r[pos_hg38] != "NA"
                && r[pos_hg19] != "NaN"
                && r[pos_hg38] != "NaN"
        })
        .collect::<Vec<_>>();
    let mut raw_data_missing = Data {
        header,
        data: raw_data_missing,
    };
    debug!(
        header = ?raw_data.header,
        len = raw_data.header.len(),
        "Raw data header"
    );
    debug!(
        header = ?raw_data_merged.header,
        len = raw_data_merged.header.len(),
        "Merged data header"
    );
    debug!(
        header = ?raw_data_missing.header,
        len = raw_data_missing.header.len(),
        "Missing data header"
    );
    debug!("Reordering columns");
    raw_data_merged.reorder(&new_order);
    for i in 0..dbsnp.header.len() {
        if !dbsnp_idxs.contains(&i) {
            debug!(i, header = dbsnp.header[i], "Adding missing column");
            raw_data_missing.header.push(dbsnp.header[i].clone());
        }
    }
    raw_data_missing.header.push("unique_id".to_string());
    let header_len = raw_data_missing.header.len();
    raw_data_missing.data.par_iter_mut().for_each(|r| {
        reserve_to(r, header_len);
        for i in 0..dbsnp.header.len() {
            if !dbsnp_idxs.contains(&i) {
                r.push("NA".to_string());
            }
        }
        r.push(format!(
            "{}_{}_{}_{}",
            r[raw_data_idxs[0]], r[raw_data_idxs[1]], r[raw_data_idxs[2]], r[raw_data_idxs[3]]
        ));
    });
    debug!(header = ?raw_data_missing.header);
    assert_eq!(
        raw_data_missing.header.len(),
        raw_data_missing.data[0].len()
    );
    raw_data_missing.reorder(&new_order);
    debug!(header = ?raw_data_merged.header);

    assert_eq!(raw_data_merged.header.len(), raw_data_merged.data[0].len());
    debug!(header = ?raw_data_missing.header);
    assert_eq!(
        raw_data_missing.header.len(),
        raw_data_missing.data[0].len()
    );
    (raw_data_merged, raw_data_missing)
}

#[tracing::instrument(skip(ctx, raw_data_merged, raw_data_missing))]
fn ref_alt_check(ctx: &Ctx, mut raw_data_merged: Data, raw_data_missing: Data) -> Data {
    let chr_hg38 = raw_data_missing.idx("chr_hg38");
    let pos_hg38 = raw_data_missing.idx("pos_hg38");
    let inputs = raw_data_missing
        .data
        .iter()
        .map(|r| format!("chr{}:{}-{}", r[chr_hg38], r[pos_hg38], r[pos_hg38]))
        .collect::<Vec<_>>();
    let num_inputs = inputs.len();
    let num_threads = ctx
        .args
        .samtools_threads
        .unwrap_or_else(|| num_cpus::get() * 4);
    let nucleotides = Mutex::new(Vec::with_capacity(num_inputs));
    nucleotides
        .lock()
        .unwrap()
        .extend((0..num_inputs).map(|_| MaybeUninit::uninit()));
    let chunk_size = 5000;
    let chunks = (num_inputs + chunk_size - 1) / chunk_size;
    let chunks = Mutex::new((0..chunks).collect::<Vec<_>>());
    debug!(
        num_threads,
        num_inputs,
        chunk_size,
        chunks = chunks.lock().unwrap().len(),
        "Running samtools"
    );
    std::thread::scope(|s| {
        for _ in 0..num_threads {
            s.spawn(|| {
                loop {
                    let chunk = {
                        let mut chunks = chunks.lock().unwrap();
                        if chunks.is_empty() {
                            return;
                        }
                        chunks.pop().unwrap()
                    };
                    let j = chunk * chunk_size;
                    let end = (j + chunk_size).min(num_inputs);
                    let input = &inputs[j..end];
                    debug!(chunk, "Got input");
                    let mut cmd = std::process::Command::new(&ctx.args.samtools);
                    cmd.arg("faidx");
                    cmd.arg(&ctx.args.fasta_ref);
                    for i in input {
                        cmd.arg(i);
                    }
                    debug!(chunk, "Constructed samtools command");
                    let output = match cmd.output() {
                        Ok(o) => o,
                        Err(e) if e.kind() == std::io::ErrorKind::OutOfMemory => {
                            chunks.lock().unwrap().push(chunk);
                            return;
                        },
                        Err(e) => {
                            error!(chunk, ?e, "Failed to run samtools");
                            return;
                        },
                    };
                    debug!(chunk, "Ran samtools");
                    let output = String::from_utf8(output.stdout).unwrap();
                    let mut nucleotides = nucleotides.lock().unwrap();
                    for (idx, l) in output.lines().filter(|x| !x.starts_with('>')).enumerate() {
                        nucleotides[idx + j].write(if l.len() > 1 {
                            "N".to_string()
                        } else {
                            l.to_uppercase()
                        });
                    }
                    debug!(chunk, "Finished samtools");
                }
            });
        }
    });
    debug!("Finished samtools");
    let nucleotides: Vec<String> =
        unsafe { std::mem::transmute(nucleotides.into_inner().unwrap()) };
    debug!("Flattened nucleotides");
    // let mut file = std::fs::File::create("nucleotides.txt.gz").unwrap();
    // for n in &nucleotides {
    //     writeln!(file, "{n}").unwrap();
    // }
    // drop(file);
    let ref_ = raw_data_merged.idx("ref");
    let alt = raw_data_merged.idx("alt");
    let effect_size = raw_data_merged.idx("effect_size");
    let eaf = raw_data_merged.idx("EAF");
    raw_data_merged.data.par_extend(
        raw_data_missing
            .data
            .into_par_iter()
            .zip(nucleotides)
            .filter_map(|(mut d, n)| {
                if d[alt] == n {
                    let (one, two) = d.split_at_mut(alt.max(ref_));
                    let min = alt.min(ref_);
                    let max = alt.max(ref_) - one.len();
                    std::mem::swap(&mut one[min], &mut two[max]);
                    let es = d[effect_size].parse::<f64>().unwrap();
                    d[effect_size] = (-es).to_string();
                    if d[eaf] != "NA" && d[eaf] != "NaN" {
                        let e = d[eaf].parse::<f64>().unwrap();
                        d[eaf] = (1.0 - e).to_string();
                    }
                    Some(d)
                } else if d[ref_] == n {
                    Some(d)
                } else {
                    None
                }
            }),
    );
    debug!("Merged missing data");
    raw_data_merged
}

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::builder()
                .with_default_directive(tracing::Level::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let args = Args::parse();
    if args.google_sheets_id.starts_with("http") {
        error!("google_sheets_id should be the ID of the Google Sheets document, not the URL. For example, if the URL is https://docs.google.com/spreadsheets/d/1a2b3c4d5e6f7g8h9i0j1k2l3m4n5o6p7q8r9s0t1u2v3w4x5y6z7/edit#gid=0, the ID is 1a2b3c4d5e6f7g8h9i0j1k2l3m4n5o6p7q8r9s0t1u2v3w4x5y6z7");
        return;
    }
    let spreadsheet = reqwest::blocking::get(format!(
        "https://sheets.googleapis.com/v4/spreadsheets/{}?key={}",
        args.google_sheets_id, GOOGLE_SHEETS_API_KEY
    ))
    .unwrap()
    .error_for_status()
    .unwrap();
    let spreadsheet = spreadsheet.text().unwrap();
    let spreadsheet: serde_json::Value = serde_json::from_str(&spreadsheet).unwrap();
    let spreadsheet = spreadsheet["sheets"].as_array().unwrap()[0]["properties"]["title"]
        .as_str()
        .unwrap();
    let data = reqwest::blocking::get(format!(
        "https://sheets.googleapis.com/v4/spreadsheets/{}/values/{}?key={}",
        args.google_sheets_id, spreadsheet, GOOGLE_SHEETS_API_KEY
    ))
    .unwrap()
    .error_for_status()
    .unwrap();
    let data = data.text().unwrap();
    let data: serde_json::Value = serde_json::from_str(&data).unwrap();
    let data = data["values"].as_array().unwrap();
    let header = data[0].as_array().unwrap();
    let header = header
        .iter()
        .map(|x| x.as_str().unwrap().to_string())
        .collect::<Vec<_>>();
    let data = data[1..]
        .iter()
        .map(|x| {
            x.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_str().unwrap().to_string())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let data = Data { header, data };
    debug!("Header: {:?}", data.header);
    let ctx = Ctx { args, sheet: data };
    info!(trait_name = %ctx.args.trait_name, "Starting pipeline");
    info!("Starting preformatting");
    let raw_data = preformat(&ctx);
    // raw_data.write("raw_data.txt.gz");
    info!("Starting liftover");
    liftover(&ctx, &raw_data);
    info!("Starting dbSNP matching");
    let (raw_data_merged, raw_data_missing) = dbsnp_matching(&ctx, raw_data);
    // raw_data_merged.write("raw_data_merged.txt.gz");
    // raw_data_missing.write("raw_data_missing.txt.gz");
    info!("Starting ref/alt check");
    let final_data = ref_alt_check(&ctx, raw_data_merged, raw_data_missing);
    info!("Writing final data to {}", ctx.args.output_file);
    final_data.write(&ctx.args.output_file);
    info!("Pipeline complete");
}
