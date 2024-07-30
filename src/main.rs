use clap::Parser;
use itertools::Itertools;
use rayon::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    io::Write,
    mem::MaybeUninit,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Mutex,
    },
};
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
    trait_name: String,
    #[arg(short = 'i', long)]
    raw_input_dir: String,
    #[arg(short, long)]
    liftover: String,
    #[arg(long)]
    liftover_dir: String,
    #[arg(short = 'r', long)]
    grs_dir: String,
    #[arg(short, long)]
    dbsnp_file: String,
    #[arg(short, long)]
    samtools: String,
    #[arg(short, long)]
    fasta_ref: String,
    #[arg(short, long)]
    output_file: String,
}

pub struct Ctx {
    args: Args,
    sheet: Data,
}

#[derive(Clone)]
pub struct Data {
    header: Vec<String>,
    data: Vec<Vec<String>>,
}

impl Data {
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
        let idx = self.idx(key);
        debug!(key, idx, "Mutating column");
        self.data.iter_mut().map(move |x| &mut x[idx])
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
    let delim = ctx.sheet.get_from_row(row, "column_delim");
    let mut contents = if delim == "\t" || delim == "tab" {
        csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(raw_input_file)
    } else if delim == "," || delim == "comma" {
        csv::ReaderBuilder::new()
            .delimiter(b',')
            .has_headers(true)
            .from_path(raw_input_file)
    } else if delim == "space" {
        csv::ReaderBuilder::new()
            .delimiter(b' ')
            .has_headers(true)
            .from_path(raw_input_file)
    } else {
        error!("Invalid column delimiter {}", delim);
        panic!();
    }
    .unwrap();
    let header = contents
        .headers()
        .unwrap()
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>();
    if header.len() <= 4 {
        error!("Raw input file has less than 5 columns, likely the column delimiter has been misspecified");
        panic!();
    }
    let data = contents.records().map(|x| x.unwrap()).collect::<Vec<_>>();
    let data = data
        .iter()
        .map(|x| x.iter().map(|x| x.to_string()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let na = "NA".to_string();
    let mut raw_data = Data { header, data };
    debug!(header = ?raw_data.header, "Header");
    for col in ASSIGN_COL_NAMES.iter() {
        let idx = raw_data.idx_opt(col);
        if idx.is_none() {
            debug!(col, "Assigning NA to missing column");
            raw_data.header.push(col.to_string());
            for r in raw_data.data.iter_mut() {
                r.push(na.clone());
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
        .into_iter()
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
        warn!("All effect sizes are positive yet effect_is_OR has been set to N. Please double check that effect estimates from the raw data file are indeed regression coefficients and not odds ratios");
    }
    if effect_is_or == "Y" && effect_sizes.iter().any(|x| *x < 0.0) {
        warn!("Some effect sizes are negative yet effect_is_OR has been set to Y. Please double check that effect estimates from the raw data file are indeed odds or hazard ratios and not regression coefficients");
    }
    if effect_is_or == "Y" {
        let data = std::mem::take(&mut raw_data.data);
        let effect_size = raw_data.idx("effect_size");
        raw_data.data = data
            .into_iter()
            .zip(effect_sizes)
            .filter_map(|(mut r, e)| {
                let l = e.ln();
                if l.is_nan() {
                    None
                } else {
                    r[effect_size] = l.to_string();
                    Some(r)
                }
            })
            .collect::<Vec<_>>();
    }
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
    // if no sample sizes indicated and gwas legend input is NA then set all three columns to NA
    for var in ["total", "case", "ctrl"] {
        if !raw_data.header.contains(&format!("N_{}", var)) {
            raw_data.header.push(format!("N_{}", var));
            for r in raw_data.data.iter_mut() {
                r.push(na.clone());
            }
        }
    }
    // compile case control or total sample sizes if inoformation is available
    let n_case = raw_data.idx("N_case");
    let n_ctrl = raw_data.idx("N_ctrl");
    let n_total = raw_data.idx("N_total");
    for r in raw_data.data.iter_mut() {
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
    }
    let pos = raw_data.idx("pos");
    let chr = raw_data.idx("chr");
    let hg_version = ctx.sheet.get_from_row(row, "hg_version");
    raw_data.header[pos] = format!("pos_{}", hg_version);
    raw_data.header[chr] = format!("chr_{}", hg_version);
    debug!(header = ?raw_data.header, "Header");
    raw_data
}

#[tracing::instrument(skip(ctx, raw_data))]
fn liftover(ctx: &Ctx, raw_data: &Data) {
    let liftover_dir = std::path::Path::new(&ctx.args.liftover_dir);
    let mut bed = std::fs::File::create(liftover_dir.join("input.bed")).unwrap();
    let pos_hg17 = raw_data.header.contains(&"pos_hg17".to_string());
    let pos_hg18 = raw_data.header.contains(&"pos_hg18".to_string());
    let pos_hg19 = raw_data.header.contains(&"pos_hg19".to_string());
    let pos_hg38 = raw_data.header.contains(&"pos_hg38".to_string());
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
                .arg(liftover_dir.join("input.bed"))
                .arg(liftover_dir.join(if pos_hg17 {
                    "hg17ToHg19.over.chain.gz"
                } else {
                    "hg18ToHg19.over.chain.gz"
                }))
                .arg(liftover_dir.join("input2.bed"))
                .arg(liftover_dir.join("1unlifted.bed"))
                .status()
                .unwrap();
            let mut hg19 = std::fs::File::create(liftover_dir.join("hg19.bed")).unwrap();
            for line in std::fs::read_to_string(liftover_dir.join("input2.bed"))
                .unwrap()
                .lines()
            {
                writeln!(hg19, "{}", line.strip_prefix("chr").unwrap_or(line)).unwrap();
            }
        } else {
            std::fs::rename(
                liftover_dir.join("input.bed"),
                liftover_dir.join("input2.bed"),
            )
            .unwrap();
        }
        std::process::Command::new(&ctx.args.liftover)
            .arg(liftover_dir.join("input2.bed"))
            .arg(liftover_dir.join(if pos_hg38 {
                "hg38ToHg19.over.chain.gz"
            } else {
                "hg19ToHg38.over.chain.gz"
            }))
            .arg(liftover_dir.join("final.bed"))
            .arg(liftover_dir.join("unlifted.bed"))
            .status()
            .unwrap();
        let hg38_input = if pos_hg38 { "input2.bed" } else { "final.bed" };
        let mut hg38 = std::fs::File::create(liftover_dir.join("hg38.bed")).unwrap();
        for line in std::fs::read_to_string(liftover_dir.join(hg38_input))
            .unwrap()
            .lines()
        {
            writeln!(hg38, "{}", line.strip_prefix("chr").unwrap_or(line)).unwrap();
        }
        std::fs::remove_file(liftover_dir.join(hg38_input)).unwrap();
        if pos_hg19 || pos_hg38 {
            let hg19_input = if pos_hg38 { "final.bed" } else { "input2.bed" };
            let mut hg19 = std::fs::File::create(liftover_dir.join("hg19.bed")).unwrap();
            for line in std::fs::read_to_string(liftover_dir.join(hg19_input))
                .unwrap()
                .lines()
            {
                writeln!(hg19, "{}", line.strip_prefix("chr").unwrap_or(line)).unwrap();
            }
            std::fs::remove_file(liftover_dir.join(hg19_input)).unwrap();
        }
    } else {
        error!("No position columns found in the raw data file");
        panic!();
    }
}

#[tracing::instrument(skip(ctx, raw_data))]
fn dbsnp_matching(ctx: &Ctx, raw_data: &mut Data) -> (Data, Data) {
    debug!("Reading hg19 and hg38 bed files");
    let mut hg19 = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(std::env::current_dir().unwrap().join("hg19.bed"))
        .unwrap();
    let hg19 = hg19.records().map(|x| x.unwrap()).collect::<Vec<_>>();
    let mut hg38 = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(std::env::current_dir().unwrap().join("hg38.bed"))
        .unwrap();
    let hg38 = hg38.records().map(|x| x.unwrap()).collect::<Vec<_>>();
    debug!(
        hg19 = hg19.len(),
        hg38 = hg38.len(),
        raw_data = raw_data.data.len(),
        "Read hg19 and hg38 bed files"
    );
    raw_data.header.push("chr_hg19".to_string());
    raw_data.header.push("pos_hg19".to_string());
    raw_data.header.push("chr_hg38".to_string());
    raw_data.header.push("pos_hg38".to_string());
    raw_data
        .data
        .par_iter_mut()
        .zip(hg19.par_iter())
        .zip(hg38.par_iter())
        .for_each(|((r, hg19), hg38)| {
            let chr_hg19 = hg19.get(0).unwrap().to_string();
            let pos_hg19 = hg19.get(2).unwrap().to_string();
            let chr_hg38 = hg38.get(0).unwrap().to_string();
            let pos_hg38 = hg38.get(2).unwrap().to_string();
            r.push(chr_hg19);
            r.push(pos_hg19);
            r.push(chr_hg38);
            r.push(pos_hg38);
        });

    debug!("Reordering columns");
    let new_headers = [
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
    ];
    let new_order = new_headers
        .iter()
        .map(|x| raw_data.idx(x))
        .collect::<Vec<_>>();
    let nrows = raw_data.data.len();
    let data = std::mem::take(&mut raw_data.data);
    let new_data: Vec<MaybeUninit<Vec<String>>> =
        (0..nrows).map(|_| MaybeUninit::uninit()).collect();
    data.into_par_iter().enumerate().for_each(|(i, r)| {
        let new_r = r
            .into_iter()
            .enumerate()
            .filter(|(i, _)| new_order.contains(i))
            .sorted_by_key(|(i, _)| new_order.iter().position(|x| x == i))
            .map(|(_, x)| x)
            .collect::<Vec<_>>();
        unsafe { &mut *new_data.as_ptr().add(i).cast_mut() }.write(new_r);
    });
    raw_data.header = new_headers
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>();
    raw_data.data =
        unsafe { std::mem::transmute::<Vec<MaybeUninit<Vec<String>>>, Vec<Vec<String>>>(new_data) };
    debug!(len = raw_data.data.len(), "Raw data after bed matching");

    debug!("Reading dbSNP file");
    let dbsnp = flate2::read::GzDecoder::new(std::fs::File::open(&ctx.args.dbsnp_file).unwrap());
    let mut dbsnp = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(dbsnp);
    let header = dbsnp
        .headers()
        .unwrap()
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>();
    let data = dbsnp
        .records()
        .par_bridge()
        .map(|x| x.unwrap().iter().map(|x| x.to_string()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let dbsnp = Data { header, data };
    debug!("Merging dbSNP data");
    let dbsnp_idxs = [
        dbsnp.idx("chr"),
        dbsnp.idx("pos_hg19"),
        dbsnp.idx("ref"),
        dbsnp.idx("alt"),
        dbsnp.idx("pos_hg38"),
    ];
    let dbsnp_map: HashMap<(&str, &str, &str, &str, &str), &Vec<String>> =
        HashMap::from_iter(dbsnp.data.iter().map(|x| {
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
    let raw_data_data = std::mem::take(&mut raw_data_merged.data);
    for i in 0..dbsnp.header.len() {
        if !dbsnp_idxs.contains(&i) {
            raw_data_merged.header.push(dbsnp.header[i].clone());
        }
    }
    raw_data_merged.header.push("unique_id".to_string());
    let unique_id_idx = raw_data_merged.idx("unique_id");
    let mut raw_data_flipped = raw_data_merged.clone();
    raw_data_merged.data = raw_data_data
        .into_par_iter()
        .filter_map(|mut r| {
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
    raw_data_flipped_data = raw_data_flipped_data
        .into_par_iter()
        .filter_map(|mut r| {
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
        .into_iter()
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
    let raw_unique_ids: HashSet<(&str, &str, &str, &str)> = HashSet::from_iter(
        raw_data_merged
            .data
            .iter()
            .map(|r| {
                (
                    r[raw_data_idxs[0]].as_str(),
                    r[raw_data_idxs[1]].as_str(),
                    r[raw_data_idxs[2]].as_str(),
                    r[raw_data_idxs[3]].as_str(),
                )
            })
            .chain(raw_data_merged.data.iter().map(|r| {
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
    let raw_data_missing = raw_data
        .data
        .iter()
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
        .cloned()
        .collect::<Vec<_>>();
    let mut raw_data_missing = Data {
        header: raw_data.header.clone(),
        data: raw_data_missing,
    };
    let data = std::mem::take(&mut raw_data_merged.data);
    let new_data: Vec<MaybeUninit<Vec<String>>> =
        (0..data.len()).map(|_| MaybeUninit::uninit()).collect();
    data.into_par_iter().enumerate().for_each(|(i, r)| {
        let new_r = r
            .into_iter()
            .enumerate()
            .filter(|(i, _)| new_order.contains(&raw_data_merged.header[*i].as_str()))
            .sorted_by_key(|(i, _)| {
                new_order
                    .iter()
                    .position(|x| x == &raw_data_merged.header[*i])
            })
            .map(|(_, x)| x)
            .collect::<Vec<_>>();
        unsafe { &mut *new_data.as_ptr().add(i).cast_mut() }.write(new_r);
    });
    raw_data_merged.header = new_order.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    raw_data_merged.data =
        unsafe { std::mem::transmute::<Vec<MaybeUninit<Vec<String>>>, Vec<Vec<String>>>(new_data) };
    debug!("Reordering columns");
    for i in new_order {
        if !raw_data_merged.header.contains(&i.to_string()) {
            raw_data_merged.header.push(i.to_string());
            for r in raw_data_merged.data.iter_mut() {
                r.push("NA".to_string());
            }
        }
    }
    let data = std::mem::take(&mut raw_data_missing.data);
    let new_data: Vec<MaybeUninit<Vec<String>>> =
        (0..data.len()).map(|_| MaybeUninit::uninit()).collect();
    data.into_par_iter().enumerate().for_each(|(i, r)| {
        let new_r = r
            .into_iter()
            .enumerate()
            .filter(|(i, _)| new_order.contains(&raw_data_missing.header[*i].as_str()))
            .sorted_by_key(|(i, _)| {
                new_order
                    .iter()
                    .position(|x| x == &raw_data_missing.header[*i])
            })
            .map(|(_, x)| x)
            .collect::<Vec<_>>();
        unsafe { &mut *new_data.as_ptr().add(i).cast_mut() }.write(new_r);
    });
    raw_data_missing.header = new_order.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    raw_data_missing.data =
        unsafe { std::mem::transmute::<Vec<MaybeUninit<Vec<String>>>, Vec<Vec<String>>>(new_data) };
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
    let max = inputs.len();
    let atomic = AtomicUsize::new(0);
    let nucleotides = Mutex::new(Vec::with_capacity(max));
    nucleotides
        .lock()
        .unwrap()
        .extend((0..max).map(|_| MaybeUninit::uninit()));
    std::thread::scope(|s| {
        for _ in 0..num_cpus::get() {
            s.spawn(|| loop {
                let j = atomic.fetch_add(1, Ordering::Relaxed);
                if j >= max {
                    break;
                }
                let input = &inputs[j];
                let output = std::process::Command::new(&ctx.args.samtools)
                    .arg("faidx")
                    .arg(&ctx.args.fasta_ref)
                    .arg(input)
                    .output()
                    .unwrap();
                let output = String::from_utf8(output.stdout).unwrap();
                let n = output
                    .lines()
                    .map(|l| {
                        if l.split(' ').next().unwrap().len() > 1 {
                            "N".to_string()
                        } else {
                            l.to_uppercase()
                        }
                    })
                    .collect::<Vec<_>>();
                nucleotides.lock().unwrap()[j].write(n);
            });
        }
    });
    let nucleotides = unsafe {
        std::mem::transmute::<Vec<MaybeUninit<Vec<String>>>, Vec<Vec<String>>>(
            nucleotides.into_inner().unwrap(),
        )
    };
    let ref_ = raw_data_merged.idx("ref");
    let alt = raw_data_merged.idx("alt");
    let effect_size = raw_data_merged.idx("effect_size");
    let eaf = raw_data_merged.idx("EAF");
    for (mut d, n) in raw_data_missing
        .data
        .into_iter()
        .zip(nucleotides.into_iter().flatten())
    {
        if d[alt] == n {
            let (one, two) = d.split_at_mut(alt.max(ref_));
            let min = alt.min(ref_);
            let max = alt.max(ref_);
            std::mem::swap(&mut one[min], &mut two[max]);
            let es = d[effect_size].parse::<f64>().unwrap();
            d[effect_size] = (-es).to_string();
            let e = d[eaf].parse::<f64>().unwrap();
            d[eaf] = (1.0 - e).to_string();
        }
        raw_data_merged.data.push(d);
    }
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
    let mut raw_data = preformat(&ctx);
    info!("Starting liftover");
    liftover(&ctx, &raw_data);
    info!("Starting dbSNP matching");
    let (raw_data_merged, raw_data_missing) = dbsnp_matching(&ctx, &mut raw_data);
    info!("Starting ref/alt check");
    let final_data = ref_alt_check(&ctx, raw_data_merged, raw_data_missing);
    info!("Writing final data to {}", ctx.args.output_file);
    let file = std::fs::File::create(&ctx.args.output_file).unwrap();
    let mut writer = flate2::write::GzEncoder::new(&file, flate2::Compression::default());
    debug!(len = final_data.data.len(), "Writing rows",);
    writeln!(writer, "{}", final_data.header.join("\t")).unwrap();
    for r in final_data.data {
        writeln!(writer, "{}", r.join("\t")).unwrap();
    }
    writer.finish().unwrap();
    info!("Pipeline complete");
}
