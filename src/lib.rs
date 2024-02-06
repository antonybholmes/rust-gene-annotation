use std::{
    cmp, collections::{BTreeMap, HashMap}, error::Error, num
};

use dna::Location;
use loctogene::TSSRegion;
use serde::Serialize;

mod tests;

const NA: &str = "n/a";
const PROMOTER: &str = "promoter";
const EXONIC: &str = "exonic";
const INTRONIC: &str = "intronic";

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

#[derive(Serialize)]
pub struct GeneAnnotation {
    pub gene_ids: Vec<String>,
    pub gene_symbols: Vec<String>,
    pub relative_to_gene: Vec<String>,
    pub tss_distance: Vec<i32>,
}

struct GeneProm {
    is_exon: bool,
    is_promoter: bool,
    abs_d: i32,
    d: i32,
}

pub struct Annotate {
    genesdb: loctogene::Loctogene,
    tss_region: TSSRegion,
    n: u16,
}

impl Annotate {
    pub fn new(genesdb: loctogene::Loctogene, tss_region: TSSRegion, n: u16) -> Annotate {
        return Annotate {
            genesdb,
            tss_region,
            n,
        };
    }

    pub fn annotate(&self, location: &Location) -> Result<GeneAnnotation, Box<dyn Error>> {
        let mid: i32 = location.mid();

        let search_loc: Location = Location::new(
            &location.chr,
            location.start + self.tss_region.offset_5p,
            location.end - self.tss_region.offset_5p,
        )?;

        let genes_within: Vec<loctogene::GenomicFeature> = self
            .genesdb
            .get_genes_within(&search_loc, loctogene::Level::Transcript)?;

        // we need the unique ids to symbols
        let mut id_map: HashMap<&str, &str> = HashMap::new();
        let mut is_promoter_map: HashMap<&str, GeneProm> = HashMap::new();
        //let mut dist_map: HashMap<&str, bool> = HashMap::new();

        for (i, gene) in genes_within.iter().enumerate() {
            let s: &str = &gene.gene_symbol;

            println!("{} {}", gene.gene_id, gene.gene_symbol);
            id_map.insert(s, &gene.gene_id);

            let is_promoter = (gene.strand == "+"
                && mid >= gene.start + self.tss_region.offset_5p
                && mid <= gene.start - self.tss_region.offset_3p)
                || (gene.strand == "-"
                    && mid >= gene.end + self.tss_region.offset_3p
                    && mid <= gene.end - self.tss_region.offset_5p);

            let d: i32 = (gene.start - mid);

            // update using or so if we become true, we stay true
            is_promoter_map
                .entry(s)
                .and_modify(|v| {
                    v.is_promoter = v.is_promoter || is_promoter;
                    v.d = cmp::min(v.d, d);
                    v.abs_d = cmp::min(v.abs_d, d.abs())
                })
                .or_insert(GeneProm {
                    is_exon: false,
                    is_promoter: false,
                    d: i32::MAX,
                    abs_d: i32::MAX,
                });
        }

        let gene_symbols: Vec<&str> = id_map.keys().map(|k| *k).collect();

        print!("{}", gene_symbols.join(";"));

        let gene_ids: Vec<String> = if genes_within.len() > 0 {
            genes_within
                .iter()
                .map(|g| g.gene_id.clone())
                .collect::<Vec<String>>()
        } else {
            vec![]
        };

        let gene_symbols: Vec<String> = if genes_within.len() > 0 {
            genes_within
                .iter()
                .map(|g| g.gene_symbol.clone())
                .collect::<Vec<String>>()
        } else {
            vec![]
        };

        // need to know if location is in an exon
        let exons: Vec<loctogene::GenomicFeature> =
            self.genesdb
                .get_closest_genes(&location, self.n, loctogene::Level::Exon)?;

        let mut relative_to_gene: Vec<String> = Vec::with_capacity(genes_within.len());

        for gene_id in genes_within.iter().map(|g| g.gene_id.clone()) {
            let mut labels: Vec<String> = Vec::with_capacity(2);

            let promoters: Vec<loctogene::GenomicFeature> =
                self.genesdb
                    .in_promoter(&location, &gene_id, &self.tss_region)?;

            if promoters.len() > 0 {
                labels.push(PROMOTER.to_owned())
            }

            // see if this gene id is in an exon
            let exons: Vec<loctogene::GenomicFeature> =
                match self.genesdb.in_exon(&location, &gene_id) {
                    Ok(exons) => exons,
                    Err(err) => return Err(err),
                };

            if exons.len() > 0 {
                labels.push(EXONIC.to_owned())
            } else {
                labels.push(INTRONIC.to_owned())
            }

            let label: String = labels.join(",");

            relative_to_gene.push(label);
        }

        let tss_distance: Vec<i32> = Vec::with_capacity(genes_within.len());

        let closest_genes: Vec<loctogene::GenomicFeature> =
            self.genesdb
                .get_closest_genes(&location, self.n, loctogene::Level::Gene)?;

        Ok(GeneAnnotation {
            gene_ids,
            gene_symbols,
            relative_to_gene,
            tss_distance,
        })
    }
}
