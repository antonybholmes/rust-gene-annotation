use std::{
    cmp,
    collections::{BTreeMap, BTreeSet, HashMap},
    error::Error,
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
pub struct ClosestGene {
    pub gene_id: String,
    pub gene_symbol: String,
    pub prom_label: String,
    pub tss_dist: i32,
}

#[derive(Serialize)]
pub struct GeneAnnotation {
    pub gene_ids: String,
    pub gene_symbols: String,
    pub prom_labels: String,
    pub tss_dists: String,
    pub closest_genes: Vec<ClosestGene>,
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

        // extend search area to account  for promoter
        // let search_loc: Location = Location::new(
        //     &location.chr,
        //     location.start - self.tss_region.offset_5p.abs(),
        //     location.end + self.tss_region.offset_5p.abs(),
        // )?;

        println!("{}", self.tss_region);
        println!("{}", location);

        let genes_within: Vec<loctogene::GenomicFeature> = self.genesdb.get_genes_within_promoter(
            &location,
            loctogene::Level::Transcript,
            cmp::max(
                self.tss_region.offset_5p.abs(),
                self.tss_region.offset_3p.abs(),
            ),
        )?;

        // we need the unique ids to symbols
        let mut id_map: HashMap<&str, &str> = HashMap::new();
        let mut promoter_map: HashMap<&str, GeneProm> = HashMap::new();
        //let mut dist_map: HashMap<&str, bool> = HashMap::new();

        for gene in genes_within.iter() {
            let id: &str = &gene.gene_id;

            println!(
                "{} {} {} {} {}",
                gene.gene_id, gene.gene_symbol, gene.start, gene.end, gene.strand
            );

            id_map.insert(id, &gene.gene_symbol);

            let exons: Vec<loctogene::GenomicFeature> = self.genesdb.in_exon(&location, id)?;

            let is_promoter: bool = (gene.strand == "+"
                && mid >= gene.start + self.tss_region.offset_5p
                && mid <= gene.start + self.tss_region.offset_3p)
                || (gene.strand == "-"
                    && mid >= gene.end - self.tss_region.offset_3p
                    && mid <= gene.end - self.tss_region.offset_5p);

            let d: i32 = if gene.strand == "+" {
                gene.start - mid
            } else {
                gene.end - mid
            };

            println!("{} {} {}", gene.end - mid, gene.end, mid);

            // update by inserting default case and then updating
            promoter_map
                .entry(id)
                .and_modify(|v| {
                    v.is_promoter = v.is_promoter || is_promoter;
                    v.is_exon = v.is_exon || exons.len() > 0;

                    let abs_d: i32 = d.abs();

                    if abs_d < v.abs_d {
                        v.d = d;
                        v.abs_d = abs_d;
                    }
                })
                .or_insert(GeneProm {
                    is_exon: false,
                    is_promoter: false,
                    d: i32::MAX,
                    abs_d: i32::MAX,
                });
        }

        // sort the ids by distance
        let mut dist_map: BTreeMap<i32, BTreeSet<&str>> = BTreeMap::new();

        for id in id_map.keys().map(|k| *k) {
            dist_map
                .entry(promoter_map.get(id).unwrap().abs_d)
                .and_modify(|v| {
                    v.insert(id);
                })
                .or_insert({
                    let mut dids: BTreeSet<&str> = BTreeSet::new();
                    dids.insert(id);
                    dids
                });
        }

        // now put the ids in distance order
        let mut ids: Vec<&str> = Vec::with_capacity(id_map.len());

        for d in dist_map.keys() {
            for id in dist_map.get(d).unwrap() {
                ids.push(*id);
            }
        }

        let gene_symbols: String = ids
            .iter()
            .map(|id: &&str| *id_map.get(id).unwrap())
            .collect::<Vec<&str>>()
            .join(",");

        let labels = ids
            .iter()
            .map(|id: &&str| {
                let p = promoter_map.get(id).unwrap();
                let mut labels: Vec<&str> = Vec::with_capacity(2);

                if p.is_promoter {
                    labels.push(PROMOTER);
                }

                if p.is_exon {
                    labels.push(EXONIC);
                } else {
                    labels.push(INTRONIC);
                }

                labels.join(",")
            })
            .collect::<Vec<String>>()
            .join(";");

        let prom_labels: String = ids
            .iter()
            .map(|id: &&str| {
                let p = promoter_map.get(id).unwrap();
                let mut labels: Vec<&str> = Vec::with_capacity(2);

                if p.is_promoter {
                    labels.push(PROMOTER);
                }

                if p.is_exon {
                    labels.push(EXONIC);
                } else {
                    labels.push(INTRONIC);
                }

                labels.join(",")
            })
            .collect::<Vec<String>>()
            .join(";");

        let tss_dists: String = ids
            .iter()
            .map(|id: &&str| promoter_map.get(id).unwrap().d.to_string())
            .collect::<Vec<String>>()
            .join(";");

        println!("{}", ids.join(";"));
        println!("{}", gene_symbols);
        println!("{}", labels);

        println!(
            "{}",
            ids.iter()
                .map(|id| *id_map.get(id).unwrap())
                .collect::<Vec<&str>>()
                .join(";")
        );
        println!(
            "{}",
            dist_map
                .keys()
                .map(|k| k.to_string())
                .collect::<Vec<String>>()
                .join(";")
        );
        //print!("{}", gene_symbols.iter().collect::<Vec<&str>>().sort());

        let closest_genes: Vec<loctogene::GenomicFeature> =
            self.genesdb
                .get_closest_genes(&location, self.n, loctogene::Level::Gene)?;

        let annotation: GeneAnnotation = GeneAnnotation {
            gene_ids: ids.join(";"),
            gene_symbols,
            prom_labels,
            tss_dists,
            closest_genes: closest_genes
                .iter()
                .map(|cg| {
                    ClosestGene {
                    gene_id: cg.gene_id.to_owned(),
                    gene_symbol: cg.gene_symbol.to_owned(),
                    tss_dist: cg.dist,
                    prom_label: "".to_owned(),
                }})
                .collect(),
        };

        Ok(annotation)
    }
}
