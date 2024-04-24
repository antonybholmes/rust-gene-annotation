use std::{
    cmp,
    collections::{BTreeMap, BTreeSet, HashMap},
};

use crate::loctogene::{GenesResult, GenomicFeature, Level, LoctogeneDb, TSSRegion};
use csv::{Writer, WriterBuilder};
use dna::Location;
use serde::Serialize;

pub const NA: &str = "n/a";
pub const PROMOTER: &str = "promoter";
pub const EXONIC: &str = "exonic";
pub const INTRONIC: &str = "intronic";
pub const INTERGENIC: &str = "intergenic";

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

#[derive(Serialize, Clone)]
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
    is_promoter: bool,
    is_intronic: bool,
    is_exon: bool,
    abs_d: i32,
    d: i32,
}

pub struct Annotate {
    genesdb: LoctogeneDb,
    tss_region: TSSRegion,
    n: u16,
}

impl Annotate {
    pub fn new(genesdb: LoctogeneDb, tss_region: TSSRegion, n: u16) -> Self {
        return Annotate {
            genesdb,
            tss_region,
            n,
        };
    }

    pub async fn annotate(&self, location: &Location) -> GenesResult<GeneAnnotation> {
        let mid: u32 = location.mid();

        // extend search area to account  for promoter
        // let search_loc: Location = Location::new(
        //     &location.chr,
        //     location.start - self.tss_region.offset_5p.abs(),
        //     location.end + self.tss_region.offset_5p.abs(),
        // )?;

        let genes_within: Vec<GenomicFeature> = self
            .genesdb
            .get_genes_within_promoter(
                &location,
                &Level::Transcript,
                cmp::max(self.tss_region.offset_5p(), self.tss_region.offset_3p()),
            )
            .await?;

        // we need the unique ids to symbols
        let mut id_map: HashMap<&str, &str> = HashMap::new();
        let mut promoter_map: HashMap<&str, GeneProm> = HashMap::new();
        //let mut dist_map: HashMap<&str, bool> = HashMap::new();

        for gene in genes_within.iter() {
            let id: &str = &gene.gene_id;

            // println!(
            //     "{} {} {} {} {}",
            //     gene.gene_id, gene.gene_symbol, gene.start, gene.end, gene.strand
            // );

            id_map.insert(id, &gene.gene_symbol);

            //let labels = self.classify_location(location, gene);

            let exons: Vec<GenomicFeature> = self.genesdb.in_exon(&location, &id).await?;

            let is_exon: bool = exons.len() > 0;

            let is_promoter: bool = (gene.strand == "+"
                && mid >= gene.start - self.tss_region.offset_5p()
                && mid <= gene.start + self.tss_region.offset_3p())
                || (gene.strand == "-"
                    && mid >= gene.end - self.tss_region.offset_3p()
                    && mid <= gene.end + self.tss_region.offset_5p());

            let is_intronic = mid >= gene.start && mid <= gene.end;

            let d: i32 = if gene.strand == "+" {
                (gene.start as i32) - (mid as i32)
            } else {
                (gene.end as i32) - (mid as i32)
            };

            //println!("{} {} {}", gene.end - mid, gene.end, mid);

            // update by inserting default case and then updating
            promoter_map
                .entry(id)
                .and_modify(|v: &mut GeneProm| {
                    v.is_intronic = v.is_intronic || is_intronic;
                    v.is_promoter = v.is_promoter || is_promoter;
                    v.is_exon = v.is_exon || exons.len() > 0;

                    let abs_d: i32 = d.abs();

                    if abs_d < v.abs_d {
                        v.d = d;
                        v.abs_d = abs_d;
                    }
                })
                .or_insert(GeneProm {
                    is_promoter,
                    is_intronic,
                    is_exon,
                    d,
                    abs_d: d.abs(),
                });
        }

        // sort the ids by distance
        let mut dist_map: BTreeMap<i32, BTreeSet<&str>> = BTreeMap::new();

        for id in id_map.keys() {
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
                ids.push(id);
            }
        }

        println!("within d {}", ids.len());

        // make a list of the symbols in distance order
        let mut gene_symbols: Vec<&str> = ids.iter().map(|id| id_map[id]).collect::<Vec<&str>>();

        let prom_labels: Vec<String> = ids
            .iter()
            .map(|id| {
                let p = &promoter_map[id];
                make_label(p.is_promoter, p.is_exon, p.is_intronic)
            })
            .collect::<Vec<String>>();

        let mut tss_dists: Vec<String> = ids
            .iter()
            .map(|id| promoter_map[id].d.to_string())
            .collect::<Vec<String>>();

        if ids.len() == 0 {
            ids.push(NA);
            gene_symbols.push(NA);
            tss_dists.push(NA.to_string());
        }

        println!("{} geneids", ids.join(";"));
        println!("{}", gene_symbols.join(";"));
        println!("{}", prom_labels.join(";"));
        println!("{}", tss_dists.join(";"));

        let features: Vec<GenomicFeature> = self
            .genesdb
            .get_closest_genes(&location, self.n, Level::Gene)
            .await?;

        let mut closest_genes: Vec<ClosestGene> = Vec::with_capacity(features.len());

        for i in 0..features.len() {
            let feature = &features[i];

            let prom_label = self.classify_location(location,  feature).await;

            let closest = ClosestGene {
                gene_id: feature.gene_id.to_owned(),
                gene_symbol: feature.gene_symbol.to_owned(),
                tss_dist: feature.dist,
                prom_label,
            };

            closest_genes.push(closest);
        }

        let annotation: GeneAnnotation = GeneAnnotation {
            gene_ids: ids.join(";"),
            gene_symbols: gene_symbols.join(";"),
            prom_labels: prom_labels.join(";"),
            tss_dists: tss_dists.join(";"),
            closest_genes,
        };

        Ok(annotation)
    }

    async fn classify_location(&self, location: &Location, feature: &GenomicFeature) -> String {
        let mid: u32 = location.mid();

        let s: u32 = if feature.strand == "+" {
            feature.start - self.tss_region.offset_5p()
        } else {
            feature.start
        };

        let e: u32 = if feature.strand == "-" {
            feature.end + self.tss_region.offset_5p()
        } else {
            feature.end
        };

        if location.start > e || location.end < s {
            return INTERGENIC.to_string();
        }

        let is_promoter: bool = (feature.strand == "+"
            && mid >= s
            && mid <= feature.start + self.tss_region.offset_3p())
            || (feature.strand == "-"
                && mid >= feature.end - self.tss_region.offset_3p()
                && mid <= e);

        let exons: Vec<GenomicFeature> =
            match self.genesdb.in_exon(&location, &feature.gene_id).await {
                Ok(exons) => exons,
                Err(_) => vec![],
            };

        let is_exon = exons.len() > 0;

        let is_intronic = mid >= feature.start && mid <= feature.end;

        return make_label(is_promoter, is_exon, is_intronic);
    }

    pub async fn make_gene_table(
        &self,
        locations: &Vec<Location>,
        closest_n: u16,
        ts: &TSSRegion,
    ) -> GenesResult<String> {
        let mut wtr: Writer<Vec<u8>> = WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);

        let capacity: usize = 6 + closest_n as usize;

        let mut headers: Vec<String> = Vec::with_capacity(capacity);

        headers.push("Location".to_owned());
        headers.push("ID".to_owned());
        headers.push("Gene Symbol".to_owned());
        headers.push(format!(
            "Relative To Gene (prom=-{}/+{}kb)",
            ts.offset_5p() / 1000,
            ts.offset_3p() / 1000
        ));
        headers.push("TSS Distance".to_owned());

        for i in 1..(closest_n + 1) {
            headers.push(format!("#{} Closest ID", i));
            headers.push(format!("#{} Closest Gene Symbols", i));
            headers.push(format!(
                "#{} Relative To Closet Gene (prom=-{}/+{}kb)",
                i,
                ts.offset_5p() / 1000,
                ts.offset_3p() / 1000
            ));
            headers.push(format!("#{} TSS Closest Distance", i));
        }

        wtr.write_record(&headers)?;

        for location in locations {
            let annotation: GeneAnnotation = self.annotate(location).await?;

            let mut row: Vec<String> = Vec::with_capacity(capacity);

            row.push(location.to_string());
            row.push(annotation.gene_ids);
            row.push(annotation.gene_symbols);
            row.push(annotation.prom_labels);
            row.push(annotation.tss_dists);

            for closest_gene in annotation.closest_genes.iter() {
                row.push(closest_gene.gene_id.clone());
                row.push(closest_gene.gene_symbol.clone());
                row.push(closest_gene.prom_label.clone());
                row.push(closest_gene.tss_dist.to_string());
            }

            wtr.write_record(&row)?;
        }

        let inner: Vec<u8> = wtr.into_inner()?;
        let data: String = String::from_utf8(inner)?;

        Ok(data)
    }
}

fn make_label(is_promoter: bool, is_exon: bool, is_intronic: bool) -> String {
    let mut labels: Vec<&str> = Vec::with_capacity(2);

    if is_promoter {
        labels.push(PROMOTER);
    }

    if is_exon {
        labels.push(EXONIC);
    } else {
        if is_intronic {
            labels.push(INTRONIC);
        }
    }

    return labels.join(",");
}
