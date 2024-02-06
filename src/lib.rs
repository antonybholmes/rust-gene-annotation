use std::collections::HashSet;

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

    pub fn annotate(&self, location: &Location) -> Result<GeneAnnotation, String> {
        let genes_within: Vec<loctogene::GenomicFeature> = match self
            .genesdb
            .get_genes_within(&location, loctogene::Level::Gene)
        {
            Ok(genes_within) => genes_within,
            Err(err) => return Err(err),
        };

        let gene_ids: Vec<String> = if genes_within.len() > 0 {
            genes_within
                .iter()
                .map(|g| g.gene_id.clone())
                .collect::<Vec<String>>()
        } else {
            vec![NA.to_string()]
        };

        let gene_symbols: Vec<String> = if genes_within.len() > 0 {
            genes_within
                .iter()
                .map(|g| g.gene_symbol.clone())
                .collect::<Vec<String>>()
        } else {
            vec![NA.to_string()]
        };

        // need to know if location is in an exon
        let exons: Vec<loctogene::GenomicFeature> =
            match self
                .genesdb
                .get_closest_genes(&location, self.n, loctogene::Level::Exon)
            {
                Ok(exons) => exons,
                Err(err) => return Err(err),
            };

        let mut exonic: Vec<String> = Vec::with_capacity(genes_within.len());

        for gene_id in genes_within.iter().map(|g| g.gene_id.clone()) {
            let mut labels: Vec<String> = Vec::with_capacity(2);

            let promoters: Vec<loctogene::GenomicFeature> =
                match self
                    .genesdb
                    .in_promoter(&location, &gene_id, &self.tss_region)
                {
                    Ok(exons) => exons,
                    Err(err) => return Err(err),
                };

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

            exonic.push(label);
        }

        if exonic.len() == 0 {
            exonic.push(NA.to_owned())
        }

        let closest_genes: Vec<loctogene::GenomicFeature> =
            match self
                .genesdb
                .get_closest_genes(&location, self.n, loctogene::Level::Gene)
            {
                Ok(genes_within) => genes_within,
                Err(err) => return Err(err),
            };

        Ok(GeneAnnotation {
            gene_ids,
            gene_symbols,
        })
    }
}
