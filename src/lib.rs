use dna::Location;
use serde::Serialize;

mod tests;

//const NO_FEATURES: [Features; 0] = [] .to_vec();

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

#[derive(Serialize)]
pub struct TSSRegion {
    pub offset_5p: i32,
    pub offset_3p: i32,
}

impl TSSRegion {
    pub fn new(offset_5p: i32, offset_3p: i32) -> TSSRegion {
        return TSSRegion {
            offset_5p,
            offset_3p,
        };
    }
}

pub const DEFAULT_TSS_REGION: TSSRegion = TSSRegion{offset_5p:-2000, offset_3p:1000};


#[derive(Serialize)]
pub struct GeneAnnotation {
    pub genes_within: Vec<loctogene::GenomicFeature>,
    pub closest_genes: Vec<loctogene::GenomicFeature>
}

pub struct Annotate {
    genesdb: loctogene::Loctogene,
    tss_region: TSSRegion,
    n: u16,
}

impl Annotate {
    pub fn new(genesdb: loctogene::Loctogene, tss_region: TSSRegion, n:u16) -> Annotate {
        return Annotate { genesdb, tss_region, n };
    }

    pub fn annotate(&self, location: Location ) ->Result<GeneAnnotation, String> {
        let genes_within: Vec<loctogene::GenomicFeature> = match self.genesdb.get_genes_within(&location, loctogene::Level::Gene) {
            Ok(genes_within)=>genes_within,
            Err(err)=>return Err(err)
        };

        let closest_genes: Vec<loctogene::GenomicFeature> = match self.genesdb.get_closest_genes(&location, self.n,loctogene::Level::Gene) {
            Ok(genes_within)=>genes_within,
            Err(err)=>return Err(err)
        };

        Ok(GeneAnnotation{genes_within, closest_genes})
    }
}
