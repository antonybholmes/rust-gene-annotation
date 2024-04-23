

#[cfg(test)]
use std::error::Error;
#[cfg(test)]
use dna::Location;
 
#[cfg(test)]
use crate::annotate::Annotate;
#[cfg(test)]
use crate::annotate::GeneAnnotation;
#[cfg(test)]
use crate::loctogene::GenomicFeature;
#[cfg(test)]
use crate::loctogene::Level;
#[cfg(test)]
use crate::loctogene::LoctogeneDb;
#[cfg(test)]

#[cfg(test)]
use crate::loctogene::TSSRegion;

#[test]
fn test_annotation() ->Result<(), Box<dyn Error>>{
    

    //let loc: Location = Location::parse("chr3:187721370-187733550")?;

    


    let loc: Location = Location::parse("chr3:187745448-187745468")?;

    let genesdb: LoctogeneDb = LoctogeneDb::new("../docker-rust-edb-api/data/loctogene/grch38.db")?;


    let annotatedb: Annotate = Annotate::new(genesdb, TSSRegion::default(), 10);

    let annotation: GeneAnnotation = annotatedb.annotate(&loc)?;

    //let js: serde_json::Value = json!(records);

    //println!("{}", js);

    Ok(())
}

#[test]
fn test_within() {
 
    let loc: Location = match Location::parse("chr3:187721370-187733550") {
        Ok(loc)=>loc,
        Err(err)=>panic!("{}", err)
    };

    let genesdb: LoctogeneDb = match LoctogeneDb::new("../docker-rust-edb-api/data/loctogene/grch38.db") {
        Ok(db)=>db,
        Err(err)=>panic!("{}", err)
    };

    let records:Vec<GenomicFeature>  =  match genesdb.get_genes_within(&loc, &Level::Gene) {
        Ok(records)=>records,
        Err(err)=>panic!("{}", err)
    };

    let js: serde_json::Value = json!(records);

    println!("{}", js);

}

#[test]
fn test_closest() {
    let loc: Location = match Location::parse("chr3:187721370-187733550") {
        Ok(loc)=>loc,
        Err(err)=>panic!("{}", err)
    };

    let genesdb: LoctogeneDb = match LoctogeneDb::new("../docker-rust-edb-api/data/loctogene/grch38.db") {
        Ok(db)=>db,
        Err(err)=>panic!("{}", err)
    };

    let records:Vec<GenomicFeature>  =  match genesdb.get_closest_genes(&loc, 10, Level::Gene) {
        Ok(records)=>records,
        Err(err)=>panic!("{}", err)
    };

    let js: serde_json::Value = json!(records);

    println!("{}", js);

}