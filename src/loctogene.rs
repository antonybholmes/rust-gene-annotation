use std::{
    error::Error,
    fmt::{self, Display},
    string::FromUtf8Error,
};

use csv::IntoInnerError;
use dna::Location;

use serde::Serialize;
use sqlx::{FromRow, Pool, Sqlite};

const WITHIN_GENE_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, $1 - stranded_start 
    FROM genes 
    WHERE level = $2 AND chr = $3 AND ((start <= $4 AND end >= $4) OR (start <= $5 AND end >= $5)) 
    ORDER BY start ASC"#;

const WITHIN_GENE_AND_PROMOTER_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, $1 - stranded_start 
    FROM genes 
    WHERE level = $2 AND chr = $3 AND ((start - $4 <= $5 AND end + $4 >= $5) OR (start - $4 <= $6 AND end + $4 >= $5)) 
    ORDER BY start ASC"#;

const IN_EXON_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, ? - $1 
    FROM genes 
    WHERE level=3 AND gene_id=$2 AND chr=$3 AND ((start <= $4 AND end >= $4) OR (start <= $5 AND end >= $5)) 
    ORDER BY start ASC"#;

const IN_PROMOTER_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, ? - stranded_start 
    FROM genes 
    WHERE level=2 AND gene_id=? AND chr=? AND ? >= stranded_start - ? AND ? <= stranded_start + ? 
    ORDER BY start ASC"#;

const CLOSEST_GENE_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, $1 - stranded_start 
	FROM genes
	WHERE level=$2 AND chr=$3
	ORDER BY ABS(stranded_start - $1) 
	LIMIT $4"#;

#[derive(Serialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Plus = 1,
    Neg = 2,
}

impl From<&str> for Strand {
    fn from(level: &str) -> Self {
        match level {
            "-" => Strand::Plus,
            _ => Strand::Plus,
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Strand::Neg => write!(f, "-"),
            _ => write!(f, "+"),
        }
    }
}

#[derive(Serialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Level {
    Gene = 1,
    Transcript = 2,
    Exon = 3,
}

impl From<&str> for Level {
    fn from(level: &str) -> Self {
        match level {
            "transcript" => Level::Transcript,
            "exon" => Level::Exon,
            "2" => Level::Transcript,
            "3" => Level::Exon,
            _ => Level::Gene,
        }
    }
}

impl From<u8> for Level {
    fn from(level: u8) -> Self {
        match level {
            2 => Level::Transcript,
            3 => Level::Exon,
            _ => Level::Gene,
        }
    }
}

impl fmt::Display for Level {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Level::Gene => write!(f, "Gene"),
            Level::Transcript => write!(f, "Transcript"),
            Level::Exon => write!(f, "Exon"),
        }
    }
}

#[derive(Serialize, Debug, PartialEq, Eq, Clone, Copy)]
pub struct TSSRegion {
    offset_5p: u32,
    offset_3p: u32,
}

impl TSSRegion {
    pub fn new(offset_5p: u32, offset_3p: u32) -> Self {
        return TSSRegion {
            offset_5p: offset_5p,
            offset_3p: offset_3p,
        };
    }

    pub fn offset_5p(self) -> u32 {
        return self.offset_5p;
    }

    pub fn offset_3p(self) -> u32 {
        return self.offset_3p;
    }
}

impl Default for TSSRegion {
    fn default() -> Self {
        Self {
            offset_5p: 2000,
            offset_3p: 1000,
        }
    }
}

impl fmt::Display for TSSRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{},{}]", self.offset_5p, self.offset_3p)
    }
}

//
//    offset_5p: 2000,
//    offset_3p: 1000,
//};

#[derive(Serialize, Debug, PartialEq, Eq, Clone, FromRow)]
pub struct GenomicFeature {
    pub id: u32,
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub strand: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub dist: i32,
}

// #[derive(Serialize)]
// pub struct GenomicFeatures {
//     pub level: Level,
//     pub features: Vec<GenomicFeature>,
// }

//const NO_FEATURES: [Features; 0] = [] .to_vec();

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

#[derive(Debug, Clone)]
pub enum GenesError {
    DatabaseError(String),
    FormatError(String),
}

impl Error for GenesError {}

impl Display for GenesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GenesError::DatabaseError(error) => write!(f, "{}", error),
            GenesError::FormatError(error) => write!(f, "{}", error),
        }
    }
}

impl From<csv::Error> for GenesError {
    fn from(e: csv::Error) -> GenesError {
        return GenesError::FormatError(e.to_string());
    }
}

impl From<FromUtf8Error> for GenesError {
    fn from(e: FromUtf8Error) -> GenesError {
        return GenesError::FormatError(e.to_string());
    }
}

impl<W> From<IntoInnerError<W>> for GenesError {
    fn from(e: IntoInnerError<W>) -> GenesError {
        return GenesError::FormatError(e.to_string());
    }
}

impl From<sqlx::Error> for GenesError {
    fn from(e: sqlx::Error) -> GenesError {
        return GenesError::FormatError(e.to_string());
    }
}

pub type GenesResult<T> = Result<T, GenesError>;
pub type FeaturesResult = GenesResult<Vec<GenomicFeature>>;

pub struct LoctogeneDb {
    pool: Pool<Sqlite>,
}

impl LoctogeneDb {
    pub fn new(pool: Pool<Sqlite>) -> Self {
        // let db: Connection = match Connection::open(file) {
        //     Ok(db) => db,
        //     Err(err) => return Err(format!("{}", err)),
        // };

        // let manager: SqliteConnectionManager = SqliteConnectionManager::file(file);

        // let pool: r2d2::Pool<SqliteConnectionManager> = match r2d2::Pool::builder().build(manager) {
        //     Ok(pool) => pool,
        //     Err(_) => return Err(GenesError::DatabaseError(format!("{} not found", file))),
        // };

        Self { pool }
    }

    // pub fn get_genes_within_stranded(
    //     &self,
    //     location: &Location,
    //     strand: Strand,
    //     level: Level,
    // ) -> Result<Vec<GenomicFeature>, Box<dyn Error>> {
    //     let mid: u32 = location.mid();

    //     let pool: r2d2::PooledConnection<SqliteConnectionManager> = self.pool.get()?;

    //     let mut stmt: rusqlite::Statement<'_> = pool.prepare(WITHIN_GENE_STRANDED_SQL)?;

    //     let mapped_rows = stmt.query_map(
    //         rusqlite::params![
    //             mid,
    //             level as u8,
    //             location.chr,
    //             strand.to_string(),
    //             location.start,
    //             location.start,
    //             location.end,
    //             location.end
    //         ],
    //         |row: &rusqlite::Row<'_>| row_to_feature(row),
    //     )?;

    //     let features: Vec<GenomicFeature> = mapped_rows
    //         .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
    //         .collect::<Vec<GenomicFeature>>();

    //     Ok(features)
    // }

    // pub fn conn(&self) -> GenesResult<PooledConnection<SqliteConnectionManager>> {
    //     match self.pool.get() {
    //         Ok(pool) => Ok(pool),
    //         Err(_) => Err(GenesError::DatabaseError(format!("error getting pool"))),
    //     }
    // }

    pub async fn get_genes_within(&self, location: &Location, level: &Level) -> FeaturesResult {
        let mid: u32 = location.mid();

        //let pool = self.conn()?;

        //let mut stmt = stmt(&pool, WITHIN_GENE_SQL)?;

        let features = sqlx::query_as::<_, GenomicFeature>(WITHIN_GENE_SQL)
            .bind(mid)
            .bind(level.to_string())
            .bind(&location.chr)
            .bind(location.start)
            //.bind(location.start)
            .bind(location.end)
            //.bind(location.end)
            .fetch_all(&self.pool)
            .await?;

        // let mapped_rows = match stmt.query_map(
        //     rusqlite::params![
        //         mid,
        //         *level as u8,
        //         location.chr,
        //         location.start,
        //         location.start,
        //         location.end,
        //         location.end
        //     ],
        //     |row| row_to_feature(row),
        // ) {
        //     Ok(mapped_rows) => mapped_rows,
        //     Err(_) => return Err(GenesError::DatabaseError(format!("error getting rows"))),
        // };

        // let features: Vec<GenomicFeature> = mapped_rows
        //     .filter_map(|x| x.ok())
        //     .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    pub async fn get_genes_within_promoter(
        &self,
        location: &Location,
        level: &Level,
        pad: u32,
    ) -> FeaturesResult {
        let mid: u32 = location.mid();

        //let pool = self.conn()?;

        //let mut stmt = stmt(&pool, WITHIN_GENE_AND_PROMOTER_SQL)?;

        let features = sqlx::query_as::<_, GenomicFeature>(WITHIN_GENE_AND_PROMOTER_SQL)
            .bind(mid)
            .bind(level.to_string())
            .bind(&location.chr)
            .bind(pad)
            .bind(location.start)
            //.bind(location.start)
            .bind(location.end)
            //.bind(location.end)
            .fetch_all(&self.pool)
            .await?;

        // let mapped_rows = match stmt.query_map(
        //     rusqlite::params![
        //         mid,
        //         *level as u8,
        //         location.chr,
        //         pad,
        //         location.start,
        //         pad,
        //         location.start,
        //         pad,
        //         location.end,
        //         pad,
        //         location.end
        //     ],
        //     |row| row_to_feature(row),
        // ) {
        //     Ok(mapped_rows) => mapped_rows,
        //     Err(_) => return Err(GenesError::DatabaseError(format!("error getting rows"))),
        // };

        // let features: Vec<GenomicFeature> = mapped_rows
        //     .filter_map(|x| x.ok())
        //     .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    // Returns the exons that a location is in within a particular gene. Useful
    // for determining if a gene is exonic or not.
    pub async fn in_exon(&self, location: &Location, gene_id: &str) -> FeaturesResult {
        let mid: u32 = location.mid();

        //let pool = self.conn()?;

        let features = sqlx::query_as::<_, GenomicFeature>(IN_EXON_SQL)
            .bind(mid)
            .bind(gene_id)
            .bind(&location.chr)
            .bind(location.start)
            //.bind(location.start)
            .bind(location.end)
            //.bind(location.end)
            .fetch_all(&self.pool)
            .await?;

        // let mut stmt = stmt(&pool, IN_EXON_SQL)?;

        // let mapped_rows = match stmt.query_map(
        //     rusqlite::params![
        //         mid,
        //         gene_id,
        //         location.chr,
        //         location.start,
        //         location.start,
        //         location.end,
        //         location.end
        //     ],
        //     |row| row_to_feature(row),
        // ) {
        //     Ok(mapped_rows) => mapped_rows,
        //     Err(_) => return Err(GenesError::DatabaseError(format!("error getting rows"))),
        // };

        // let features: Vec<GenomicFeature> = mapped_rows
        //     .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
        //    .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    // Returns a list of features if location is in tss of specific gene
    pub async fn in_promoter(
        &self,
        location: &Location,
        gene_id: &str,
        tss_region: &TSSRegion,
    ) -> FeaturesResult {
        let mid: u32 = location.mid();

        let features_pos = sqlx::query_as::<_, GenomicFeature>(IN_PROMOTER_SQL)
            .bind(mid)
            .bind(gene_id)
            .bind(&location.chr)
            .bind(tss_region.offset_5p)
            //.bind(location.start)
            .bind(tss_region.offset_3p)
            //.bind(location.end)
            .fetch_all(&self.pool)
            .await?;

        let features_neg = sqlx::query_as::<_, GenomicFeature>(IN_PROMOTER_SQL)
            .bind(mid)
            .bind(gene_id)
            .bind(&location.chr)
            .bind(tss_region.offset_3p)
            //.bind(location.start)
            .bind(tss_region.offset_5p)
            //.bind(location.end)
            .fetch_all(&self.pool)
            .await?;

        // let pool = self.conn()?;

        // let mut stmt1 = stmt(&pool, IN_PROMOTER_SQL)?;

        // let mapped_rows_1 = match stmt1.query_map(
        //     rusqlite::params![
        //         mid,
        //         gene_id,
        //         location.chr,
        //         mid,
        //         tss_region.offset_5p,
        //         mid,
        //         tss_region.offset_3p,
        //     ],
        //     |row| row_to_feature(row),
        // ) {
        //     Ok(mapped_rows) => mapped_rows,
        //     Err(_) => return Err(GenesError::DatabaseError(format!("error getting rows"))),
        // };

        // let features_pos = mapped_rows_1.filter_map(|x| x.ok());

        // let mut stmt2 = stmt(&pool, IN_PROMOTER_SQL)?;

        // // negative strand so flip tss region
        // let mapped_rows_2 = match stmt2.query_map(
        //     rusqlite::params![
        //         mid,
        //         gene_id,
        //         location.chr,
        //         mid,
        //         tss_region.offset_3p,
        //         mid,
        //         tss_region.offset_5p,
        //     ],
        //     |row| row_to_feature(row),
        // ) {
        //     Ok(mapped_rows) => mapped_rows,
        //     Err(_) => return Err(GenesError::DatabaseError(format!("error getting rows"))),
        // };

        // let features_neg =
        //     mapped_rows_2.filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok());

        let features = features_pos
            .into_iter()
            .chain(features_neg.into_iter())
            .collect();

        Ok(features)
    }

    pub async fn get_closest_genes(
        &self,
        location: &dna::Location,
        n: u16,
        level: Level,
    ) -> FeaturesResult {
        let mid: u32 = location.mid();

        let features = sqlx::query_as::<_, GenomicFeature>(CLOSEST_GENE_SQL)
            .bind(mid)
            .bind(&level.to_string())
            .bind(&location.chr)
            .bind(n)
            //.bind(location.start)
            //.bind(location.end)
            //.bind(location.end)
            .fetch_all(&self.pool)
            .await?;

        //let pool = self.conn()?;

        //let mut stmt = stmt(&pool, CLOSEST_GENE_SQL)?;

        // query_map converts rusqlite into a standard iterator
        // let mapped_rows = match stmt.query_map(
        //     rusqlite::params![mid, level as u8, location.chr, mid, n],
        //     |row| row_to_feature(row),
        // ) {
        //     Ok(mapped_rows) => mapped_rows,
        //     Err(_) => return Err(GenesError::DatabaseError(format!("error getting rows"))),
        // };

        // filter map because the query returns an iterator of results
        // and if element is ok, the data is the feature record. Use
        // filter map to keep only the valid records and convert them to
        // actual data by removing the Ok wrapper
        // let features: Vec<GenomicFeature> = mapped_rows
        //     .filter_map(|x| x.ok())
        //     .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    // Returns element
}

// fn stmt<'a>(
//     conn: &'a r2d2::PooledConnection<SqliteConnectionManager>,
//     sql: &'a str,
// ) -> GenesResult<rusqlite::CachedStatement<'a>> {
//     match conn.prepare_cached(sql) {
//         Ok(stmt) => Ok(stmt),
//         Err(_) => {
//             return Err(GenesError::DatabaseError(format!(
//                 "error preparing statement"
//             )))
//         }
//     }
// }

// pub fn unwrap_stmt<'a>(
//     stmt: Result<rusqlite::CachedStatement<'a>, rusqlite::Error>
// ) -> GenesResult<rusqlite::CachedStatement<'a>> {

//     match stmt {
//         Ok(stmt) => Ok(stmt),
//         Err(_) => {
//             Err(GenesError::DatabaseError(format!(
//                 "error preparing statement"
//             )))
//         }
//     }
// }

// fn row_to_feature(row: &rusqlite::Row<'_>) -> Result<GenomicFeature, rusqlite::Error> {
//     let id: u32 = row.get(0)?;
//     let chr: String = row.get(1)?;
//     let start: u32 = row.get(2)?;
//     let end: u32 = row.get(3)?;
//     let strand: String = row.get(4)?;
//     let gene_id: String = row.get(5)?;
//     let gene_symbol: String = row.get(6)?;
//     let dist: i32 = row.get(7)?;

//     Ok(GenomicFeature {
//         id,
//         chr,
//         start,
//         end,
//         strand,
//         gene_id,
//         gene_symbol,
//         dist,
//     })
// }
