use flashmap::index::{MinimizerIndex, SeedType};
use flashmap::index::twobit_minimizer::{hash_code, seed_key};
use std::{collections::HashMap, io::Read};
fn build(refs: Vec<(String, Vec<u8>)>, cap: usize) -> MinimizerIndex {
    MinimizerIndex::build(refs,19,1,cap,SeedType::Minimizer,9,3,1,1,None).unwrap()
}
fn decode(mut code:u64)->Vec<u8> { let mut seq=vec![0;19]; for b in seq.iter_mut().rev(){*b=b"ACGT"[(code&3) as usize];code>>=2;}seq }
fn check_rs(idx: &MinimizerIndex, sequences: &[Vec<u8>]) {
    use rs_lra::SeedIndex;
    let path="/tmp/rs-lra-review-20260906/probe.fmi";
    idx.save_to_file(path).unwrap();
    let reader=rs_lra::MinimizerIndex::open(path).unwrap();
    for seq in sequences { for seed in reader.query_seeds(seq) { println!("rs-lra lookup={:?} hits={:?}",reader.lookup(&seed),rs_lra::collect_hits(&reader,&seed)); } }
}
fn main(){
    let idx=build(vec![("repeat".into(),vec![b'A';40])],1);
    let hits=idx.packed_index.as_ref().unwrap().lookup(hash_code(0),0);
    println!("CAP: raw=22 stored={} capped={} (expected true)",hits.len(),hits.is_capped());
    check_rs(&idx,&[vec![b'A';19]]);
    let mut seen=HashMap::new();
    for i in 1u64..1_000_000 {
        let code=i<<14;
        if let Some(old)=seen.insert(seed_key(hash_code(code)),code){
            println!("COLLISION codes={old},{code} key={}",seed_key(hash_code(code)));
            let idx=build(vec![("first".into(),decode(old)),("second".into(),decode(code))],200);
            let packed=idx.packed_index.as_ref().unwrap();
            println!("first lookup={:?}",packed.lookup(hash_code(old),old));
            println!("second lookup={:?}",packed.lookup(hash_code(code),code));
            check_rs(&idx,&[decode(old),decode(code)]);
            break;
        }
    }
    use rs_lra::{Alignment,Cigar,CigarOp,ContigId,MappedRead,MappingResult,Strand};
    let cigar=Cigar::new((0..65536).map(|i|if i%2==0 {CigarOp::Match(1)} else {CigarOp::Ins(1)})).unwrap();
    let a=Alignment::new(ContigId(0),0,Strand::Forward,0,cigar,0,60,32768).unwrap();
    let mapped=MappedRead{name:"long".into(),sequence:vec![b'A';65536],qualities:None,tags:None,aux:None,mapping:MappingResult{primary:Some(a),..Default::default()}};
    let encoder=rs_lra::bam::BamRecordEncoder::from_contigs([(ContigId(0),"ref",100000)]);
    let mut compressed=vec![];encoder.encode_batch(&[mapped],&mut compressed);
    let mut raw=vec![];flate2::read::MultiGzDecoder::new(&compressed[..]).read_to_end(&mut raw).unwrap();
    println!("BAM actual_ops=65536 encoded_n_cigar={} (expected placeholder 2 + CG)",u16::from_le_bytes(raw[16..18].try_into().unwrap()));
}
