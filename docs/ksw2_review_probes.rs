use ksw2rs::*;
fn main() {
    let q = b"ACGTACGTAAAAAAAA";
    let t = b"ACGTACGTCCCCCCCC";
    println!("local suffix: {:?}", rs_lra::align_local_dual_affine(q,t,32));
    let mut mat = [-4i8;25]; for i in 0..4 {mat[i*5+i]=2;} mat[24]=0;
    let q = [0u8,1,2,3,0,1,2,3];
    let mut a = Aligner::new();
    for (open,ext,open2,ext2) in [(6,2,24,1),(6,2,24,2),(6,1,24,2)] {
        let i = Extd2Input {query:&q,target:&q,m:5,mat:&mat,q:open,e:ext,q2:open2,e2:ext2,w:32,zdrop:-1,end_bonus:0,flag:0};
        println!("mte {open}/{ext}/{open2}/{ext2}: {:?}",a.align_extd2(&i));
    }
    for (open,ext,open2,ext2) in [(6,2,24,1),(6,2,24,2),(6,1,24,2)] {
        let qq=[0u8;1]; let tt=[0u8;4];
        let i=Extd2Input {query:&qq,target:&tt,m:5,mat:&mat,q:open,e:ext,q2:open2,e2:ext2,w:32,zdrop:-1,end_bonus:0,flag:0};
        let expected=2-std::cmp::min(open as i32+3*ext as i32,open2 as i32+3*ext2 as i32);
        println!("boundary {open}/{ext}/{open2}/{ext2} expected {expected}: {:?}",a.align_extd2(&i));
    }
    let qn=[0u8,4,0];
    let i = Extd2Input {query:&qn,target:&qn,m:5,mat:&mat,q:6,e:2,q2:24,e2:1,w:32,zdrop:-1,end_bonus:0,flag:0};
    println!("N default {:?}",a.align_extd2(&i));
    println!("N generic {:?}",a.align_extd2(&Extd2Input {flag:KSW_EZ_GENERIC_SC,..i}));
}
