# Design: Ultra-Fast Placement & Endpoint Priority Ranking in RS-LRA

## 1. Problemformulering & Motivation

### A. Repetitiva regioner orsakar kandidat-explosion
Långa HiFi-reads (15–25 kb) som spänner över transposoner (LINE-1, Alu), mikrosatelliter eller segmentella duplikationer har ofta unika flanker i ändarna men hundratals träffar i mitten:
- **Nuläget:** Varje intern repeat-träff ackumulerar stödjande segment ($N \times 100$), vilket genererar upp till 20 konkurrerande `CandidateRegion`s.
- **Konsekvens:** RS-LRA tvingas köra ankar-sökning och Minimap-DP chaining på alla 20 kandidater, trots att de två unika ändarna redan entydigt bestämmer var readen hör hemma. Detta sänker dessutom MAPQ i onödan.

### B. Behov av $O(1)$-placering för genomisk binning / komprimering
För streaming, extern binning (sortering av reads i ~100 kb genomiska block) eller snabb filtrering behövs en metod som omedelbart returnerar readens position utan att köra fullständig DP eller CIGAR-generering.

---

## 2. Komponent 1: Ultra-Fast Placement API (`fast_placement`)

### 2.1 Datastrukturer

```rust
/// Snabb, approximativ placering av en read baserad på änd- och flankprobes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastPlacement {
    pub contig: ContigId,
    pub strand: Strand,
    pub ref_start: u64,
    pub ref_end: u64,
    pub confidence: PlacementConfidence,
    pub left_probe_hits: u8,
    pub right_probe_hits: u8,
    pub diagonal_delta: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PlacementConfidence {
    /// Båda ändarna matchar unikt och diagonalt konsistent med read-längden.
    High,
    /// En ände matchar unikt tillsammans med 1–2 glesa interna stöd.
    Medium,
    /// Svagt eller flertydigt stöd; kräver fullständig mappning.
    Ambiguous,
}
```

### 2.2 Algoritm & Geometri

```text
Read (15 000 bp):
[ L1  L2 ] ------------------- [ M1  M2 ] ------------------- [ R1  R2 ]
 0..1000 bp                     7..8 kb                        14..15 kb
```

1. **Probe Selection (6 probes totalt):**
   - $L_1, L_2$: De två mest sällsynta minimizers i readens första 1 000 bp.
   - $R_1, R_2$: De två mest sällsynta minimizers i readens sista 1 000 bp.
   - $M_1, M_2$: Två sällsynta minimizers i mitten ($0.45 \times L .. 0.55 \times L$).
2. **Lookup & Diagonalkoppling:**
   - Slå upp referensträffar i `SeedIndex`.
   - För varje par $(L_i, R_j)$ på samma `(ContigId, Strand)`:
     $$\text{read\_span} = \text{pos}(R_j) - \text{pos}(L_i)$$
     $$\text{ref\_span} = \text{ref}(R_j) - \text{ref}(L_i)$$
     $$\Delta_{\text{diagonal}} = | \text{ref\_span} - \text{read\_span} |$$
3. **Kriterium för `PlacementConfidence::High`:**
   - $\Delta_{\text{diagonal}} \le \text{tolerance}$ (t.ex. $\le 1\,000$ bp).
   - Både $L$ och $R$ har låg frekvens ($\le 10$ träffar i referensen).
   - Inget annat kontig har ett konsistent par.
   - **Tidsåtgång:** $< 5 \ \mu\text{s}$ per read (endast 6 index-lookups och en enkel par-matchning).

---

## 3. Komponent 2: Endpoint Priority Ranking & Repeat Pruning

### 3.1 Hierarkisk Kandidat-Klassificering

I `src/candidates.rs` klassificeras varje kluster utifrån sina änd-probes:

```rust
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub enum EndpointTier {
    /// Båda ändarna har konsistenta unika probes.
    BothEnds = 3,
    /// Endast en ände har konsistenta probes.
    SingleEnd = 2,
    /// Inga änd-probes; endast interna repeat-probes.
    InternalOnly = 1,
    None = 0,
}
```

### 3.2 Dominant Poängsättning

Vi ändrar poängsättningen i `add_cluster` från en liten additiv bonus (+250) till en **dominant tier-viktning**:

```rust
let base_score = (supporting_segments.len() as i32 * 100)
    + (unique_probes.len() as i32 * 10)
    + cluster.len() as i32;

let score = match endpoint_tier {
    // BothEnds: Får en dominant 2 000-poängs bas + 1.5x multiplikator
    EndpointTier::BothEnds => 2000 + (base_score * 3 / 2),

    // SingleEnd: Normal poängsättning + moderat bonus
    EndpointTier::SingleEnd => 300 + base_score,

    // InternalOnly: Halveras och straffas kraftigt på långa reads
    EndpointTier::InternalOnly if read_len >= 2000 => (base_score / 2) - 400,

    EndpointTier::None => base_score,
};
```

### 3.3 Repeat-Bypass Pruning i `Aligner::map_read`

När kandidaterna sorterats i `Aligner::map_read`:

```rust
// Om kandidat #1 är en bekräftad BothEnds med hög konfidens:
if let Some(top) = candidates.first() {
    if top.endpoint_support == EndpointSupport::BothEnds {
        // Rensa bort alla Tier 3 (InternalOnly) kandidater direkt!
        candidates.retain(|c| c.endpoint_support != EndpointSupport::InternalOnly);
        
        // Om toppen är överlägsen: begränsa sökningen till max 2 kandidater
        if candidates.get(1).map_or(true, |c2| top.score >= c2.score + 1000) {
            candidates.truncate(1);
        }
    }
}
```

---

## 4. Integration i Aligner-arkitekturen (Tvåvägs Execution Path)

```
                       Read Inmatning
                             │
            ┌────────────────┴────────────────┐
            ▼                                 ▼
   [Snabb Binning / Pre-sort]       [Standard Aligner Pipeline]
            │                                 │
     fast_placement()                 extract_read_probes()
            │                                 │
   Returnera Contig/Start/End         cluster_probe_hits()
   på <5 µs per read                  (med BothEnds dominant tier)
                                              │
                                      Endpoint-Bypass Prune
                                      (1 kandidat istället för 20)
                                              │
                                      find_anchors() & Chaining
                                              │
                                      KSW2 Gap DP & CIGAR
```

---

## 5. Förväntade Effekter

1. **Prestanda (Genomströmning):**
   - I repetitiva regioner (t.ex. centromerer och LINE-1-rika kromosomarmar) minskar antalet utvärderade kandidatregioner från 8–20 till **1 kandidat**.
   - Detta eliminerar upp till **60–80% av CPU-cyklerna** i ankar-sökning och DP för dessa reads.
2. **MAPQ & Variant Calling Precision:**
   - Förhindrar falskt låg MAPQ (`mapq = 0` eller `mapq < 20`) när en read med unika flanker korsar en intern repeat.
3. **Binning-funktion:**
   - Gör det möjligt att implementera `rs-lra bin` eller snabbsortering direkt.
