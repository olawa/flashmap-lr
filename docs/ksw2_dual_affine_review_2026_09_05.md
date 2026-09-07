Granskning av rs-lra och ksw2rs_extd2, 2026-09-05
================================================

Uppföljning: [implementerade fixar och mätresultat](/Users/olwal516/dev/projects/rs-lra/docs/ksw2_implementation_2026_09_05.md). Nedan beskrivs tillståndet före fixarna.

Granskade revisioner: rs-lra `66eb8ea`, ksw2rs_extd2 `56e688a`. Fokus är dual-affine-kärnan, adaptern, gap assembly, STR-normalisering, endpoint clipping och relevanta prestandavägar. Detta är en riktad kodgranskning, inte en validering av hela alignerns biologiska precision. Inga produktionskällor har ändrats.

Samtliga befintliga körda tester passerar: rs-lra `cargo test --all-targets --quiet` (182 tester), ksw2rs `cargo test --quiet` (8 vanliga tester och en doctest). Körningen sker på macOS/aarch64 och verifierar därför inte x86-backenderna genom exekvering. rs-lra har en varning om dubblerat `#[test]` i main.rs:2582.

De viktigaste fynden
-------------------

**1. [P1] Dual-affine-scoring följer inte med genom efterbearbetningen.**

Referenser: [prepare.rs:233](/Users/olwal516/dev/projects/rs-lra/src/alignment/prepare.rs:233), [assembly.rs:405](/Users/olwal516/dev/projects/rs-lra/src/alignment/assembly.rs:405), [endpoint.rs:253](/Users/olwal516/dev/projects/rs-lra/src/alignment/endpoint.rs:253), [normalize.rs:62](/Users/olwal516/dev/projects/rs-lra/src/alignment/normalize.rs:62).

DP använder `min(6+2*k, 24+k)`, men `score_cigar_ops` använder endast `6+2*k`. Samma begränsning finns i poängen för endpoint clipping, vars gränssnitt inte ens tar emot den andra gapmodellen. Dessa poäng styr faktiska beslut: om interna STR-ankare ska lösas upp och om readändar ska klippas.

För ett gap på 100 bp är rätt dual-affine-kostnad 124, medan efterbearbetningen räknar 206. Skillnaden växer med gaplängden. `merge_score_diff` antar dessutom att sammanslagning alltid sparar gap_open=6. Två 20-bp-gap kostar 44+44 och ett 40-bp-gap kostar 64: besparingen är 24, inte 6. Den befintliga regeln kan därmed avvisa förbättringar som den nya modellen faktiskt föredrar.

Åtgärd: en gemensam `ScoringPolicy::gap_cost(len)`, använd i CIGAR-poäng, clipping och normalisering. Vid merge beräknas `g(d1)+g(d2)-g(d1+d2)` före mismatchjustering. Clipping måste också hålla isär insertion och deletion när en gap-run räknas. Testa gaplängder runt brytpunkten 18 samt 40/100 bp och säkerställ att scoringbesluten använder samma modell. Detta är verifierat i kod och aritmetik; effekten på ett biologiskt benchmark har inte mätts.

**2. [P1] Traceback-workspace bryter mot säkerhetskontraktet för Vec::set_len.**

Referenser: [extd2/mod.rs:103](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/mod.rs:103), motsvarande befintliga mönster i [extz2/mod.rs:174](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extz2/mod.rs:174).

Efter clear/reserve sätts längden till `p_len` innan bytes har initierats. Att aktiva DP-celler skrivs senare uppfyller inte kontraktet. Padding och rader efter z-drop behöver dessutom aldrig skrivas. Workspace har publika prepare-metoder och härledd Debug/Clone, så hela bufferten kan exponeras även utanför tracebackens avsedda läsmönster.

[Rusts dokumentation för set_len](https://doc.rust-lang.org/std/vec/struct.Vec.html#method.set_len) kräver att de nytillkomna elementen redan är initierade. Detta är ett konstaterat kontraktsbrott, inte en observerad krasch i testkörningen.

Åtgärd: använd initialiserad lagring som referenslösning. För att undvika upprepad nollställning kan bufferten behålla sin initierade längd och endast växa med resize vid behov. Alternativt krävs en genomgående korrekt MaybeUninit-design, inklusive SIMD-skrivningar och begränsning av vilka bytes som får läsas. Att bara flytta set_len till slutet räcker inte om hela den deklarerade regionen inte skrivits.

**3. [P2] Lokal extension placeras felaktigt på sekvensernas suffix.**

Referens: [dp.rs:202](/Users/olwal516/dev/projects/rs-lra/src/dp.rs:202), även äldre single-affine-vägen vid rad 104.

Reproducerat med query `ACGTACGTAAAAAAAA`, target `ACGTACGTCCCCCCCC`, band 32. `align_local_dual_affine` ger score=16, CIGAR=8M, query/ref_start=8, end=16 och NM=8. DP-träffen är i själva verket det perfekta prefixet `[0,8)` med NM=0. Adaptern härleder start som `len-consumed` trots att forward extension är förankrad vid början, och ignorerar DP:s endpointkoordinater.

Åtgärd: returnera koordinater för den extension som faktiskt körts. Om suffix-extension önskas måste input reverseras och CIGAR/koordinater transformeras konsekvent. Trimning av terminala indels kräver också att poäng och endpoints förblir förenliga. Nuvarande mapper använder full/banded-vägarna, så detta är främst ett fel i det publika local-API:t, inte belägg för att standardmappningen placerar alla reads fel.

**4. [P2] mte_q räknas med SIMD-padding i stället för sista verkliga targetposition.**

Referenser: [extd2/core.rs:1619](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/core.rs:1619), [extd2/core.rs:1854](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/core.rs:1854).

`ez.mte_q = r - env` använder den avrundade 16-byte-gränsen. För en perfekt 8-basersmatch blir mte=16 men mte_q=-1; korrekt queryindex är 7. För tre baser fås -11. Använd `r - en0 as i32` i båda implementationerna, precis som extz2 redan gör.

Den lokala C-referensen innehåller samma uttryck i c/ksw2_extd2_sse.c:361. Differentialtesterna kan därför vara gröna trots ett semantiskt fel. rs-lra använder för närvarande inte detta resultatfält, vilket begränsar den omedelbara mapperpåverkan.

**5. [P2] Publika scoringparametrar accepterar gapmodeller som randinitieringen inte hanterar.**

Referenser: [extd2/core.rs:1485](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/core.rs:1485), [extd2/core.rs:1706](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/core.rs:1706).

Reproducerat utan z-drop, query=A, target=AAAA, match=2:

| q/e/q2/e2 | Förväntad global poäng | Faktisk poäng |
| --- | ---: | ---: |
| 6/2/24/1 | -10 | -10 |
| 6/2/24/2 | -10 | -28 |
| 6/1/24/2 | -7 | -9 |

När den ena affina funktionen dominerar den andra är detta matematiskt vanlig affine alignment. Den särskilda long_thres/long_diff-initieringen hanterar inte dessa konfigurationer korrekt. Standardparametrarna i rs-lra träffas inte av exemplen.

Åtgärd: validera och dokumentera tillåtna scoringparametrar eller reducera dominerade modeller till single-affine. Validera även bytearitmetikens tillåtna intervall och matrisdimensioner vid det säkra API:ts gräns. Parametertester måste omfatta lika slopes och omkastade/dominerade modeller.

**6. [P2] N-poängen skiljer sig från adapterns dokumenterade policy.**

Referenser: [dp.rs:520](/Users/olwal516/dev/projects/rs-lra/src/dp.rs:520), [extd2/core.rs:1739](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/core.rs:1739).

Adaptern anger att ambigua jämförelser är neutrala och sätter mat[24]=0. Dual-affine-kärnans snabba scoringväg tolkar dock noll som `-e2`. Reproducerat för ANA mot ANA: score=3 i snabbvägen och 4 med GENERIC_SC, trots samma matris. Dessutom använder rs-lra:s CIGAR-ompoängsättning exakt bytejämförelse, vilket kan belöna N/N som en match.

Åtgärd: välj en uttrycklig policy för N/N och N/ACGT, implementera den konsekvent i kärna och efterbearbetning och testa den. Att bara aktivera GENERIC_SC ändrar också prestandan och matrisens N/ACGT-värden är för närvarande -4; det är därför ingen fullständig lösning på den dokumenterade neutraliteten.

Prestandapotential, i föreslagen ordning
--------------------------------------

1. **Vektorisera exakt H-uppdatering och max-reduktion.** [extd2/core.rs:347](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/core.rs:347) används av Scalar, SSE4.1, AVX2 och NEON. Alla rs-lra-anrop använder exakt uppdatering, inte APPROX_MAX. extz2 har redan arkitekturspecifika implementationer som kan tjäna som mall. Skillnaden är viktig: extd2 kräver signed i8-utvidgning och saknar extz2:s qe-subtraktion. Bevara tie-breaking för max_t, annars kan z-drop och extension-CIGAR ändras. Detta är den tydligaste första kärnoptimeringen; ingen procentuell vinst har mätts.

2. **Minska antalet DP-celler med försiktig bandstrategi.** Adaptern kan omedelbart avstå om längdskillnaden överskrider bandet och full consumption krävs. Försök därefter med ett smalt band och eskalera vid misslyckande eller misstänkt bandkontakt. En lyckad full CIGAR bevisar inte att den globala optimumvägen ligger inom bandet: kompenserande I/D kan ge liten nettodrift men stor lokal drift. Behåll referensfall med långa insertioner/deletioner och tandem repeats. De kvadratiska 16M-vakterna i gapvägen avvisar också vissa billiga bandade problem; en framtida policy bör budgetera verkligt bandarbete och tracebackminne, inte bara längdprodukten. Det är en kvalitets-/kostnadsavvägning, inte en garanterad hastighetsvinst.

3. **Spara kopiering och allokering i adaptern.** [dp.rs:477](/Users/olwal516/dev/projects/rs-lra/src/dp.rs:477) klonar packad CIGAR före konvertering till en ny vektor. Konvertera direkt medan resultatet är lånat och reservera antal operationer. Behåll DNA5-representation per orienterad read/återanvänt referensfönster där samma span provas flera gånger. TLS-workspace och resolve-once-dispatch finns redan; prioritera därför inte att återimplementera dem.

4. **AVX2-namnet betyder ännu inte 32 DP-celler per steg.** [extd2/core.rs:1030](/Users/olwal516/dev/projects/ksw2rs_extd2/src/extd2/core.rs:1030) delegerar DP till SSE4.1. Även substitutionsfyllningens AVX2-wrapper delegerade till SSE4.1 vid den granskade revisionen (korrigerat efter kontroll vid implementeringen). En riktig 256-bitars-DP-kärna är en senare x86-optimering som kräver korrekt carry mellan 128-bitarslanes och separat benchmark. Den hjälper inte denna Mac/NEON-körning.

5. **Mät dubbelarbete i STR-reparation.** prepare.rs bygger split- och continuous-vägar med `append_gap_with_policy(..., None)` och den slutliga assemblyn kan sedan beräkna gap igen. Dessa försök saknas i de vanliga DP-diagnostiktiderna. Instrumentera först och återanvänd sedan godkända gapresultat för oförändrad geometri/scoring. Överväg inte score-only följt av traceback generellt: för accepterade fall kan två pass bli dyrare än ett.

6. **Behåll packad lokal k-mer-lagring där representationen tillåter.** [anchors.rs:375](/Users/olwal516/dev/projects/rs-lra/src/anchors.rs:375) sorterar packade poster men expanderar dem sedan till två u64-vektorer. En iterator över packade offsets kan ta bort kopian och minska slutlig lagring från 16 till 8 byte per post i packable-fallet. Nuvarande API returnerar `&[u64]`, så det kräver en faktisk API-ändring och fallback för opackbara fall. Tidigare fynd om borttappade oparade indexträffar och obegränsad parprodukt är inte återrapporterade: Stage A2 och begränsad sökning finns nu i koden.

Verifiering före optimering
--------------------------

De befintliga extd2-differentialtesterna använder framför allt korta oberoende slumpsekvenser och två fasta scoringkonfigurationer. Score-only-testet beräknar ett skalärt resultat utan att jämföra det; traceback-randomtestet jämför inte skalärbackend alls. Bench-targeten heter extz2 och mäter inte den nya dual-affine-vägen.

Komplettera med en liten oberoende i32-DP-orakel för dual-affine global scoring utan z-drop. Kontrollera CIGAR-consumption, CIGAR-rescoring, koordinatintervall och scalar/SIMD/C separat. Testa 15/16/17 och 31/32/33 baser, band nära dessa gränser, workspace reuse efter olika storlekar/z-drop, långa gap kring 18-bp-brytpunkten, N och reverse CIGAR. Den medföljande C-koden är ett kompatibilitetsorakel, inte en garanti för semantisk korrekthet.

Benchmarka sedan releasebyggen på korrelerade 1/10/100 kb-sekvenser med realistiska fel, band 32/64/128/256, både score-only och traceback samt återanvänd workspace. Rapportera CPU-tid, allokeringar, peak RSS per worker och faktiskt besökta DP-celler. End-to-end behövs samma reads/index/trådantal och kontroll av CIGAR/NM/MAPQ samt variantutfall, inklusive SV/STR. Inget fullskaligt sådant datasetbenchmark har utförts i denna review.

Sonderna gav samma resultat i debug- och releasebygge. Reproducerbara sonder finns i [ksw2_review_probes.rs](/Users/olwal516/dev/projects/rs-lra/docs/ksw2_review_probes.rs). Kör dem som main.rs i en separat temporär Cargo-bin med path-dependencies på båda projekten och samma ksw2rs-patch som rs-lra.
