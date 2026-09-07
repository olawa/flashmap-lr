Implementering efter ksw2-review, 2026-09-05
==========================================

Korrekthetsändringarna implementerades och testades först. Därefter mättes den korrigerade kärnan som baslinje, innan SIMD-optimeringen infördes. Ändringarna ligger i både rs-lra och syskonprojektet ksw2rs_extd2 och behöver följas åt; rs-lra använder crate:ns nya flagga `KSW_EZ_LITERAL_N`.

Genomförda korrekthetsfixar
--------------------------

- Gemensam `ScoringPolicy::gap_cost` och CIGAR-poäng för single/dual affine. STR-beslut, gapmerge och endpoint clipping använder båda gapfunktionerna. Merge jämför faktisk substitution- och gapkostnad, även för ambigua baser. Clipping skiljer på sammanhängande I- och D-runs. Terminalpolicyns avsiktliga matchpoäng behålls, så dess befintliga kalibrering ändras inte samtidigt med gapmodellen.
- Tracebacklagring i extz2/extd2 är helt initierad. Bufferten behåller sin initierade längd och bara nya bytes nollställs vid tillväxt. Det undviker både det tidigare set_len-kontraktsbrottet och upprepad nollställning av hela matrisen.
- Lokal forward extension returnerar rätt prefixkoordinater och räknar om poängen efter eventuell trimning av terminala indels.
- `mte_q` använder den verkliga endpointen, utan SIMD-padding. Samma korrigering finns i den medföljande C-referensen, med kommentar. Oberoende tester kontrollerar att koordinaten är meningsfull.
- Dominerade och parallella gapmodeller normaliseras till identiska affina delar. Matris, alfabet och konservativa bytearitmetikgränser valideras vid extd2-ingången; felaktiga icke-tomma inputs ger dokumenterad panic i stället för tyst felaktigt resultat. Detta begränsar vilka explicita scoringparametrar API:t accepterar.
- rs-lra väljer uttryckligen neutral scoring för alla jämförelser med ambigua baser, inklusive N/N. Matrisen och scoringbesluten följer samma regel. Crate:ns standardflaggor behåller KSW2:s traditionella specialtolkning av wildcard=0 för kompatibilitet. NM är fortsatt ett separat mått.
- Ett befintligt x86-kompileringsfel i AVX2-scorefyllningen upptäcktes och undanröjdes när den vägen implementerades och kompilerades.

Genomförda optimeringar
----------------------

- Exakt H/max-uppdatering använder NEON, SSE4.1 respektive AVX2 med signed byteutvidgning och bevarad tie-breaking. Den skalära implementationen finns kvar som referens/fallback.
- AVX2-scorefyllningen arbetar nu med 32 bytes åt gången. Den gamla AVX2-wrappern delegerade även detta till SSE4.1; den första reviewrapportens formulering om detta har korrigerats. Själva AVX2-DP-rekurrensen använder fortfarande 128-bitarssteg.
- Adaptern konverterar lånad packad CIGAR direkt och reserverar utrymme för resultatet. Mellanklonen försvinner i alla fyra anropsvägarna.
- Full alignment avvisar omedelbart längdskillnader större än bandet, innan DNA5-kodning och DP.
- Packbara lokala k-mer-kartor behåller en u64 per träff, i stället för att expandera till två u64-vektorer. Iterationen återställer absoluta koordinater utan tillfällig vektor. Stora k/fönster använder oförlustgivande parlagring; 128-träffarsgränsen är bevarad.
- STR-reparationens tidigare dolda gapförsök skickar nu vidare diagnostiken och ingår i DP-räknare/tider. Detta ger en användbar grund för att mäta eventuell memoization av upprepade gap.

Mätresultat
-----------

Benchmark: [extd2_throughput.rs](/Users/olwal516/dev/projects/ksw2rs_extd2/benches/extd2_throughput.rs). Kör med `cargo bench --bench extd2_throughput` i ksw2rs_extd2. Deterministiska korrelerade sekvenser på 1/10/100 kb, en 20-bp-deletion, substitutioner var 211:e bas, band 32/128, z-drop 100, återanvänd uppvärmd workspace. Varje värde är medianen av tre batcher om minst 120 ms, mätt som väggtid per alignment i releasebygge på Apple M4/ARM64. Den lokala Cargo-konfigurationen använder target-cpu=native.

Med traceback blev throughput **1,48–2,06× högre** än den korrigerade baslinjen. En separat upprepning gav **1,58–2,03×**. Score-only gav **1,53–2,48×**, respektive **1,72–2,42×** i upprepningen. Poäng och CIGAR-hash är identiska i samtliga benchmarkfall. [Rådata från båda körningarna](/Users/olwal516/dev/projects/rs-lra/docs/ksw2_throughput_2026_09_05.csv).

Mätningen gäller crate:ns DP-kärna. Ingen procentuell end-to-end-vinst för hela mappern eller förbättring i biologisk variantprecision påstås. Kartlagringens minskning från 16 till 8 bytes per packbar träff följer av representationen; total RSS för en verklig mappningskörning är inte uppmätt.

Verifiering
-----------

- rs-lra: 156 bibliotekstester, 28 binärtester och fyra CLI-tester passerar i release. Debugtester passerade också under arbetet, inklusive det nya CLI-testet.
- Nytt CLI-test använder --dual-affine och två workers och verifierar 40-bp-deletion/25-bp-insertion på både forward och reverse strand, koordinater, CIGAR-gap och NM.
- ksw2rs_extd2: 13 vanliga tester och en doctest passerar på ARM64; både debug och release har körts. x86_64-testerna passerar under Rosetta med `RUSTFLAGS='-C target-cpu=x86-64' cargo test --target x86_64-apple-darwin`. Rosetta exponerar SSE4.1 men inte AVX2, så AVX2 är kompilerad men inte exekverad här.
- Ett oberoende konventionellt i32-DP-orakel verifierar 400 genererade fall: global poäng, CIGAR-ompoängsättning och consumption, scalar/SIMD, omkastade/dominerade/lika gapmodeller, N och reverse CIGAR. Ytterligare tester täcker SIMD-/bandgränser och workspace reuse efter z-drop.
- 2 000 direkta H-uppdateringsfall per tillgänglig SIMD-backend jämför signed värden, koordinater, tie-breaking och hela bufferten med scalar.
- C-differentialtesternas saknade scalar-assertioner är tillagda. `git diff --check` passerar. Clippy hittade två nya onödiga casts som togs bort; redan befintliga varningar om bland annat dubblerat #[test] och AlignerConfig-enumens storlek ligger utanför ändringen.

Större experiment
-----------------

Adaptiv bandbredd med ändrade sökbeslut, lättade storleksbudgetar, återanvändning av DNA5 mellan kandidater, gapcache över reparationsfaser och en full 256-bitars-DP-rekurrens är inte aktiverade. De behöver egen profilering och/eller datasetvalidering innan de kan motiveras: en fullständig CIGAR inom ett smalare band bevisar inte att optimum ligger där, och cachekostnaden kan överstiga vinsten. Den här implementationen ändrar inte dessa sökheuristiker.
