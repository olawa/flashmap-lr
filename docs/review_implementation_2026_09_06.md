# Implementering efter review 2026-09-06

Korrekthetsfixarna genomfördes före optimeringarna. Ändringarna ligger i arbetskopiorna för rs-lra och flashmap-standalone; tidigare ändringar i ksw2rs_extd2 är bevarade. Inga commits skapade.

## Genomfört

- **Seedkollisioner:** en gles tabell med full canonical code per berörd range disambiguerar referensseeds som delar residualhash och fingerprint. Både FlashMap och rs-lra fortsätter till rätt range. Även kollision mellan ett borttaget seed och ett kvarvarande seed täcks. Detta fixar förlusten av rätt referensträfflista; den vanliga fingerprintrepresentationen är fortfarande ett probabilistiskt filter för queries som inte finns i referensen.
- **Capped singleton:** inline endast för verkligt råantal ett; samplade enstaka träffar lagras out-of-line med capped=true. First/sample/adaptive/isolated är testade. Isolated med noll kvarvarande träffar behåller capped metadata. Metadata sorteras efter full hash, som dess binärsökning kräver.
- **Format v14:** explicit collision-sektion även om tabellen är tom. Primary, hybrid primary/secondary, SJDB och packed cache är uppdaterade; packed cache är v3. Äldre index avvisas och måste byggas om, eftersom singletonfelets förlorade information inte kan rekonstrueras generellt. Gamla filer har inte skrivits över eller byggts om automatiskt.
- **Gemensam fmi-crate:** FlashMap använder nu rs-lra/crates/fmi för cap-policy inklusive isolated, formatversion och validering av collision-poster. Full metadata/POD-konsolidering återstår. FlashMap har därmed en lokal path-dependency till denna crate.
- **BAM:** long-CIGAR skrivs som två placeholder-operationer och CG:B:I. Namn, koordinater och operationslängder valideras innan batchen skriver output. `BamRecordEncoder::encode_batch` returnerar nu io::Result; CLI vidarebefordrar felen. Remapping filtrerar uppercase CG.
- **AS/NM:** beräknas tillsammans över färdig CIGAR enligt vald scoring. Kedjepoäng och placement-rank behålls separat. Exakt snabbväg och normaliserad bandad helread följer samma AS-kontrakt. CLI-regressioner kontrollerar single och dual, inklusive reverse strand.
- **Worker-pool:** global kreditgräns över queued/running/pending batcher; kredit återlämnas först vid ordnad emission. Antalet begränsas till max(16,4*workers). Readerns egen inläsningsbatch tillkommer; detta är en batchgräns, inte ett absolut byte/RSS-tak för godtyckligt stora reads. Pending-peak och gräns finns i profilen. Sinkfel och panic i source/mapper/sink väcker väntare före scoped join.

## Optimeringar som behållits

- Begränsad DP-cache inom en CIGAR-montering. Endast reparationsprober fyller cachen; slutmonteringen återanvänder dem utan att lagra varje nytt engångsgap. Nyckeln har båda sekvensintervallen, band och scoring; flaggor och z-drop är fasta för denna full-DP-väg. Högst 32 poster, högst 4096 CIGAR-operationer per sparat resultat. Scope lånar de oföränderliga sekvenserna och återställer TLS även vid tidig retur/panic.
- Riktad STR-regression bekräftar samma CIGAR med färre DP-anrop mellan reparation och slutlig gapfyllning. Cacheträffar inkluderar återanvändning av DNA5-arbetet för det identiska anropet.
- DP-adaptern skiljer internt på tom input, budget, omöjligt band, z-drop, ofullständig traceback och ogiltig CIGAR. Befintliga publika Option-API:n är bevarade. Felorsaker och cacheträffar exponeras i profileringen; återförsök med ändrat band har inte tagits bort utan kvalitetsunderlag.
- Forward-query lånas som Cow i interna CIGAR-/helread-vägar. Reverse byggs fortfarande till en orienterad kopia. Det publika build_chain_cigar behåller sitt tidigare ägda returkontrakt.
- Indexbyggaren använder iterator för All/First/Spaced i stället för en offsets-Vec per distinkt seed. Iteratorn testas mot tidigare urval. Isolated behåller sin positionsberoende vektor.
- Indexbyggaren reserverar inte längre en hitplats per råträff innan inline/capping. Tillfällig fingerprinttabell används bara inom residualgrupper med flera seeds. Nya formatgränser valideras innan bygge.

## Mätning

Deterministisk end-to-end-jämförelse av 1000 syntetiska reads på cirka 10 kb, 400 kb referens med korta repeats, blandade insertioner/deletioner/substitutioner och båda orienteringarna. Dual-affine, --profile, 1 respektive 4 workers. Tre växelvisa körningar per variant på samma maskin; inkluderar in-memory-indexbygge och SAM-skrivning. Baslinjen är releasebinären efter korrekthetsfixarna och före DP-cache/query-Cow.

| Workers | Korrigerad baslinje, median | Slutlig variant, median | Kvot baslinje/ny |
|---|---:|---:|---:|
| 1 | 1,927 s | 1,880 s | 1,025 |
| 4 | 0,550 s | 0,554 s | 0,993 |

Alla SAM-poster är byte-identiska mellan körningarna, SHA256 `6d0f35ecd96c4c19c443ff4b4c75b384b014cf820a513b991e4b690e688c81ce`. Headerns programrad exkluderas från jämförelsen. 30 gap-DP-anrop återanvänds; cirka 22 500 återstår. Datasetet är inte tillräckligt för att påstå generell hastighetsvinst eller biologisk kvalitetsförbättring. Tidigare mellanvariant gav i princip oförändrad tid. Indexbyggarens end-to-end-vinst/RSS är inte mätt.

Körbara benchmarkskriptet är `scripts/review_benchmark.py`; ange --baseline, --optimized och --out. Rådata: `docs/review_benchmark_2026_09_06.json`.

Ett separat experiment tog bort dubbla skrivningar i extd2:s arbetsbuffert. Crate-benchmarken visade ingen vinst och flera fall blev långsammare. Experimentet återställdes; inga nya ksw2-optimeringar från denna omgång behålls.

## Verifiering

- rs-lra release: 163 bibliotekstester, 28 binärtester, fyra CLI och två indexkontraktstester passerade. fmi: åtta tester passerade. Därefter tillkom ett test som uttryckligen avvisar v13; debugsviten kördes igen.
- FlashMap: 839 bibliotekstester passerar, tre ignorerade. cargo check --all-targets passerar.
- ksw2: 13 tester och en doctest passerade under buffertförsöket; behållen kärna är den tidigare testade implementationen.
- Samtools quickcheck accepterar long-CIGAR-BAM; samtools view återställer alla 65 536 operationer och 65 536 sekvensbaser.
- Små golden v14-fixtures, byggda av FlashMaps publika API, testar rs-lra:s lookup och visit_hits för båda tidigare indexfelen. FlashMap har egna save/load-regressioner för samtliga cap-policyer och kollisioner.
- git diff --check passerar i samtliga tre projekt. Clippy identifierade små nya test-/stilvarningar som rättades; tidigare varningar om bland annat stort AlignerConfig och dubblerat #[test] finns kvar.

## Återstår

Detta är implementationen av de konkreta korrekthetsfynden och första optimeringssteget. Det gör inte rs-lra till en färdig generell long-read-aligner.

- DNA5-kodning mellan olika DP-intervall och reverse-query mellan separata kandidater återanvänds ännu inte generellt. En helgenomskopia för detta är inte införd.
- Adaptiva band, färre retries med ändrat band och ändrade kandidatbudgetar behöver kvalitetsjämförelse innan sökbeteendet ändras.
- Secondary-utdata, empirisk MAPQ-kalibrering och teknikprofiler för HiFi/ONT återstår.
- Storskalig validering av ultralånga reads, repeats och strukturella varianter återstår, liksom mätning av indexbyggarens peak RSS och prefixstorlekar.
- Ingen ny AVX2-runtimevalidering, distributionslösning för path-crates eller splice-aware RNA/assembly/overlap har genomförts.
