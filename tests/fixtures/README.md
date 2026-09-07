# FlashMap v14 interoperability fixtures

Generated with the FlashMap public builder, then saved with `save_to_file`.

- `singleton_v14.fmi`: one 40-base A homopolymer, k=19, w=1, cap=1, policy=first. Its 22 raw occurrences must be reported as Sampled, never Complete.
- `collision_v14.fmi`: two 19-base contigs decoded from canonical codes 449396736 and 943947776, k=19, w=1, cap=200. Both hash to residual 636675604 and fingerprint 0. Queries must return contigs 0 and 1 respectively.

The reproduction program in `docs/review_probes_2026_09_06.rs` constructs these data through FlashMap's API. The `.fmi` format is owned jointly through `crates/fmi`; changes must keep the FlashMap roundtrip tests and `tests/index_contract.rs` green. Timestamps in fixture metadata are not part of lookup semantics.
