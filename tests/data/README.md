# Test data

`water64_pilot.npz` and `water64_pilot_expected.json` are copies of
`seamm_mbe`'s pilot-frame fixture (64 waters, the prototype's TinkerCliffs
labels). See `seamm_mbe/tests/data/README.md` for their provenance, units and
the one known difference from the prototype (its truncated kB -> atm factor).
Here they drive the MBE step's own path, fed through the Evaluator-shaped
results the step consumes.
