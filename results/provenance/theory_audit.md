# PRIME result-to-proof map

This map follows the reviewed manuscript and [supplementary proofs](../../manuscript/sections/supplement_proofs.tex).

| Result | Assumptions and scope | Proof |
|---|---|---|
| Shapley allocation | Nonnegative normalized allocation; ANOVA component formula under independent inputs | Main Section 2.2 |
| Selected-set equivalence | Unique minimal mean-sufficient set; its complement in a sufficient retained set is the greatest null coalition | Supplement S1.1 |
| Noise-risk decomposition | Square-integrable response with conditional-mean-zero noise; heteroskedasticity allowed | Supplement S1.2 |
| Removal of affine mean drift | Conditions L1–L3, fixed dimension, interior queries, h tending to zero and K tending to infinity; unpenalized local OLS bias, with numerical-ridge terms treated separately | Supplement S1.4 |
| Fixed-subset risk consistency | Conditions L1–L3, K tending to infinity and K/n tending to zero; includes the numerical ridge and boundary queries | Supplement S1.5 |
| Forward coverage | Stable ridge/BIC seed, positive separated active gains, uniform gain consistency and vanishing buffers | Supplement S1.6 |
| Active retention | Positive active deletion margin, uniform deletion-risk consistency and vanishing buffers on the relevant active-containing sets | Supplement S1.7 |
| Canonical-value continuity | Local at normalized monotone games with a greatest null coalition; uniform for bounded dimension and a common positive lower bound on nonnull values | Supplement S1.8 |
| Selected-set PME stability | Fixed set with positive signal variance and a greatest null coalition; uniform subset-risk consistency; availability tends to one | Supplement S1.9 |
| End-to-end PRIME consistency | Fixed square-integrable Y marginal, bounded active/selected dimensions, signal variance bounded below, forward/retention conditions with fixed margins and explicit uniform selected-subset risk consistency | Supplement S1.10 |

The query scale is empirical, not an established sampling standard error; the buffer is not a hypothesis test. The selection proofs use its vanishing bound under the stated uniform moment conditions.

Allocation consistency concerns exact projection and canonical allocation, conditional on availability. Fixed numerical tolerances are described in the supplementary implementation protocol.

The [prior audit](historical/theory_audit_before_20260916.md) is preserved unchanged as a historical record; the map above describes the current results.
