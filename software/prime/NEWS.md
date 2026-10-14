# prime 0.5.2

* Replaced the fixed ridge threshold with BIC selection along the ridge-ranked
  prefix path.
* Added paired querywise calibration to forward gains and sequential backward
  deletion losses.
* Preserved raw risk ordering when all forward candidate values lie below the
  zero-signal boundary; calibration still determines whether the leader enters.
* Constructed the final proportional allocation from noise-adjusted subset
  risks before monotone projection.
* Removed the obsolete `alpha` tuning argument.

# prime 0.4.0

* Initial release of the PRIME screening and allocation workflow.
* Added local-linear subset-value estimation, inclusive forward screening,
  dimension-adjusted backward reduction, and Data Twinning support.
* Added structured print, summary, and plot methods for fitted and screening
  objects.
