ROOT := $(CURDIR)
R_LIB := $(ROOT)/r-library
R_ENV := R_LIBS_USER="$(R_LIB)$(if $(R_LIBS_USER),:$(R_LIBS_USER)):$(HOME)/R/library"
PKG_VERSION := $(shell Rscript --vanilla -e 'cat(read.dcf("software/prime/DESCRIPTION")[1,"Version"])')
VERIFIED := python3 experiments/99_validate/verified_workflow.py
FIGURE_PYTHON ?= python3

.PHONY: deps bootstrap verify-workflow install test smoke check manual validate run-main run-concrete run-observed run-submission refresh-benchmark expand-competitors validate-main validate-submission validate-release figure-concrete artifacts provenance

deps:
	mkdir -p "$(R_LIB)"
	$(R_ENV) Rscript --vanilla -e 'install.packages(c("FNN", "MASS", "Rcpp", "RcppArmadillo", "testthat", "twinning", "ggplot2", "jsonlite"), lib="$(R_LIB)", repos="https://cloud.r-project.org")'

bootstrap:
	$(R_ENV) $(VERIFIED) bootstrap

verify-workflow:
	python3 -m unittest discover -s experiments/99_validate -p 'test_verified_workflow.py' -v

install:
	mkdir -p "$(R_LIB)"
	$(R_ENV) Rscript --vanilla -e 'Rcpp::compileAttributes("software/prime")'
	$(R_ENV) R CMD INSTALL --preclean -l "$(R_LIB)" software/prime

test: install
	$(R_ENV) Rscript --vanilla -e 'testthat::test_local("software/prime", reporter="summary", load_package="installed")'

smoke: install
	$(R_ENV) Rscript --vanilla experiments/99_validate/01_smoke.R

validate:
	Rscript --vanilla experiments/99_validate/00_repository_snapshot.R

check:
	R CMD build software/prime
	$(R_ENV) R CMD check --no-manual prime_$(PKG_VERSION).tar.gz

manual:
	mkdir -p docs tmp/rd2pdf
	R CMD Rd2pdf --no-preview --force --output=docs/prime.pdf --build-dir=tmp/rd2pdf software/prime

run-main:
	PRIME_SELECTION_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/01_selection/01_run.R
	PRIME_RULE_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/07_selection_calibration/01_buffer_ablation.R
	PRIME_LEVEL_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/07_selection_calibration/02_clipped_level.R
	PRIME_SELECTION_LINK_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/02_selection_to_pme/01_run.R
	PRIME_FRIEDMAN_DETAIL_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/02_selection_to_pme/04_friedman_detail.R
	PRIME_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/03_pme/01_run.R
	PRIME_CORRELATED_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/04_correlated_pipelines/01_run.R
	PRIME_SCALING_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/05_scaling/01_run.R
	PRIME_TWINNING_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/06_twinning/01_run.R
	PRIME_OBSERVED_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/09_observed_data/01_run.R
	PRIME_CONCRETE_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/10_concrete_application/01_run.R
	PRIME_OBSERVED_APPLICATION_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/11_observed_application/01_run.R

run-concrete:
	PRIME_CONCRETE_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/10_concrete_application/01_run.R

run-observed:
	PRIME_OBSERVED_APPLICATION_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace main --script experiments/11_observed_application/01_run.R

run-submission:
	PRIME_CORES=4 $(R_ENV) Rscript --vanilla experiments/12_submission_update/01_run_matched_k.R

refresh-benchmark:
	PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace refresh --script experiments/01_selection/03_refresh_ackley3.R
	PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace refresh --script experiments/02_selection_to_pme/03_refresh_ackley3.R
	PRIME_TWINNING_PROFILE=main PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace refresh --script experiments/06_twinning/01_run.R

expand-competitors:
	PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace expand --script experiments/04_correlated_pipelines/03_expand_competitors.R
	PRIME_CORES=6 $(R_ENV) $(VERIFIED) run --namespace expand --script experiments/05_scaling/03_expand_competitors.R

validate-main:
	@test -f manuscript/main.tex || { printf "%s\n" "NOT_IN_PUBLIC_EXPORT: validate-main requires the private historical manuscript; use make validate-release for supported saved code/data checks."; exit 2; }
	PRIME_SELECTION_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/01_selection/02_validate.R
	PRIME_SELECTION_CALIBRATION_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/07_selection_calibration/03_validate.R
	PRIME_SELECTION_LINK_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/02_selection_to_pme/02_validate.R
	PRIME_FRIEDMAN_DETAIL_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/02_selection_to_pme/05_validate_friedman_detail.R
	PRIME_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/03_pme/02_validate.R
	PRIME_CORRELATED_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/04_correlated_pipelines/02_validate.R
	PRIME_SCALING_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/05_scaling/02_validate.R
	PRIME_TWINNING_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/06_twinning/02_validate.R
	PRIME_OBSERVED_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/09_observed_data/02_validate.R
	PRIME_CONCRETE_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/10_concrete_application/02_validate.R
	PRIME_OBSERVED_APPLICATION_PROFILE=main $(R_ENV) $(VERIFIED) run --namespace main --script experiments/11_observed_application/02_validate.R
	$(R_ENV) $(VERIFIED) run --namespace main --script experiments/99_validate/03_pme_reduction_and_knn_parity.R
	$(R_ENV) $(VERIFIED) run --namespace main --script experiments/99_validate/04_manuscript_claims.R

artifacts: artifacts-current-20261011

provenance:
	$(R_ENV) Rscript --vanilla experiments/99_validate/05_write_provenance.R

validate-submission:
	python3 experiments/14_current_numerical/validate.py --historical-data-only

validate-release: validate-submission
	python3 experiments/15_current_numerical_20261011/validate.py

figure-concrete:
	mkdir -p build/research_graphics
	$(FIGURE_PYTHON) experiments/13_concrete_submission/build_concrete_p8_figure.py --output build/research_graphics/concrete_p8_distinct_subsamples.pdf --receipt build/research_graphics/CONCRETE_GRAPHIC_RECEIPT.json

.PHONY: deps-current replay-current
deps-current: deps
	$(R_ENV) Rscript --vanilla -e 'install.packages(c("digest", "sensitivity", "RANN", "gtools"), lib="$(R_LIB)", repos="https://cloud.r-project.org")'

replay-current:
	$(R_ENV) python3 experiments/14_current_numerical/run_reproduction.py smoke

.PHONY: validate-current-20261011 artifacts-current-20261011 deps-current-20261011
validate-current-20261011:
	python3 experiments/15_current_numerical_20261011/validate.py

artifacts-current-20261011:
	$(R_ENV) python3 experiments/15_current_numerical_20261011/build.py

deps-current-20261011:
	mkdir -p "$(R_LIB)"
	$(R_ENV) Rscript --vanilla -e 'install.packages(c("ggplot2", "patchwork", "jsonlite"), lib="$(R_LIB)", repos="https://cloud.r-project.org")'
