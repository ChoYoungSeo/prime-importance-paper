.PHONY: install test lint build

install:
	python -m pip install -e '.[dev]'

test:
	python -m pytest

lint:
	python -m ruff check src tests

build:
	python -m build
