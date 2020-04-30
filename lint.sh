#!/usr/bin/env bash

echo "## FLAKE8"
python3 -m flake8

echo ""
echo ""
echo ""

echo "## PYLINT"
python3 -m pylint **/*.py \
--extension-pkg-whitelist=gemmi,PySide2 \
--disable=bad-continuation \
--disable=missing-module-docstring \
--disable=missing-class-docstring \
--disable=missing-function-docstring \
--disable=not-an-iterable \
--disable=unsubscriptable-object
