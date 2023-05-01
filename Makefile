
.PHONY: install
install:
	pip install .

codestyle:
	pyupgrade --exit-zero-even-if-changed --py37-plus **/*.py
	isort --settings-file pyproject.toml ./
	black --config pyproject.toml ./

.PHONY: lint
lint: codestyle
