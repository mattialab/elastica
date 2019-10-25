all: black fmt localcheck flake8

black:
	@black --version
	@find . -maxdepth 3 -name '*.py'\
		| while read -r src; do black "$$src"; done

flake8:
	@flake8 --version
	@find . -maxdepth 3 -name '*.py'\
		| while read -r src; do flake8 "$$src"; done

pylint:
	@pylint --version
	@find . -maxdepth 3 -name '*.py'\
		| while read -r src; do pylint -rn "$$src"; done

# The last regexp in find matches bash and zsh files
fmt:
	@shfmt --version
	@find . -type f -maxdepth 3 -name '*.sh' | while read -r src; do shfmt -w "$$src"; done

localcheck:
	@shellcheck --version
	@find . -type f -maxdepth 3 -name '*.sh' | while read -r src; do shellcheck -s bash -e SC1071 "$$src"; done
